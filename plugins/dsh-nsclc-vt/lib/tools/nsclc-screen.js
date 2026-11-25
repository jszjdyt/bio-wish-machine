/**
 * `nsclc_screen` — run mono and/or pairwise perturbation screens over one or
 * more mutational backgrounds and rank them the way the manuscript does.
 * @module dsh-nsclc-vt/tools/nsclc-screen
 */
import { writeFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { defineTool } from '../host.js';
import { enumerateJobs, runJobs, scoreScreen, rankCombinations, estimateCost, defaultWorkers, combinationKey, pairKey } from '../screen.js';
import { num, signed, table, ADDITIVE_NOTE, DIFFERENCE_NOTE, MODEL_PREDICTION_NOTE } from '../format.js';

const PER_BACKGROUND_SCHEMA = {
    type: 'object',
    additionalProperties: false,
    properties: {
        background: { type: 'string', required: true },
        additive: { type: 'number', description: ADDITIVE_NOTE },
        difference: { type: 'number', required: true },
        proliferation: { type: 'number', required: true },
        death: { type: 'number', required: true },
    },
};

const OUTPUT_SCHEMA = {
    type: 'object',
    additionalProperties: false,
    properties: {
        evidence_kind: { type: 'string', required: true, enum: ['model_prediction'], description: 'Every number here is a simulation output, not a measurement.' },
        evidence_note: { type: 'string', required: true, description: MODEL_PREDICTION_NOTE },
        mode: { type: 'string', required: true, enum: ['mono', 'pair', 'both'] },
        backgrounds: { type: 'array', required: true, items: { type: 'string' } },
        levels: { type: 'array', required: true, items: { type: 'integer' } },
        candidateCount: { type: 'integer', required: true, description: 'Nodes that were perturbed.' },
        analysesRun: { type: 'integer', required: true },
        elapsedSeconds: { type: 'number', required: true },
        rankBy: { type: 'string', required: true },
        baselines: {
            type: 'array',
            required: true,
            items: {
                type: 'object',
                additionalProperties: false,
                properties: {
                    background: { type: 'string', required: true },
                    proliferation: { type: 'number', required: true },
                    death: { type: 'number', required: true },
                    difference: { type: 'number', required: true },
                },
            },
        },
        combinations: {
            type: 'array',
            description: 'Present for pair/both modes, ranked best-first by the requested criterion.',
            items: {
                type: 'object',
                additionalProperties: false,
                properties: {
                    combination: { type: 'string', required: true },
                    synergyMean: { type: 'number', required: true, description: 'Mean additive score across the screened backgrounds; negative = synergistic.' },
                    synergySum: { type: 'number', required: true },
                    improved: { type: 'integer', required: true, description: 'Backgrounds where the pair beat its better single component.' },
                    worse: { type: 'integer', required: true },
                    equal: { type: 'integer', required: true },
                    backgroundsScored: { type: 'integer', required: true },
                    perBackground: { type: 'array', required: true, items: PER_BACKGROUND_SCHEMA },
                },
            },
        },
        singles: {
            type: 'array',
            description: 'Present for mono/both modes, ranked best-first.',
            items: {
                type: 'object',
                additionalProperties: false,
                properties: {
                    node: { type: 'string', required: true },
                    meanDifference: { type: 'number', required: true, description: `Mean across backgrounds of ${DIFFERENCE_NOTE}` },
                    meanProliferation: { type: 'number', required: true },
                    meanDeath: { type: 'number', required: true },
                    deltaVsBaseline: { type: 'number', required: true, description: 'Mean difference change versus the untreated background; negative is better.' },
                    perBackground: { type: 'array', required: true, items: PER_BACKGROUND_SCHEMA },
                },
            },
        },
        csvPath: { type: 'string', description: 'Absolute path of the written result table, when output_csv was set.' },
        notes: { type: 'array', items: { type: 'string' } },
    },
};

/** Escape one CSV field. */
const csvField = (v) => {
    const s = v === undefined || v === null ? '' : String(v);
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
};

/**
 * @param {ReturnType<import('../analysis.js').createAnalyst>} analyst
 * @param {{maxAnalyses?: number, workers?: number}} [options]
 */
export function makeNsclcScreenTool(analyst, options = {}) {
    const model = analyst.model;
    const defaultMaxAnalyses = Number.isFinite(options.maxAnalyses) ? options.maxAnalyses : 60000;
    const workers = Number.isFinite(options.workers) && options.workers > 0 ? options.workers : defaultWorkers();

    return defineTool({
        name: 'nsclc_screen',
        description: `Screen the NSCLC virtual tumor for therapeutic perturbations. In "mono" mode it knocks down every candidate node in every background and ranks the single agents by how far they move the phenotype readouts; in "pair" mode it additionally evaluates every pair and scores each combination with the manuscript's additivity score. ${ADDITIVE_NOTE} ${DIFFERENCE_NOTE} Use candidates to keep a pair screen small.`,
        parameters: {
            backgrounds: { type: 'array', items: { type: 'string' }, description: 'Backgrounds to screen: published cell-line names ("a549", "pc9", "healthy"), explicit genotypes ("RAS:3,STK11:0"), or a mix. Defaults to all seven published backgrounds.' },
            mode: { type: 'string', enum: ['mono', 'pair', 'both'], description: 'Which screens to run (default mono).' },
            candidates: { type: 'array', items: { type: 'string' }, description: 'Restrict perturbations to these nodes/genes. Strongly recommended for pair screens — the full pairwise space is ~13k pairs per background.' },
            exclude: { type: 'array', items: { type: 'string' }, description: 'Nodes to leave unperturbed.' },
            levels: { type: 'array', items: { type: 'integer' }, description: 'Levels each candidate is forced to (default [0], i.e. knockout — the manuscript\'s convention). Use [0,3] to also screen full activation.' },
            rank_by: { type: 'string', enum: ['synergy', 'benefit', 'window', 'death'], description: 'Combination ranking (pair mode). synergy = most negative predicted difference-vs-best-single first; benefit = best mean predicted benefit; window = tumour-selective benefit, mean(tumour benefit) - benefit(healthy), requires the healthy background; death = most predicted cell death. Default synergy in pair mode, benefit in mono mode.' },
            top: { type: 'integer', description: 'Rows to return per ranked list (default 25).' },
            min_backgrounds_improved: { type: 'integer', description: 'Keep only combinations that beat their better single component in at least this many backgrounds. The manuscript uses 2 of 7.' },
            output_csv: { type: 'string', description: 'Write the full per-background result table here (path relative to the working directory, or absolute).' },
            max_analyses: { type: 'integer', description: 'Refuse to start if the enumeration exceeds this many model runs (default 60000).' },
        },
        output: {
            schema: OUTPUT_SCHEMA,
            render: (_args, value) => {
                const lines = [`${value.mode} screen over ${value.backgrounds.length} background(s) [${value.backgrounds.join(', ')}]; `
                    + `${value.analysesRun} model runs in ${value.elapsedSeconds.toFixed(1)}s.`];
                lines.push('');
                lines.push('Untreated baselines:');
                lines.push(table(['background', 'prolif', 'death', 'difference'], value.baselines.map((b) => [
                    b.background, num(b.proliferation), num(b.death), num(b.difference),
                ]), { align: ['l', 'r', 'r', 'r'] }));

                if (Array.isArray(value.combinations) && value.combinations.length > 0) {
                    lines.push('', `Top combinations by ${value.rankBy} (Δ vs best single < 0 = predicted better than its better single agent; model prediction, not experimental synergy):`);
                    lines.push(table(['combination', 'Δ vs best single', 'improved', 'worse', 'per-background Δ vs best single'], value.combinations.map((c) => [
                        c.combination,
                        signed(c.synergyMean),
                        `${c.improved}/${c.backgroundsScored}`,
                        c.worse,
                        c.perBackground.map((p) => `${p.background}:${num(p.additive, 1)}`).join(' '),
                    ]), { align: ['l', 'r', 'r', 'r', 'l'] }));
                }
                if (Array.isArray(value.singles) && value.singles.length > 0) {
                    lines.push('', 'Top single agents (negative delta = better than untreated):');
                    lines.push(table(['node', 'prolif', 'death', 'difference', 'delta'], value.singles.map((s) => [
                        s.node, num(s.meanProliferation), num(s.meanDeath), num(s.meanDifference), signed(s.deltaVsBaseline),
                    ]), { align: ['l', 'r', 'r', 'r', 'r'] }));
                }
                if (value.csvPath) lines.push('', `Full table written to ${value.csvPath}`);
                if (Array.isArray(value.notes) && value.notes.length > 0) lines.push('', ...value.notes);
                return [{ type: 'text', text: lines.join('\n') }];
            },
        },
        timeoutMs: 900000,
        isConcurrencySafe: () => true,
        async execute(args, exec) {
            const backgrounds = Array.isArray(args.backgrounds) && args.backgrounds.length > 0
                ? args.backgrounds
                : analyst.tables.publishedBackgrounds.map((b) => b.name);
            const mode = args.mode ?? 'mono';
            const levels = Array.isArray(args.levels) && args.levels.length > 0 ? args.levels.map(Number) : [0];
            const maxAnalyses = Number.isInteger(args.max_analyses) ? args.max_analyses : defaultMaxAnalyses;

            const enumeration = enumerateJobs(analyst, {
                backgrounds,
                mode,
                candidates: args.candidates,
                exclude: args.exclude,
                levels,
            });
            const { jobs } = enumeration;
            const notes = [...enumeration.notes];
            if (jobs.length > maxAnalyses) {
                throw new Error(`this screen needs ${jobs.length} model runs, above the ${maxAnalyses} cap. `
                    + 'Narrow it with `candidates`, fewer `backgrounds`, or raise `max_analyses`.');
            }
            notes.push(`Estimated ${estimateCost(jobs.length, 1.45, workers).toFixed(1)}s across ${workers} worker threads.`);

            const started = Date.now();
            const results = await runJobs(analyst, jobs, {
                workers,
                onProgress: exec?.onProgress,
            });
            const scored = scoreScreen(analyst, results, {
                minBackgroundsImproved: Number.isInteger(args.min_backgrounds_improved)
                    ? args.min_backgrounds_improved : undefined,
            });
            const elapsed = (Date.now() - started) / 1000;

            const rankBy = args.rank_by ?? (mode === 'mono' ? 'benefit' : 'synergy');
            const top = Number.isInteger(args.top) ? Math.max(1, Math.min(200, args.top)) : 25;

            const backgroundsSeen = [...new Set(results.map((r) => r.background))];
            const baselineBy = new Map(scored.baselines.map((b) => [b.background, b]));

            const combos = rankCombinations(scored.combinations, rankBy).slice(0, top).map((c) => ({
                // `node@level` tokens keep different intervention strengths distinct and
                // the pair symmetric; this is the same key the aggregator and CSV use.
                combination: combinationKey(c.combination, c.levels),
                synergyMean: c.additiveMean,
                synergySum: c.additiveSum,
                improved: c.improved,
                worse: c.worse,
                equal: c.equal,
                backgroundsScored: c.backgrounds,
                perBackground: c.perBackground.map((r) => ({
                    background: r.background,
                    ...(r.additive !== undefined ? { additive: r.additive } : {}),
                    difference: r.difference,
                    proliferation: r.proliferation,
                    death: r.death,
                })),
            }));

            // Aggregate singles per node across backgrounds.
            /** @type {Map<string, {node:string, rows:Array<any>}>} */
            const byNode = new Map();
            for (const s of scored.singles) {
                if (!byNode.has(s.node)) byNode.set(s.node, { node: s.node, rows: [] });
                byNode.get(s.node).rows.push(s);
            }
            const singles = [...byNode.values()].map(({ node, rows }) => {
                const n = rows.length || 1;
                const meanDifference = rows.reduce((a, r) => a + r.difference, 0) / n;
                const meanProliferation = rows.reduce((a, r) => a + r.proliferation, 0) / n;
                const meanDeath = rows.reduce((a, r) => a + r.death, 0) / n;
                const baselineMean = rows.reduce((a, r) => a + (baselineBy.get(r.background)?.difference ?? 0), 0) / n;
                return {
                    node,
                    meanDifference,
                    meanProliferation,
                    meanDeath,
                    deltaVsBaseline: meanDifference - baselineMean,
                    perBackground: rows.map((r) => ({
                        background: r.background,
                        difference: r.difference,
                        proliferation: r.proliferation,
                        death: r.death,
                    })),
                };
            }).sort((a, b) => a.meanDifference - b.meanDifference).slice(0, top);

            const response = {
                evidence_kind: 'model_prediction',
                evidence_note: MODEL_PREDICTION_NOTE,
                mode,
                backgrounds: backgroundsSeen,
                levels,
                candidateCount: enumeration.candidateIds.length,
                analysesRun: jobs.length,
                elapsedSeconds: elapsed,
                rankBy,
                baselines: scored.baselines.map((b) => ({
                    background: b.background,
                    proliferation: b.proliferation,
                    death: b.death,
                    difference: b.difference,
                })),
                ...(mode === 'mono' ? {} : { combinations: combos }),
                ...(mode === 'pair' ? {} : { singles }),
                notes,
            };

            if (typeof args.output_csv === 'string' && args.output_csv.trim() !== '') {
                const path = resolve(process.cwd(), args.output_csv);
                const header = ['background', 'kind', 'a', 'level_a', 'b', 'level_b',
                    'proliferation', 'death', 'difference', 'additive', 'uncertainty', 'stabilizing', 'iterations'];
                const additiveBy = new Map();
                for (const row of scored.rows) {
                    additiveBy.set(pairKey(row.background, row.combination, row.levels), row.additive);
                }
                const lines = [header.join(',')];
                const csvAdditive = new Map();
                for (const r of results) {
                    const nameA = r.a === undefined ? '' : model.names.get(r.a);
                    const nameB = r.b === undefined ? '' : model.names.get(r.b);
                    const levels = [r.levelA, r.levelB];
                    const key = r.kind === 'pair' ? pairKey(r.background, [nameA, nameB], levels) : undefined;
                    const additive = key === undefined ? undefined : additiveBy.get(key);
                    if (key !== undefined) csvAdditive.set(key, additive);
                    lines.push([
                        r.background, r.kind, nameA, r.levelA ?? '', nameB, r.levelB ?? '',
                        r.proliferation, r.death, r.difference, additive, r.uncertainty, r.stabilizing, r.iterations,
                    ].map(csvField).join(','));
                }
                // Same-source check: the CSV `additive` column and the in-memory scores must
                // be keyed identically (background|target@level, symmetric). Any drift is a bug.
                const mismatches = [];
                for (const [key, want] of additiveBy) {
                    if (!csvAdditive.has(key)) {
                        mismatches.push(`missing key ${key}`);
                        continue;
                    }
                    const got = csvAdditive.get(key);
                    if (typeof want === 'number' && typeof got === 'number') {
                        if (Math.abs(want - got) > 1e-12) mismatches.push(`${key}: csv=${got} memory=${want}`);
                    } else if (want !== got) {
                        mismatches.push(`${key}: csv=${got} memory=${want}`);
                    }
                }
                if (additiveBy.size !== csvAdditive.size) {
                    mismatches.push(`key count csv=${csvAdditive.size} memory=${additiveBy.size}`);
                }
                if (mismatches.length > 0) {
                    throw new Error(`in-memory vs CSV additive keys diverged: ${mismatches.slice(0, 5).join('; ')}`);
                }
                writeFileSync(path, `${lines.join('\n')}\n`, 'utf8');
                response.csvPath = path;
            }

            return response;
        },
    });
}
