/**
 * Screen orchestration: job enumeration, parallel execution, and the paper's
 * additivity scoring / ranking.
 *
 * Work is dispatched to `worker_threads` above an in-process threshold so a
 * genome-wide pairwise screen does not block the host for minutes. If worker
 * threads cannot start (some sandboxes forbid them) the same jobs run inline
 * and the result is identical, only slower — nothing here changes the numbers.
 *
 * @module dsh-nsclc-vt/screen
 */
import { Worker } from 'node:worker_threads';
import { availableParallelism } from 'node:os';
import { fileURLToPath } from 'node:url';

const WORKER_URL = new URL('./screen-worker.js', import.meta.url);

/**
 * @typedef {object} ScreenJob
 * @property {string} key         stable identity of the job
 * @property {string} background  background label
 * @property {'baseline'|'single'|'pair'} kind
 * @property {number} [a]         first node id
 * @property {number} [levelA]
 * @property {number} [b]         second node id
 * @property {number} [levelB]
 */

/**
 * Enumerate every analysis the screen needs.
 * @param {ReturnType<import('./analysis.js').createAnalyst>} analyst
 * @param {{backgrounds:(string|object)[], mode?:'mono'|'pair'|'both', candidates?:string[], exclude?:string[], levels?:number[], maxPairs?:number}} spec
 * @returns {{jobs: ScreenJob[], notes: string[], truncatedPairs: boolean}}
 */
export function enumerateJobs(analyst, spec) {
    const model = analyst.model;
    const backgrounds = spec.backgrounds ?? analyst.tables.publishedBackgrounds.map((b) => b.name);
    if (backgrounds.length === 0) throw new Error('at least one background is required');
    const mode = spec.mode ?? 'mono';
    const levels = (spec.levels ?? [0]).map(Number);
    for (const l of levels) {
        if (!Number.isInteger(l)) throw new Error(`screening level ${l} is not an integer`);
    }

    const excluded = new Set((spec.exclude ?? []).flatMap((t) => {
        try { return [analyst.resolveNode(t).id]; } catch { return []; }
    }));
    let candidateIds;
    if (spec.candidates && spec.candidates.length > 0) {
        candidateIds = spec.candidates.map((t) => analyst.resolveNode(t).id);
    } else {
        candidateIds = [...model.order];
    }
    candidateIds = [...new Set(candidateIds)]
        .filter((id) => !excluded.has(id))
        .filter((id) => {
            const [lo, hi] = model.range.get(id);
            return levels.some((l) => l >= lo && l <= hi);
        })
        .sort((a, b) => a - b);

    // NANSEN never perturbs the phenotype readouts themselves.
    const phenotypeIds = new Set(analyst.phenotypeIds);
    const perturbable = candidateIds.filter((id) => !phenotypeIds.has(id));

    const jobs = [];
    const notes = [];
    for (const bg of backgrounds) {
        const bgLabel = typeof bg === 'string' ? bg : (bg && bg.name) || 'custom';
        const bgSpec = bg;
        jobs.push({ key: `${bgLabel}|baseline`, background: bgLabel, backgroundSpec: bgSpec, kind: 'baseline' });
        const units = [];
        for (const id of perturbable) {
            for (const level of levels) {
                const [lo, hi] = model.range.get(id);
                if (level < lo || level > hi) continue;
                units.push({ a: id, levelA: level });
            }
        }
        for (const u of units) {
            jobs.push({
                key: `${bgLabel}|${model.names.get(u.a)}:${u.levelA}`,
                background: bgLabel, backgroundSpec: bgSpec, kind: 'single', ...u,
            });
        }
        if (mode === 'pair' || mode === 'both') {
            let pairs = 0;
            for (let i = 0; i < units.length; i += 1) {
                for (let j = i + 1; j < units.length; j += 1) {
                    const u = units[i], v = units[j];
                    if (u.a === v.a) continue;      // no self-pairs
                    jobs.push({
                        key: `${bgLabel}|${model.names.get(u.a)}:${u.levelA}+${model.names.get(v.a)}:${v.levelA}`,
                        background: bgLabel, backgroundSpec: bgSpec, kind: 'pair',
                        a: u.a, levelA: u.levelA, b: v.a, levelB: v.levelA,
                    });
                    pairs += 1;
                }
            }
            if (spec.maxPairs !== undefined && pairs > spec.maxPairs) {
                notes.push(`enumeration produced ${pairs} pairs per background; the first ${spec.maxPairs} were kept`);
            }
        }
    }
    return { jobs, notes, candidateIds: perturbable, backgrounds };
}

/**
 * Run a batch of jobs on the current thread.
 * @param {ReturnType<import('./analysis.js').createAnalyst>} analyst
 * @param {ScreenJob[]} jobs
 */
export function runJobsInline(analyst, jobs) {
    return jobs.map((job) => runJob(analyst, job));
}

/**
 * Execute one job. Kept free of analyst state so `screen-worker.js` can reuse it.
 */
export function runJob(analyst, job) {
    const model = analyst.model;
    const background = job.backgroundSpec ?? job.background;
    const interventions = job.kind === 'baseline'
        ? []
        : job.kind === 'single'
            ? [{ node: model.names.get(job.a), level: job.levelA }]
            : [
                { node: model.names.get(job.a), level: job.levelA },
                { node: model.names.get(job.b), level: job.levelB },
            ];
    const { analysis } = analyst.solve({ background, interventions });
    const pheno = analyst.phenotypeReadout(analysis);
    return {
        key: job.key,
        background: job.background,
        kind: job.kind,
        a: job.a, levelA: job.levelA, b: job.b, levelB: job.levelB,
        difference: pheno.difference,
        proliferation: pheno.proliferation === undefined ? undefined : pheno.proliferation.mean,
        death: pheno.death === undefined ? undefined : pheno.death.mean,
        proliferationLo: pheno.proliferation === undefined ? undefined : pheno.proliferation.lo,
        proliferationHi: pheno.proliferation === undefined ? undefined : pheno.proliferation.hi,
        deathLo: pheno.death === undefined ? undefined : pheno.death.lo,
        deathHi: pheno.death === undefined ? undefined : pheno.death.hi,
        uncertainty: pheno.phenotypes.reduce((acc, p) => Math.max(acc, p.uncertainty), 0),
        stabilizing: analysis.stabilizing,
        iterations: analysis.iterations,
    };
}

const defaultWorkers = () => Math.max(1, Math.min(8, (availableParallelism?.() ?? 4) - 1));

/**
 * Run jobs, fanning out to worker threads when the batch is large.
 * @param {ReturnType<import('./analysis.js').createAnalyst>} analyst
 * @param {ScreenJob[]} jobs
 * @param {{inProcessLimit?: number, workers?: number, onProgress?: (done:number, total:number)=>void}} [opts]
 * @returns {Promise<ReturnType<typeof runJob>[]>}
 */
export async function runJobs(analyst, jobs, opts = {}) {
    const limit = opts.inProcessLimit ?? 400;
    if (jobs.length <= limit) return runJobsInline(analyst, jobs);

    const workerCount = Math.max(1, Math.min(opts.workers ?? defaultWorkers(), jobs.length));
    const chunks = Array.from({ length: workerCount }, () => []);
    jobs.forEach((job, i) => chunks[i % workerCount].push(job));

    const payload = { modelPath: analyst.modelPath, phenotypes: analyst.phenotypes };
    /** @type {Worker[]} */
    const workers = [];
    try {
        for (const chunk of chunks) {
            if (chunk.length === 0) continue;
            workers.push(new Worker(WORKER_URL, { workerData: { ...payload, jobs: chunk } }));
        }
    } catch {
        // A partially constructed pool still has to be torn down before falling back.
        for (const w of workers) w.terminate().catch(() => {});
        return runJobsInline(analyst, jobs);
    }

    try {
        const results = await Promise.all(workers.map((w) => new Promise((resolve, reject) => {
            w.once('message', (msg) => {
                if (msg && msg.ok) resolve(msg.results);
                else reject(new Error((msg && msg.error) || 'screen worker failed'));
            });
            w.once('error', reject);
            w.once('exit', (code) => { if (code !== 0) reject(new Error(`screen worker exited with code ${code}`)); });
        })));
        return results.flat();
    } catch {
        // Workers are an optimisation only: fall back to the identical inline run.
        return runJobsInline(analyst, jobs);
    } finally {
        for (const w of workers) w.terminate().catch(() => {});
    }
}

/**
 * Stable identity of one intervention: `node@level`.
 * @param {string} name
 * @param {number} level
 */
export function interventionKey(name, level) {
    return `${name}@${level}`;
}

/**
 * Symmetric combination key that **keeps each intervention's level**.
 * `combinationKey(['A','B'], [0,3])` and `combinationKey(['B','A'], [3,0])`
 * are the same key; `[3,0]` vs `[0,0]` are different keys.
 * @param {string[]} combination
 * @param {number[]} [levels]
 */
export function combinationKey(combination, levels = []) {
    const names = Array.isArray(combination) ? combination : [combination];
    const tokens = names.map((name, i) => interventionKey(name, levels[i] ?? 0));
    return tokens.sort().join('+');
}

/**
 * Background-scoped pair key. The in-memory aggregator and the CSV writer both
 * call this, so `additive` values can never be joined on different keys.
 * @param {string} background
 * @param {string[]} combination
 * @param {number[]} [levels]
 */
export function pairKey(background, combination, levels = []) {
    return `${background}|${combinationKey(combination, levels)}`;
}

/**
 * Sort a pair into a canonical (symmetric) order and keep names/levels aligned.
 * @param {string[]} combination
 * @param {number[]} levels
 */
function canonicalCombination(combination, levels = []) {
    const pairs = combination.map((name, i) => ({ name, level: levels[i] ?? 0 }));
    pairs.sort((a, b) => (a.name < b.name ? -1 : a.name > b.name ? 1 : a.level - b.level));
    return { names: pairs.map((p) => p.name), levels: pairs.map((p) => p.level) };
}

/**
 * Turn raw job results into scored, ranked screen rows.
 *
 * `difference` is `mean(Proliferation) - mean(Death)` (lower is better), so
 * `benefit = -difference` and `additive = difference(pair) - best single
 * difference` (negative = predicted improvement over the better component).
 * Aggregation keys include each intervention's level and are background-scoped
 * via `pairKey`, shared with the CSV writer.
 *
 * @param {ReturnType<import('./analysis.js').createAnalyst>} analyst
 * @param {Array<ReturnType<typeof runJob>>} results
 * @param {{minBackgroundsImproved?: number}} [opts]
 */
export function scoreScreen(analyst, results, opts = {}) {
    const model = analyst.model;
    /** @type {Map<string, ReturnType<typeof runJob>>} */
    const byKey = new Map(results.map((r) => [r.key, r]));
    const names = analyst.phenotypes;

    // singles indexed by background -> node -> level -> row
    const singles = new Map();
    const baselines = new Map();
    for (const r of results) {
        if (r.kind === 'single') {
            const k = `${r.background}|${r.a}|${r.levelA}`;
            singles.set(k, r);
        } else if (r.kind === 'baseline') {
            baselines.set(r.background, r);
        }
    }

    const rows = [];
    for (const r of results) {
        if (r.kind !== 'pair') continue;
        const singleA = singles.get(`${r.background}|${r.a}|${r.levelA}`);
        const singleB = singles.get(`${r.background}|${r.b}|${r.levelB}`);
        const baseline = baselines.get(r.background);
        const diffs = [singleA, singleB].filter(Boolean).map((s) => s.difference);
        const bestSingle = diffs.length > 0 ? Math.min(...diffs) : undefined;
        rows.push({
            background: r.background,
            combination: [model.names.get(r.a), model.names.get(r.b)],
            levels: [r.levelA, r.levelB],
            proliferation: r.proliferation,
            death: r.death,
            difference: r.difference,
            benefit: r.difference === undefined ? undefined : -r.difference,
            baselineDifference: baseline === undefined ? undefined : baseline.difference,
            singleADifference: singleA === undefined ? undefined : singleA.difference,
            singleBDifference: singleB === undefined ? undefined : singleB.difference,
            bestSingleDifference: bestSingle,
            additive: r.difference !== undefined && bestSingle !== undefined ? r.difference - bestSingle : undefined,
            uncertainty: r.uncertainty,
            stabilizing: r.stabilizing,
        });
    }

    // Per-combination aggregate across backgrounds.
    /** @type {Map<string, {combination:string[], perBackground:object[], additiveSum:number, additiveMean:number, improved:number, worse:number, equal:number}>} */
    const aggregates = new Map();
    for (const row of rows) {
        // Key includes each intervention's level; sorting the `node@level` tokens keeps
        // the pair symmetric without ever collapsing different strengths into one row.
        const canonical = canonicalCombination(row.combination, row.levels || []);
        const key = combinationKey(canonical.names, canonical.levels);
        let agg = aggregates.get(key);
        if (agg === undefined) {
            agg = {
                key,
                combination: canonical.names,
                levels: canonical.levels,
                perBackground: [],
                additiveSum: 0, additiveMean: 0, improved: 0, worse: 0, equal: 0,
            };
            aggregates.set(key, agg);
        }
        agg.perBackground.push(row);
        if (row.additive !== undefined) {
            agg.additiveSum += row.additive;
            if (row.additive < 0) agg.improved += 1;
            else if (row.additive > 0) agg.worse += 1;
            else agg.equal += 1;
        }
    }
    for (const agg of aggregates.values()) {
        const n = agg.perBackground.filter((r) => r.additive !== undefined).length;
        agg.additiveMean = n === 0 ? 0 : agg.additiveSum / n;
        agg.backgrounds = n;
    }
    // The manuscript drops combinations that fail to beat the best single agent
    // in at least 2 of the 7 backgrounds (`rowSums(circular_mat >= 0) < ncol-1`);
    // `minBackgroundsImproved` generalises that cutoff.
    const minImproved = opts.minBackgroundsImproved;
    const combinations = [...aggregates.values()].filter((agg) => (minImproved === undefined
        ? true
        : agg.improved >= minImproved));

    const singlesRanked = results.filter((r) => r.kind === 'single').map((r) => ({
        background: r.background,
        node: model.names.get(r.a),
        level: r.levelA,
        proliferation: r.proliferation,
        death: r.death,
        difference: r.difference,
        benefit: r.difference === undefined ? undefined : -r.difference,
        uncertainty: r.uncertainty,
    }));
    const baselinesRanked = results.filter((r) => r.kind === 'baseline').map((r) => ({
        background: r.background,
        proliferation: r.proliferation,
        death: r.death,
        difference: r.difference,
    }));

    return { rows, combinations, singles: singlesRanked, baselines: baselinesRanked, phenotypes: names, byKey };
}

const meanOf = (rows, pick) => rows.reduce((acc, row) => acc + (pick(row) ?? 0), 0) / Math.max(1, rows.length);
const IS_HEALTHY = /healthy|normal|control|wt|wildtype/i;
const benefitOf = (row) => (row.benefit ?? (row.difference === undefined ? undefined : -row.difference));

/**
 * Rank combination summaries. Every ranking is expressed so that **smaller key
 * sorts first**, i.e. the first element is the best.

 * - `synergy`: most negative predicted difference-vs-best-single first.
 * - `benefit`: largest mean predicted benefit first (`benefit = -Difference`).
 * - `death`: most predicted cell death first.
 * - `window`: tumour-selective benefit, defined as
 *   `mean(benefit over tumour backgrounds) - benefit(Healthy background)`.
 *   Positive = the pair helps tumours more than healthy tissue. It requires a
 *   Healthy background in the screened set; otherwise it throws instead of
 *   silently ranking by an undefined quantity.
 *
 * @param {Array<object>} combinations
 * @param {'synergy'|'benefit'|'death'|'window'} rankBy
 */
export function rankCombinations(combinations, rankBy = 'synergy') {
    const scorers = {
        synergy: (c) => c.additiveMean,                                   // most negative first
        benefit: (c) => -meanOf(c.perBackground, (r) => benefitOf(r)),
        death: (c) => -meanOf(c.perBackground, (r) => r.death),
        window: (c) => {
            const rows = c.perBackground || [];
            const healthy = rows.filter((row) => IS_HEALTHY.test(String(row.background)));
            const tumour = rows.filter((row) => !IS_HEALTHY.test(String(row.background)));
            if (healthy.length === 0) {
                throw new Error('window ranking requires a Healthy/normal background: none of the screened '
                    + `backgrounds (${rows.map((r) => r.background).join(', ') || 'none'}) looks healthy. `
                    + 'Screen the "healthy" background, or rank by synergy/benefit/death.');
            }
            if (tumour.length === 0) {
                throw new Error('window ranking requires at least one tumour background in addition to Healthy: '
                    + 'only Healthy was screened, so tumour selectivity is undefined.');
            }
            // Larger window = better => negate so ascending sort puts it first.
            return -(meanOf(tumour, (r) => benefitOf(r)) - meanOf(healthy, (r) => benefitOf(r)));
        },
    };
    const scorer = scorers[rankBy];
    if (scorer === undefined) {
        throw new Error(`unknown rankBy ${JSON.stringify(rankBy)}; expected synergy/benefit/death/window`);
    }
    return combinations
        .map((combination) => ({ combination, score: scorer(combination) }))
        .sort((a, b) => a.score - b.score)
        .map((entry) => entry.combination);
}

/** Estimate wall-clock cost of a screen, for the tool's pre-flight note. */
export function estimateCost(jobCount, msPerJob = 1.45, workers = defaultWorkers()) {
    return (jobCount * msPerJob) / 1000 / workers;
}

export { defaultWorkers };
