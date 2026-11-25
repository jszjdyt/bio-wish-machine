/**
 * The NSCLC "virtual tumor" analysis layer: genotype resolution, steady-state
 * prediction, mono/pair screens and the paper's additivity score.
 *
 * Screening protocol reproduced from the manuscript's own analysis scripts
 * (`src/plot/single_screen_plot.R`, `scripts/figures/main/figure_5_fancy_circle.R`):
 *
 *  - a "drug" is a knockout: the node is fixed to `RangeFrom` (level 0 for
 *    every node of this model);
 *  - `Difference = mean(Proliferation) - mean(Death)`, so **lower is better**;
 *  - per background, `additive(pair) = Difference(pair) - min(Difference(a),
 *    Difference(b))` over the two single knockouts, i.e. how much the
 *    combination improves on the *better* of its components.
 *    **Negative additive = synergy, positive = antagonism**;
 *  - the paper drops combinations that fail to beat the best single agent in at
 *    least `n-1` of the 7 published backgrounds, which is exposed here as
 *    `minBackgroundsImproved`.
 *
 * @module dsh-nsclc-vt/analysis
 */
import { loadModel, withKnockouts, formulaOf } from './model.js';
import { analyzeStability, readout } from './vmcai.js';
import { loadTables, DEFAULT_MODEL_PATH } from './data.js';
import { freeVars } from './expr.js';

/** Phenotypes the model reports and the paper screens on. */
export const DEFAULT_PHENOTYPES = ['Proliferation', 'Death'];

/** Normalise a user-supplied node token for lookup. */
const norm = (s) => String(s).trim().toLowerCase().replace(/[\s_-]+/g, '');

/**
 * Build an analyst bound to one model file.
 * @param {{modelPath?: string, phenotypes?: string[]}} [opts]
 */
export function createAnalyst(opts = {}) {
    const modelPath = opts.modelPath ?? DEFAULT_MODEL_PATH;
    const phenotypes = opts.phenotypes ?? DEFAULT_PHENOTYPES;
    const model = loadModel(modelPath);
    const tables = loadTables();

    // ---- annotation -------------------------------------------------------
    const geneByNode = new Map();
    const proteinByNode = new Map();
    for (const r of tables.nodesInfo) {
        if (r.node) {
            geneByNode.set(r.node, r.gene);
            proteinByNode.set(r.node, r.protein);
        }
        if (r.gene) {
            const id = model.byName.get(r.node);
            if (id !== undefined && !model.byKey.has(r.gene.toLowerCase())) {
                model.byKey.set(r.gene.toLowerCase(), id);
            }
        }
    }
    const curatedFunction = new Map();
    for (const r of tables.targetFunctions) curatedFunction.set(r.node, r);

    const phenotypeIds = phenotypes.map((p) => {
        const id = model.byName.get(p);
        if (id === undefined) throw new Error(`the model has no phenotype node named ${JSON.stringify(p)}`);
        return id;
    });

    const primary = {
        proliferation: model.byName.get('Proliferation'),
        death: model.byName.get('Death'),
    };

    /**
     * Resolve a node name / gene symbol / id to a model variable.
     * @param {string|number} token
     * @returns {{id:number, name:string, gene?:string}}
     */
    function resolveNode(token) {
        if (typeof token === 'number' || /^\d+$/.test(String(token))) {
            const id = Number(token);
            if (model.names.has(id)) return describeNode(id);
        }
        const raw = String(token).trim();
        if (model.byName.has(raw)) return describeNode(model.byName.get(raw));
        const key = raw.toLowerCase();
        if (model.byKey.has(key)) {
            const id = model.byKey.get(key);
            if (id >= 0) return describeNode(id);
        }
        const stripped = norm(raw);
        const fuzzy = model.order.filter((id) => norm(model.names.get(id)) === stripped);
        if (fuzzy.length === 1) return describeNode(fuzzy[0]);
        const byGene = model.order.filter((id) => norm(geneByNode.get(model.names.get(id)) ?? '') === stripped);
        if (byGene.length === 1) return describeNode(byGene[0]);
        throw new Error(`unknown node or gene ${JSON.stringify(raw)}`
            + (byGene.length > 1 ? ` (matches ${byGene.map((i) => model.names.get(i)).join(', ')})` : ''));
    }

    function describeNode(id) {
        return { id, name: model.names.get(id), gene: geneByNode.get(model.names.get(id)) };
    }

    /**
     * Turn `"RAS:3"`, `{node, level}`, or `"RAS"` into `[id, level]`.
     * A bare name means a knockout at the node's lower bound.
     */
    function parsePerturbation(spec) {
        let token, level;
        if (typeof spec === 'string') {
            const m = /^(.*?)\s*[:=]\s*(-?\d+)$/.exec(spec.trim());
            if (m) { token = m[1]; level = Number(m[2]); } else { token = spec; level = undefined; }
        } else if (spec && typeof spec === 'object') {
            token = spec.node ?? spec.name ?? spec.gene;
            level = spec.level ?? spec.activity;
        } else {
            throw new Error(`cannot read perturbation ${JSON.stringify(spec)}`);
        }
        if (token === undefined) throw new Error(`perturbation ${JSON.stringify(spec)} has no node`);
        const { id, name } = resolveNode(token);
        const [lo, hi] = model.range.get(id);
        const resolved = level === undefined ? lo : Number(level);
        if (!Number.isInteger(resolved) || resolved < lo || resolved > hi) {
            throw new Error(`level ${resolved} for ${name} is outside its range [${lo}, ${hi}]`);
        }
        return { id, name, level: resolved };
    }

    /**
     * Resolve a background: a published cell-line name, a map of node -> level,
     * or an array of `"NODE:level"` strings.
     * @returns {{name?:string, mutations:Array<{id:number,name:string,level:number}>}}
     */
    function resolveBackground(spec) {
        if (spec === undefined || spec === null || spec === '') {
            return { mutations: [] };
        }
        if (typeof spec === 'string') {
            const wanted = norm(spec);
            const hit = tables.publishedBackgrounds.find((b) => norm(b.name) === wanted);
            if (hit) {
                return { name: hit.name, mutations: hit.mutations.map(parsePerturbation) };
            }
            // Fall through: allow "RAS:3,STK11:0" as a background string.
            if (/[:=]/.test(spec)) {
                return { mutations: spec.split(/[,;]/).filter((s) => s.trim() !== '').map(parsePerturbation) };
            }
            throw new Error(`unknown background ${JSON.stringify(spec)}; known: `
                + tables.publishedBackgrounds.map((b) => b.name).join(', '));
        }
        if (Array.isArray(spec)) return { mutations: spec.map(parsePerturbation) };
        if (typeof spec === 'object') {
            const entries = Object.entries(spec);
            if (entries.length > 0 && entries.every(([, v]) => typeof v === 'number')) {
                return { mutations: entries.map(([node, level]) => parsePerturbation({ node, level })) };
            }
            if (Array.isArray(spec.mutations)) {
                return { name: spec.name, mutations: spec.mutations.map(parsePerturbation) };
            }
        }
        throw new Error(`cannot read background ${JSON.stringify(spec)}`);
    }

    /**
     * Run one VMCAI analysis.
     * @param {{background?: any, interventions?: any[], recordHistory?: boolean, maxIterations?: number}} req
     */
    function solve(req = {}) {
        const bg = resolveBackground(req.background);
        const perms = (req.interventions ?? []).map(parsePerturbation);
        const ko = new Map(bg.mutations.map((m) => [m.id, m.level]));
        // Perturbations win over the background — the paper's NANSEN pipeline
        // puts the perturbation flags first on the BioCheckConsole command line
        // and the first `-ko` for a node wins.
        for (const p of perms) ko.set(p.id, p.level);
        const analysis = analyzeStability(withKnockouts(model, ko), {
            recordHistory: req.recordHistory === true,
            maxIterations: req.maxIterations,
        });
        return { background: bg, interventions: perms, ko, analysis };
    }

    /** Phenotype readouts plus the paper's `Difference`. */
    function phenotypeReadout(analysis) {
        const out = phenotypeIds.map((id) => ({
            node: model.names.get(id),
            ...readout(model, analysis.bounds, id),
        }));
        const proliferation = primary.proliferation === undefined ? undefined : readout(model, analysis.bounds, primary.proliferation);
        const death = primary.death === undefined ? undefined : readout(model, analysis.bounds, primary.death);
        return {
            phenotypes: out,
            difference: proliferation !== undefined && death !== undefined
                ? proliferation.mean - death.mean
                : undefined,
            proliferation,
            death,
        };
    }

    return {
        model,
        tables,
        phenotypes,
        phenotypeIds,
        modelPath,
        resolveNode,
        describeNode,
        parsePerturbation,
        resolveBackground,
        solve,
        phenotypeReadout,
        geneOf: (id) => geneByNode.get(model.names.get(id)),
        proteinOf: (id) => proteinByNode.get(model.names.get(id)),
        curatedFunctionOf: (id) => curatedFunction.get(model.names.get(id)),
        formulaOf: (id, o) => formulaOf(model, id, o),
        isDefaultFunction: (id) => model.isDefaultExpr.get(id) === true,
        inputsOf: (id) => model.inputs.get(id),
        activatorsOf: (id) => model.activators.get(id),
        inhibitorsOf: (id) => model.inhibitors.get(id),
        outputsOf: (id) => [...model.outputs.get(id)],
        freeVarsOf: (id) => [...freeVars(model.exprs.get(id))],
    };
}

/**
 * Rank helper: `benefit` mirrors `Difference` with the sign flipped so that
 * higher is therapeutically better.
 * @param {number|undefined} difference
 */
export const benefitOf = (difference) => (difference === undefined ? undefined : -difference);
