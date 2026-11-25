/**
 * Loader for BMA (BioModelAnalyzer) qualitative-network JSON models, ported from
 * `src/BioCheckAnalyzer/Marshal.fs`.
 *
 * Two details of the original importer are load-bearing and reproduced here:
 *
 *  - A variable whose `Formula` is the empty string gets a **synthesised**
 *    default target function built from its relationship edges:
 *      * no activators, no inhibitors -> `Const(rangeMin)`
 *      * inhibitors only             -> `Const(rangeMin)`
 *      * activators only             -> `min(rangeMax, avg(activators))`
 *      * both                        -> `max(rangeMin, avg(activators) - avg(inhibitors))`
 *  - A node's **inputs** are its relationship activators followed by its
 *    relationship inhibitors — *not* the free variables of its formula. A
 *    formula input that is not a declared relationship input therefore
 *    evaluates to the node's lower bound during analysis.
 *
 * @module dsh-nsclc-vt/model
 */
import { readFileSync } from 'node:fs';
import { parseFormula, formatFormula } from './expr.js';

/**
 * @typedef {object} NetworkModel
 * @property {string} name
 * @property {number[]} order            variable ids, ascending
 * @property {Map<number,string>} names  id -> node name
 * @property {Map<number,[number,number]>} range id -> [lo, hi]
 * @property {Map<number,Array>} exprs   id -> target-function AST (original formula)
 * @property {Map<number,Array>} effectiveExprs id -> AST actually used (incl. defaults)
 * @property {Map<number,boolean>} isDefaultExpr id -> true when the formula was synthesised
 * @property {Map<number,number[]>} inputs id -> relationship input ids (activators then inhibitors)
 * @property {Map<number,number[]>} activators
 * @property {Map<number,number[]>} inhibitors
 * @property {Map<number,Set<number>>} outputs id -> ids that list it as an input
 * @property {Map<string,number>} byName  exact node name -> id
 * @property {Map<string,number>} byKey   lowercase name / gene symbol -> id
 * @property {Map<number,{gene?:string,protein?:string}>} annotation
 */

/** Coerce a JSON scalar that may be a number or a numeric string. */
const toInt = (v) => (typeof v === 'number' ? v : Number(String(v).trim()));

/**
 * Parse a BMA model JSON document.
 * @param {string|object} source path to the JSON file, or the parsed document
 * @returns {NetworkModel}
 */
export function loadModel(source) {
    const doc = typeof source === 'string' ? JSON.parse(readFileSync(source, 'utf8')) : source;
    if (!doc || typeof doc !== 'object' || !doc.Model || !Array.isArray(doc.Model.Variables)) {
        throw new Error('not a BMA model document: expected a top-level "Model" with a "Variables" array');
    }
    const variables = doc.Model.Variables;
    const relationships = Array.isArray(doc.Model.Relationships) ? doc.Model.Relationships : [];

    /** @type {Map<number,[number,number]>} */ const range = new Map();
    /** @type {Map<number,string>} */ const names = new Map();
    /** @type {Map<number,Array>} */ const exprs = new Map();
    /** @type {Map<number,boolean>} */ const isDefaultExpr = new Map();
    /** @type {Map<number,number[]>} */ const activators = new Map();
    /** @type {Map<number,number[]>} */ const inhibitors = new Map();
    const order = [];

    for (const v of variables) {
        const id = toInt(v.Id);
        if (!Number.isInteger(id)) throw new Error(`variable without a usable Id: ${JSON.stringify(v.Name)}`);
        if (range.has(id)) throw new Error(`duplicate variable id ${id}`);
        const lo = toInt(v.RangeFrom), hi = toInt(v.RangeTo);
        if (!Number.isInteger(lo) || !Number.isInteger(hi) || hi < lo) {
            throw new Error(`variable ${v.Name} (id ${id}) has an invalid range [${v.RangeFrom}, ${v.RangeTo}]`);
        }
        range.set(id, [lo, hi]);
        names.set(id, String(v.Name));
        activators.set(id, []);
        inhibitors.set(id, []);
        order.push(id);
    }
    order.sort((a, b) => a - b);

    for (const r of relationships) {
        const from = toInt(r.FromVariable), to = toInt(r.ToVariable);
        if (!range.has(from) || !range.has(to)) {
            throw new Error(`relationship ${r.Id} references unknown variable (${from} -> ${to})`);
        }
        const kind = String(r.Type);
        if (kind === 'Activator') activators.get(to).push(from);
        else if (kind === 'Inhibitor') inhibitors.get(to).push(from);
        else throw new Error(`relationship ${r.Id} has an unknown type ${JSON.stringify(kind)}`);
    }

    /** @type {Map<number,number[]>} */ const inputs = new Map();
    /** @type {Map<number,Set<number>>} */ const outputs = new Map();
    for (const id of order) {
        inputs.set(id, [...activators.get(id), ...inhibitors.get(id)]);
        outputs.set(id, new Set());
    }
    for (const id of order) {
        for (const i of inputs.get(id)) outputs.get(i).add(id);
    }

    for (const v of variables) {
        const id = toInt(v.Id);
        const [lo, hi] = range.get(id);
        const formula = typeof v.Formula === 'string' ? v.Formula.trim() : '';
        if (formula !== '') {
            exprs.set(id, parseFormula(formula));
            isDefaultExpr.set(id, false);
        } else {
            const pos = activators.get(id), neg = inhibitors.get(id);
            let ast;
            if (pos.length === 0) ast = ['const', lo];
            else if (neg.length === 0) ast = ['min', ['const', hi], ['ave', pos.map((p) => ['var', p])]];
            else {
                ast = ['max', ['const', lo], ['minus',
                    ['ave', pos.map((p) => ['var', p])],
                    ['ave', neg.map((p) => ['var', p])]]];
            }
            exprs.set(id, ast);
            isDefaultExpr.set(id, true);
        }
    }

    /** @type {Map<string,number>} */ const byName = new Map();
    /** @type {Map<string,number>} */ const byKey = new Map();
    for (const id of order) {
        const n = names.get(id);
        byName.set(n, id);
        byKey.set(n.toLowerCase(), id);
    }

    return {
        name: String(doc.Model.Name ?? 'model'),
        order,
        names,
        range,
        exprs,
        effectiveExprs: exprs,
        isDefaultExpr,
        inputs,
        activators,
        inhibitors,
        outputs,
        byName,
        byKey,
        annotation: new Map(),
    };
}

/**
 * Return a copy of `model` in which the given variables have been knocked out,
 * i.e. their target function replaced by a constant — exactly `QN.ko` in
 * `BioCheckAnalyzer/QN.fs`.
 * @param {NetworkModel} model
 * @param {Map<number,number>|Iterable<[number,number]>} ko
 * @returns {NetworkModel}
 */
export function withKnockouts(model, ko) {
    const entries = ko instanceof Map ? [...ko.entries()] : [...ko];
    if (entries.length === 0) return model;
    const effectiveExprs = new Map(model.exprs);
    for (const [id, level] of entries) {
        if (!model.range.has(id)) throw new Error(`cannot perturb unknown variable id ${id}`);
        const [lo, hi] = model.range.get(id);
        if (!Number.isInteger(level) || level < lo || level > hi) {
            throw new Error(`level ${level} for ${model.names.get(id)} is outside its range [${lo}, ${hi}]`);
        }
        effectiveExprs.set(id, ['const', level]);
    }
    return { ...model, effectiveExprs };
}

/**
 * Human-readable target function of a node, with variable ids replaced by names.
 * @param {NetworkModel} model
 * @param {number} id
 * @param {{effective?: boolean}} [opts]
 * @returns {string}
 */
export function formulaOf(model, id, opts = {}) {
    const ast = (opts.effective ? model.effectiveExprs : model.exprs).get(id);
    if (ast === undefined) return '';
    return formatFormula(ast, (v) => model.names.get(v) ?? `var(${v})`);
}

/**
 * Pretty multi-line dump of a model, used by the `describe` tool.
 * @param {NetworkModel} model
 * @returns {string}
 */
export function summarizeModel(model) {
    const nodeCount = model.order.length;
    let edges = 0, defaults = 0;
    for (const id of model.order) {
        edges += model.activators.get(id).length + model.inhibitors.get(id).length;
        if (model.isDefaultExpr.get(id)) defaults += 1;
    }
    const ranges = new Map();
    for (const id of model.order) {
        const k = model.range.get(id).join('..');
        ranges.set(k, (ranges.get(k) ?? 0) + 1);
    }
    const rangeText = [...ranges.entries()].sort().map(([k, n]) => `${n}x [${k}]`).join(', ');
    return `Model "${model.name}": ${nodeCount} nodes, ${edges} relationships, `
        + `${defaults} node(s) using a synthesised default target function. Ranges: ${rangeText}.`;
}
