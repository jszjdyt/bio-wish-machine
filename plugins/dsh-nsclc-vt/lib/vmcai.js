/**
 * The VMCAI stability engine — a faithful port of
 * `src/BioCheckAnalyzer/GenLemmas.fs` (`stabilize_lazy`) together with
 * `src/BioCheckAnalyzer/FNewLemmas.fs` (`tighten_fast`/`tighten_internal`) and
 * `src/BioCheckAnalyzer/Stabilize.fs` (`check_stability_lazy`), the algorithm of
 * Cook, Fisher, Krepska & Piterman, *Proving Stabilization of Biological
 * Models*, VMCAI 2011.
 *
 * What the engine does
 * --------------------
 * A BMA model is a synchronous multi-valued network: every variable takes an
 * integer in its own range, and one tick replaces the whole state vector by its
 * target functions. VMCAI asks whether such a network has a *unique attractor*
 * and, if it does, what interval each variable is confined to once the
 * transient has been eliminated.
 *
 * The analysis is a chaotic interval iteration. Bounds start at each variable's
 * full range (constants start at their value, which is a fixpoint by
 * construction). A frontier holds the variables whose bounds may still tighten.
 * One step picks the *largest* variable id in the frontier, re-evaluates its
 * target function over the current bounds of its inputs, and — if that narrows
 * either endpoint — records the tighter interval and pushes the variable's
 * dependents back onto the frontier. Terminating with the frontier empty means
 * no bound can tighten any further: either every variable is pinned to a single
 * value (the network has a unique fixpoint — "stabilizing") or some interval
 * survived ("not stabilizing", where BMA then hunts for a counter-example).
 *
 * For monotone inputs `tighten_fast` evaluates the function at the endpoint
 * that attains the extremum instead of at the whole box; for non-monotone
 * inputs it enumerates the input's values. On the NSCLC model the largest
 * enumeration is 4 combinations, so a full screen is cheap.
 *
 * Reported numbers
 * ----------------
 * `iterations` is BMA's `Time` field: the number of frontier pops before the
 * frontier drained. It is reproduced exactly because the frontier is drained in
 * descending variable-id order.
 *
 * @module dsh-nsclc-vt/vmcai
 */
import { evalExpr, freeVars, isIncreasing, isDecreasing } from './expr.js';

/** Hard stop so a pathological model cannot hang the tool call. */
const DEFAULT_MAX_ITERATIONS = 200000;

/**
 * `Expr.is_a_const`: an expression is constant when it mentions no variables,
 * or is a bare reference to a variable whose range is a single point.
 * @param {Map<number,[number,number]>} range
 * @param {Array} expr
 * @returns {boolean}
 */
export function isConstantExpr(range, expr) {
    if (expr[0] === 'var') {
        const [lo, hi] = range.get(expr[1]);
        return lo === hi;
    }
    return freeVars(expr).size === 0;
}

/** Cartesian product of `[id, lo, hi]` triples, in F# `all_inputs` spirit. */
function* allInputs(items) {
    const n = items.length;
    if (n === 0) { yield new Map(); return; }
    const idx = items.map((it) => it[1]);
    for (;;) {
        const env = new Map();
        for (let k = 0; k < n; k += 1) env.set(items[k][0], idx[k]);
        yield env;
        let k = n - 1;
        while (k >= 0 && idx[k] === items[k][2]) { idx[k] = items[k][1]; k -= 1; }
        if (k < 0) return;
        idx[k] += 1;
    }
}

/**
 * `FNewLemmas.tighten_internal`: fold the target function over every input
 * combination, seeded with the swapped current interval exactly as the original
 * does, then intersect the result with the current interval on each side.
 */
function tightenInternal(node, range, expr, lower, upper, inputBounds) {
    let loAcc = upper, hiAcc = lower;
    for (const env of allInputs(inputBounds)) {
        const v = evalExpr(node, range, expr, env);
        if (v < loAcc) loAcc = v;
        if (v > hiAcc) hiAcc = v;
    }
    return [loAcc, hiAcc];
}

/** `FNewLemmas.find_max` — pin each input to the endpoint that maximises `expr`. */
function findMax(expr, varId, bounds) {
    const [lo, hi] = bounds.get(varId);
    if (isIncreasing(expr, varId)) return [varId, hi, hi];
    if (isDecreasing(expr, varId)) return [varId, lo, lo];
    return [varId, lo, hi];
}

/** `FNewLemmas.find_min` — pin each input to the endpoint that minimises `expr`. */
function findMin(expr, varId, bounds) {
    const [lo, hi] = bounds.get(varId);
    if (isIncreasing(expr, varId)) return [varId, lo, lo];
    if (isDecreasing(expr, varId)) return [varId, hi, hi];
    return [varId, lo, hi];
}

/**
 * `FNewLemmas.tighten_fast`.
 * @returns {[number, number]} the tightened interval (may equal the input)
 */
export function tightenFast(node, range, expr, lower, upper, inputs, bounds) {
    const minBounds = inputs.map((i) => findMin(expr, i, bounds));
    const maxBounds = inputs.map((i) => findMax(expr, i, bounds));
    const [newLower] = tightenInternal(node, range, expr, lower, upper, minBounds);
    const [, newUpper] = tightenInternal(node, range, expr, lower, upper, maxBounds);
    return [newLower, newUpper];
}

/**
 * @typedef {object} StabilityResult
 * @property {Map<number,[number,number]>} bounds  final interval per variable id
 * @property {number} iterations                   BMA's `Time`: frontier pops
 * @property {boolean} stabilizing                 true when every interval is a point
 * @property {boolean} converged                   false when the iteration cap was hit
 * @property {Array<{tick:number, bounds:Map<number,[number,number]>}>} [history]
 */

/**
 * Run the VMCAI stability analysis on a (possibly knocked-out) model.
 * @param {import('./model.js').NetworkModel} model
 * @param {{maxIterations?: number, recordHistory?: boolean, constantsAsFixpoints?: boolean}} [opts]
 * @returns {StabilityResult}
 */
export function analyzeStability(model, opts = {}) {
    const maxIterations = opts.maxIterations ?? DEFAULT_MAX_ITERATIONS;
    const recordHistory = opts.recordHistory === true;
    const exprs = model.effectiveExprs;
    const range = model.range;
    const inputs = model.inputs;
    const outputs = model.outputs;
    const order = model.order;

    /** @type {Map<number,[number,number]>} */
    const bounds = new Map();
    const frontier = new Set();
    for (const id of order) {
        if (isConstantExpr(range, exprs.get(id))) {
            const c = evalExpr(id, range, exprs.get(id), new Map());
            bounds.set(id, [c, c]);
        } else {
            const [lo, hi] = range.get(id);
            bounds.set(id, [lo, hi]);
            frontier.add(id);
        }
    }

    /** @type {Array<{tick:number, bounds:Map<number,[number,number]>}>|undefined} */
    const history = recordHistory ? [] : undefined;
    const snapshot = () => new Map(bounds);

    if (frontier.size === 0) {
        if (history) history.push({ tick: 0, bounds: snapshot() });
        return {
            bounds: snapshot(),
            iterations: 0,
            stabilizing: true,
            converged: true,
            ...(history ? { history } : {}),
        };
    }

    let pops = 0;
    let converged = false;
    while (pops < maxIterations) {
        if (frontier.size === 0) { converged = true; break; }
        let node = -Infinity;
        for (const f of frontier) if (f > node) node = f;
        frontier.delete(node);
        const [lower, upper] = bounds.get(node);
        const [nl, nu] = tightenFast(node, range, exprs.get(node), lower, upper, inputs.get(node), bounds);
        if (nu < upper || nl > lower) {
            bounds.set(node, [nl, nu]);
            for (const o of outputs.get(node)) frontier.add(o);
        }
        if (history) history.push({ tick: pops, bounds: snapshot() });
        pops += 1;
    }
    if (frontier.size === 0) {
        converged = true;
        if (history) history.push({ tick: pops, bounds: snapshot() });
    }

    let stabilizing = true;
    for (const [, [lo, hi]] of bounds) if (lo !== hi) { stabilizing = false; break; }
    return {
        bounds: snapshot(),
        iterations: pops,
        stabilizing,
        converged,
        ...(history ? { history } : {}),
    };
}

/**
 * The scalar readouts NANSEN computes from a VMCAI tick: `mean` is the interval
 * midpoint and `uncertainty` is the interval width normalised by the variable's
 * full range (`(hi - lo) / (range_to - range_from)`); see `process_results` in
 * NANSEN's `combo_funcs.r`.
 * @param {import('./model.js').NetworkModel} model
 * @param {Map<number,[number,number]>} bounds
 * @param {number} id
 * @returns {{lo:number, hi:number, mean:number, uncertainty:number, width:number}}
 */
export function readout(model, bounds, id) {
    const [lo, hi] = bounds.get(id);
    const [rlo, rhi] = model.range.get(id);
    const span = rhi - rlo;
    return {
        lo,
        hi,
        mean: (lo + hi) / 2,
        uncertainty: span === 0 ? 0 : (hi - lo) / span,
        width: hi - lo,
    };
}
