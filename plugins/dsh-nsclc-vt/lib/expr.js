/**
 * Expression layer of the BioModelAnalyzer (BMA) qualitative-network language —
 * a faithful port of `src/BioCheckAnalyzer/Expr.fs` and `ParsecExpr.fs`.
 *
 * The AST is a plain array-tagged tuple so it can be JSON-serialised and
 * compared structurally:
 *
 *   ['const', 3]                 ['var', 184]
 *   ['plus'|'minus'|'times'|'div'|'max'|'min', a, b]
 *   ['ceil'|'floor'|'abs', a]    ['ave', [..]]   ['sum', [..]]
 *
 * The two behaviours that matter for reproducing published results are:
 *
 *  1. `evalExpr` performs *range conversion* on every variable reference, not
 *     a plain lookup: a value `x` of a variable with range `[vLo,vHi]` becomes
 *     `(x - vLo) * (nHi - nLo)/(vHi - vLo) + nLo` in the range of the node being
 *     updated, then `trunc(v + 0.5)` and a clamp into the node's own range.
 *  2. A variable that is *absent* from the environment evaluates to the node's
 *     own lower bound — this is how BMA treats a formula input that is not a
 *     declared relationship input of that node.
 *
 * @module dsh-nsclc-vt/expr
 */

const TOKEN = /\s+|(?<num>\d+)|(?<name>[A-Za-z_][A-Za-z_0-9]*)\s*\(|(?<op>[-+*/(),])/y;

/**
 * Split a formula string into `[kind, text]` tokens.
 * @param {string} src
 * @returns {Array<[string, string]>}
 */
export function tokenize(src) {
    const toks = [];
    TOKEN.lastIndex = 0;
    let i = 0;
    while (i < src.length) {
        TOKEN.lastIndex = i;
        const m = TOKEN.exec(src);
        if (m === null || m.index !== i) {
            throw new SyntaxError(`unexpected character ${JSON.stringify(src[i])} in ${JSON.stringify(src)}`);
        }
        i = TOKEN.lastIndex;
        // For `name` the match includes the trailing "(", which the parser
        // consumes implicitly; only the identifier itself becomes the token.
        const kind = m.groups.num !== undefined ? 'num'
            : m.groups.name !== undefined ? 'name'
                : m.groups.op !== undefined ? 'op' : null;
        if (kind !== null) {
            const text = kind === 'name' ? m.groups.name : m[0].trim();
            toks.push([kind, text]);
        }
    }
    return toks;
}

class Parser {
    /** @param {Array<[string, string]>} toks */
    constructor(toks) {
        this.toks = toks;
        this.i = 0;
    }

    peek() {
        return this.i < this.toks.length ? this.toks[this.i][1] : null;
    }

    take() {
        const t = this.toks[this.i];
        this.i += 1;
        return t;
    }

    want(v) {
        if (this.peek() !== v) {
            throw new SyntaxError(`expected ${JSON.stringify(v)}, found ${JSON.stringify(this.peek())}`);
        }
        this.take();
    }

    parse() {
        const e = this.expr();
        if (this.i !== this.toks.length) throw new SyntaxError('trailing tokens in formula');
        return e;
    }

    expr() {
        let n = this.term();
        while (this.peek() === '+' || this.peek() === '-') {
            const op = this.take()[1];
            n = [op === '+' ? 'plus' : 'minus', n, this.term()];
        }
        return n;
    }

    term() {
        let n = this.unary();
        while (this.peek() === '*' || this.peek() === '/') {
            const op = this.take()[1];
            n = [op === '*' ? 'times' : 'div', n, this.unary()];
        }
        return n;
    }

    unary() {
        if (this.peek() === '-') {
            this.take();
            return ['minus', ['const', 0], this.unary()];
        }
        return this.atom();
    }

    atom() {
        const [kind, text] = this.take();
        if (text === '(') {
            const e = this.expr();
            this.want(')');
            return e;
        }
        if (kind === 'num') return ['const', Number(text)];
        if (kind === 'name') {
            const args = [];
            if (this.peek() !== ')') {
                args.push(this.expr());
                while (this.peek() === ',') {
                    this.take();
                    args.push(this.expr());
                }
            }
            this.want(')');
            switch (text) {
                case 'var': return ['var', args[0][1]];
                case 'avg': return ['ave', args];
                case 'sum': return ['sum', args];
                case 'min':
                case 'max': {
                    let e = args[0];
                    for (const a of args.slice(1)) e = [text, e, a];
                    return e;
                }
                case 'ceil':
                case 'floor':
                case 'abs': return [text, args[0]];
                default: throw new SyntaxError(`unknown function ${JSON.stringify(text)}`);
            }
        }
        throw new SyntaxError(`unexpected token ${JSON.stringify(text)}`);
    }
}

/**
 * Parse a BMA target function into an AST. Embedded newlines are folded away
 * first (the NSCLC model file contains multi-line formulas).
 * @param {string} src
 * @returns {Array} AST
 */
export function parseFormula(src) {
    return new Parser(tokenize(src.replace(/\s+/g, ' ').trim())).parse();
}

/**
 * Free variables of an AST.
 * @param {Array} e
 * @param {Set<number>} [out]
 * @returns {Set<number>}
 */
export function freeVars(e, out = new Set()) {
    switch (e[0]) {
        case 'var': out.add(e[1]); break;
        case 'plus': case 'minus': case 'times': case 'div': case 'max': case 'min':
            freeVars(e[1], out); freeVars(e[2], out); break;
        case 'ceil': case 'floor': case 'abs': freeVars(e[1], out); break;
        case 'ave': case 'sum': for (const a of e[1]) freeVars(a, out); break;
        default: break;
    }
    return out;
}

/**
 * Render an AST back to a readable formula, substituting node names.
 *
 * Parenthesisation follows real precedence — `a - (b - c)` keeps its brackets —
 * because these strings are shown to the user as the model's mechanism, and a
 * silently re-associated rendering would misdescribe the target function.
 *
 * @param {Array} e
 * @param {(id: number) => string} nameOf
 * @returns {string}
 */
export function formatFormula(e, nameOf) {
    /** Precedence: additive binds loosest, then multiplicative, then unary/atomic. */
    const PREC = { plus: 1, minus: 1, times: 2, div: 2, neg: 3, atom: 4 };
    const go = (x) => {
        switch (x[0]) {
            case 'const': return { text: String(x[1]), prec: PREC.atom };
            case 'var': return { text: nameOf(x[1]), prec: PREC.atom };
            case 'plus': return binary(x, ' + ', PREC.plus);
            case 'minus': return binary(x, ' - ', PREC.minus);
            case 'times': return binary(x, ' * ', PREC.times);
            case 'div': return binary(x, ' / ', PREC.div);
            case 'max': return call('max', [x[1], x[2]]);
            case 'min': return call('min', [x[1], x[2]]);
            case 'ceil': return call('ceil', [x[1]]);
            case 'floor': return call('floor', [x[1]]);
            case 'abs': return call('abs', [x[1]]);
            case 'ave': return call('avg', x[1]);
            case 'sum': return call('sum', x[1]);
            default: throw new SyntaxError(`unknown node ${x[0]}`);
        }
    };
    const binary = (x, op, prec) => {
        const left = go(x[1]);
        const right = go(x[2]);
        // A same-precedence right operand of `-` or `/` needs brackets; so does a
        // looser operand on either side.
        const needRight = right.prec < prec || (right.prec === prec && (op === ' - ' || op === ' / '));
        const leftText = left.prec < prec ? `(${left.text})` : left.text;
        const rightText = needRight ? `(${right.text})` : right.text;
        return { text: `${leftText}${op}${rightText}`, prec };
    };
    const call = (name, args) => ({ text: `${name}(${args.map((a) => go(a).text).join(', ')})`, prec: PREC.atom });
    return go(e).text;
}

/**
 * Evaluate a formula exactly as `Expr.eval_expr_int` does.
 * @param {number} node variable id being updated (controls scaling/clamping)
 * @param {Map<number,[number,number]>} range variable id -> [lo, hi]
 * @param {Array} e AST
 * @param {Map<number,number>|Record<number,number>} env input bindings
 * @returns {number} integer in the node's range
 */
export function evalExpr(node, range, e, env) {
    const [nodeMin, nodeMax] = range.get(node);
    const get = typeof env.get === 'function' ? (k) => env.get(k) : (k) => env[k];

    const ev = (x) => {
        switch (x[0]) {
            case 'var': {
                const v = x[1];
                const [vMin, vMax] = range.get(v);
                if (vMin === vMax) return vMin;
                const scale = (nodeMax - nodeMin) / (vMax - vMin);
                const raw = get(v);
                if (raw === undefined) return nodeMin;
                return (raw - vMin) * scale + nodeMin;
            }
            case 'const': return x[1];
            case 'plus': return ev(x[1]) + ev(x[2]);
            case 'minus': return ev(x[1]) - ev(x[2]);
            case 'times': return ev(x[1]) * ev(x[2]);
            case 'div': return ev(x[1]) / ev(x[2]);
            case 'max': return Math.max(ev(x[1]), ev(x[2]));
            case 'min': return Math.min(ev(x[1]), ev(x[2]));
            case 'ceil': return Math.ceil(ev(x[1]));
            case 'floor': return Math.floor(ev(x[1]));
            case 'abs': return Math.abs(ev(x[1]));
            case 'ave': {
                if (x[1].length === 0) return 0;
                let total = 0;
                for (const a of x[1]) total += ev(a);
                return total / x[1].length;
            }
            case 'sum': {
                let total = 0;
                for (const a of x[1]) total += ev(a);
                return total;
            }
            default: throw new SyntaxError(`unknown node ${x[0]}`);
        }
    };

    // `int (convert v)` with convert = (+ 0.5): truncate toward zero after the shift.
    const res = Math.trunc(ev(e) + 0.5);
    if (res < nodeMin) return nodeMin;
    if (res > nodeMax) return nodeMax;
    return res;
}

const POS = 'Pos', NEG = 'Neg', ZERO = 'Zero', UNK = 'Unk';
const cnt = (s, w) => (s === w ? 1 : 0);

const signPlus = (a, b) => {
    if (a === POS || b === POS) return POS;
    if (a === NEG && b === NEG) return NEG;
    if (a === UNK || b === UNK) return UNK;
    return ZERO;
};
const signMinus = (a, b) => {
    if (a === POS || b === NEG) return POS;
    if (a === NEG && b === POS) return NEG;
    if (a === UNK || b === UNK) return UNK;
    return ZERO;
};
const signTimes = (a, b) => {
    if (a === ZERO || b === ZERO) return ZERO;
    if (a === UNK || b === UNK) return UNK;
    if (a === b) return POS;
    if ((a === POS && b === NEG) || (a === NEG && b === POS)) return NEG;
    return UNK;
};
const signMax = (a, b) => {
    const p = cnt(a, POS) + cnt(b, POS), n = cnt(a, NEG) + cnt(b, NEG), z = cnt(a, ZERO) + cnt(b, ZERO);
    if (p > 0) return POS;
    if (z > 0) return ZERO;
    if (n === 2) return NEG;
    return UNK;
};
const signMin = (a, b) => {
    const p = cnt(a, POS) + cnt(b, POS), n = cnt(a, NEG) + cnt(b, NEG), z = cnt(a, ZERO) + cnt(b, ZERO);
    if (n > 0) return NEG;
    if (z > 0) return ZERO;
    if (p === 2) return POS;
    return UNK;
};

const signCache = new WeakMap();
/**
 * Sign abstraction of an expression — `Expr.sign`.
 * @param {Array} e
 * @returns {'Pos'|'Neg'|'Zero'|'Unk'}
 */
export function signOf(e) {
    const hit = signCache.get(e);
    if (hit !== undefined) return hit;
    let s;
    switch (e[0]) {
        case 'const': s = e[1] > 0 ? POS : e[1] < 0 ? NEG : ZERO; break;
        case 'var': s = POS; break;
        case 'abs': s = signOf(e[1]) === ZERO ? ZERO : POS; break;
        case 'plus': s = signPlus(signOf(e[1]), signOf(e[2])); break;
        case 'minus': s = signMinus(signOf(e[1]), signOf(e[2])); break;
        case 'times': case 'div': s = signTimes(signOf(e[1]), signOf(e[2])); break;
        case 'max': s = signMax(signOf(e[1]), signOf(e[2])); break;
        case 'min': s = signMin(signOf(e[1]), signOf(e[2])); break;
        case 'ceil': case 'floor': s = signOf(e[1]); break;
        case 'ave': case 'sum': {
            let acc = ZERO;
            for (const a of e[1]) acc = signPlus(acc, signOf(a));
            s = acc;
            break;
        }
        default: throw new SyntaxError(`unknown node ${e[0]}`);
    }
    signCache.set(e, s);
    return s;
}

const incCache = new WeakMap();
const decCache = new WeakMap();

/**
 * Monotonicity of an expression in one variable — `Expr.is_increasing`.
 * @param {Array} e
 * @param {number} v
 * @returns {boolean}
 */
export function isIncreasing(e, v) {
    let byVar = incCache.get(e);
    if (byVar === undefined) { byVar = new Map(); incCache.set(e, byVar); }
    const hit = byVar.get(v);
    if (hit !== undefined) return hit;
    let r;
    switch (e[0]) {
        case 'const': case 'var': r = true; break;
        case 'abs': r = signOf(e[1]) === POS ? isIncreasing(e[1], v) : isDecreasing(e[1], v); break;
        case 'plus': case 'min': case 'max': r = isIncreasing(e[1], v) && isIncreasing(e[2], v); break;
        case 'ceil': case 'floor': r = isIncreasing(e[1], v); break;
        case 'minus': r = isIncreasing(e[1], v) && isDecreasing(e[2], v); break;
        case 'times': {
            const s1 = signOf(e[1]), s2 = signOf(e[2]);
            const i1 = s1 === POS ? isIncreasing(e[2], v) : s1 === NEG ? isDecreasing(e[2], v) : false;
            const i2 = s2 === POS ? isIncreasing(e[1], v) : s2 === NEG ? isDecreasing(e[1], v) : false;
            r = i1 && i2;
            break;
        }
        case 'div': {
            const s1 = signOf(e[1]), s2 = signOf(e[2]);
            const i1 = s1 === POS ? isDecreasing(e[2], v) : s1 === NEG ? isIncreasing(e[2], v) : false;
            const i2 = s2 === POS ? isIncreasing(e[1], v) : s2 === NEG ? isDecreasing(e[1], v) : false;
            r = i1 && i2;
            break;
        }
        case 'ave': case 'sum': r = e[1].every((a) => isIncreasing(a, v)); break;
        default: throw new SyntaxError(`unknown node ${e[0]}`);
    }
    byVar.set(v, r);
    return r;
}

/**
 * Monotonicity of an expression in one variable — `Expr.is_decreasing`.
 * @param {Array} e
 * @param {number} v
 * @returns {boolean}
 */
export function isDecreasing(e, v) {
    let byVar = decCache.get(e);
    if (byVar === undefined) { byVar = new Map(); decCache.set(e, byVar); }
    const hit = byVar.get(v);
    if (hit !== undefined) return hit;
    let r;
    switch (e[0]) {
        case 'const': r = true; break;
        case 'var': r = e[1] !== v; break;
        case 'abs': r = signOf(e[1]) === POS ? isDecreasing(e[1], v) : isIncreasing(e[1], v); break;
        case 'plus': case 'min': case 'max': r = isDecreasing(e[1], v) && isDecreasing(e[2], v); break;
        case 'ceil': case 'floor': r = isDecreasing(e[1], v); break;
        case 'minus': r = isDecreasing(e[1], v) && isIncreasing(e[2], v); break;
        case 'times': {
            const s1 = signOf(e[1]), s2 = signOf(e[2]);
            const d1 = s1 === POS ? isDecreasing(e[2], v) : s1 === NEG ? isIncreasing(e[2], v) : false;
            const d2 = s2 === POS ? isDecreasing(e[1], v) : s2 === NEG ? isIncreasing(e[1], v) : false;
            r = d1 && d2;
            break;
        }
        case 'div': {
            const s1 = signOf(e[1]), s2 = signOf(e[2]);
            const d1 = s1 === POS ? isIncreasing(e[2], v) : s1 === NEG ? isDecreasing(e[2], v) : false;
            const d2 = s2 === POS ? isDecreasing(e[1], v) : s2 === NEG ? isIncreasing(e[1], v) : false;
            r = d1 && d2;
            break;
        }
        case 'ave': case 'sum': r = e[1].every((a) => isDecreasing(a, v)); break;
        default: throw new SyntaxError(`unknown node ${e[0]}`);
    }
    byVar.set(v, r);
    return r;
}
