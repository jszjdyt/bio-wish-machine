/**
 * Text projections shared by the tools. Tool results are rendered to plain text
 * for the transcript and returned structurally for Code Mode callers, so the
 * numbering here has to read well as both.
 * @module dsh-nsclc-vt/format
 */

/** Format a level interval as `2` or `[2,3]`. */
export function level(lo, hi) {
    return lo === hi ? String(lo) : `[${lo},${hi}]`;
}

/** Signed number with a fixed precision and an explicit plus sign. */
export function signed(x, digits = 2) {
    if (x === undefined || x === null || Number.isNaN(x)) return 'n/a';
    const v = x.toFixed(digits);
    return x > 0 ? `+${v}` : v;
}

/** Fixed-precision number, or `n/a`. */
export function num(x, digits = 2) {
    if (x === undefined || x === null || Number.isNaN(x)) return 'n/a';
    return Number(x).toFixed(digits);
}

/** Right-pad to a width, never truncating. */
export function pad(s, width) {
    const str = String(s);
    return str.length >= width ? str : str + ' '.repeat(width - str.length);
}

/** Left-pad to a width. */
export function padStart(s, width) {
    const str = String(s);
    return str.length >= width ? str : ' '.repeat(width - str.length) + str;
}

/**
 * Render a fixed-width text table.
 * @param {string[]} header
 * @param {Array<Array<string|number>>} rows
 * @param {{align?: Array<'l'|'r'>}} [opts]
 */
export function table(header, rows, opts = {}) {
    const align = opts.align ?? header.map(() => 'l');
    const cells = [header, ...rows].map((r) => r.map((c) => String(c)));
    const widths = header.map((_, i) => Math.max(...cells.map((r) => (r[i] ?? '').length)));
    const line = (r) => r.map((c, i) => (align[i] === 'r' ? padStart(c ?? '', widths[i]) : pad(c ?? '', widths[i]))).join('  ').trimEnd();
    const out = [line(cells[0]), widths.map((w, i) => (align[i] === 'r' ? '-'.repeat(w) : '-'.repeat(w))).join('  ').trimEnd()];
    for (const r of cells.slice(1)) out.push(line(r));
    return out.join('\n');
}

/**
 * Human-readable description of a genotype/intervention set.
 * @param {{name?: string, mutations: Array<{name:string, level:number}>}} background
 * @param {Array<{name:string, level:number}>} interventions
 */
export function describeCase(background, interventions) {
    const bg = background.name !== undefined
        ? `${background.name} (${background.mutations.map((m) => `${m.name}=${m.level}`).join(', ') || 'no mutations'})`
        : (background.mutations.map((m) => `${m.name}=${m.level}`).join(', ') || 'unperturbed model');
    const pert = interventions.length === 0
        ? 'none'
        : interventions.map((p) => `${p.name}=${p.level}`).join(', ');
    return { background: bg, interventions: pert };
}

/**
 * `mean(Proliferation) - mean(Death)` — the manuscript's `Difference` phenotype.
 * Lower is therapeutically better.
 */
export const DIFFERENCE_NOTE = 'Difference = mean(Proliferation) - mean(Death); lower means less proliferation and/or more death.';

/**
 * The comparison quantity, described for what it is: the **predicted difference
 * versus the best single intervention in the same in-silico model**. The
 * manuscript's `additive` name is kept only as the field name for backwards
 * compatibility; it is not a measured synergy call.
 */
export const ADDITIVE_NOTE = 'Δ vs best single intervention = Difference(pair) - min(Difference(single A), Difference(single B)) in the same in-silico model; negative = the pair beats its better single component, positive = it does not. This is a model prediction, not experimental synergy evidence.';

/** Provenance label for every model-derived number this plugin reports. */
export const MODEL_PREDICTION_NOTE = 'Model prediction (mechanistic in-silico simulation), not an experimental result; validate candidate combinations in vitro or in vivo before treating them as findings.';
