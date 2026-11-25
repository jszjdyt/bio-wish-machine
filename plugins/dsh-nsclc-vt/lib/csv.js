/**
 * Minimal RFC 4180 CSV reader — enough for the model's companion tables, which
 * use quoted fields containing commas and embedded newlines.
 * @module dsh-nsclc-vt/csv
 */

/**
 * Parse CSV text into rows of strings. A trailing newline does not produce an
 * extra empty row.
 * @param {string} text
 * @param {{delimiter?: string}} [opts]
 * @returns {string[][]}
 */
export function parseCsv(text, opts = {}) {
    const delimiter = opts.delimiter ?? ',';
    const rows = [];
    let field = '';
    let row = [];
    let quoted = false;
    let i = 0;
    const push = () => { row.push(field); field = ''; };
    const endRow = () => { push(); rows.push(row); row = []; };

    while (i < text.length) {
        const c = text[i];
        if (quoted) {
            if (c === '"') {
                if (text[i + 1] === '"') { field += '"'; i += 2; continue; }
                quoted = false; i += 1; continue;
            }
            field += c; i += 1; continue;
        }
        if (c === '"' && field === '') { quoted = true; i += 1; continue; }
        if (c === delimiter) { push(); i += 1; continue; }
        if (c === '\r') { i += 1; continue; }
        if (c === '\n') { endRow(); i += 1; continue; }
        field += c; i += 1;
    }
    if (field !== '' || row.length > 0) endRow();
    return rows;
}

/**
 * Parse CSV text into objects keyed by the header row.
 * @param {string} text
 * @returns {Array<Record<string,string>>}
 */
export function parseCsvObjects(text) {
    const rows = parseCsv(text);
    if (rows.length === 0) return [];
    const header = rows[0].map((h, idx) => (idx === 0 ? h.replace(/^\uFEFF/, '').trim() : h.trim()));
    return rows.slice(1)
        .filter((r) => r.length > 1 || (r.length === 1 && r[0] !== ''))
        .map((r) => {
            /** @type {Record<string,string>} */
            const o = {};
            header.forEach((h, idx) => { o[h] = r[idx] ?? ''; });
            return o;
        });
}
