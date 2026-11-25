/**
 * The curated companion tables that ship with the NSCLC virtual-tumor model.
 * Everything here is read-only reference data copied verbatim from
 * https://github.com/jfisher-lab/NSCLC (`model/`), so provenance and citations
 * stay with the numbers.
 * @module dsh-nsclc-vt/data
 */
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { parseCsvObjects } from './csv.js';

const HERE = dirname(fileURLToPath(import.meta.url));
/** Plugin-relative data directory. */
export const DATA_DIR = join(HERE, '..', 'data');
/** The bundled model file, as shipped. */
export const DEFAULT_MODEL_PATH = join(DATA_DIR, 'NSCLC-master.json');

const read = (name) => readFileSync(join(DATA_DIR, name), 'utf8');

let cache;

/**
 * Load (and memoise) every companion table.
 * @returns {{
 *   nodesInfo: Array<{node:string, gene:string, protein:string}>,
 *   targetFunctions: Array<{node:string, target:string, min?:number, max?:number, citation:string, reasoning:string}>,
 *   edges: Array<{from:string, to:string, type:string, tissue:string, experiment:string, citation:string, doi:string, inModel:boolean}>,
 *   publishedBackgrounds: Array<{name:string, mutations:Array<{node:string, level:number}>}>,
 * }}
 */
export function loadTables() {
    if (cache !== undefined) return cache;

    const nodesInfo = parseCsvObjects(read('nodes_info.csv')).map((r) => ({
        node: r['Node name'],
        gene: r['Gene Symbol'],
        protein: r['Protein Name (Full)'],
    }));

    const targetFunctions = parseCsvObjects(read('target_functions_master.csv')).map((r) => ({
        node: r.Node,
        target: r['Target.Function'],
        min: r['Minimum.Value'] === '' ? undefined : Number(r['Minimum.Value']),
        max: r['Maximum.Value'] === '' ? undefined : Number(r['Maximum.Value']),
        citation: r.Citation ?? '',
        reasoning: r.Reasoning ?? '',
    }));

    const edges = parseCsvObjects(read('edges_master.csv')).map((r) => ({
        from: r.From,
        to: r.To,
        type: r.Type,
        tissue: r['Cell.Line.Tissue'] ?? '',
        experiment: r.Experiment ?? '',
        citation: r.Citation ?? '',
        doi: r.DOI ?? '',
        inModel: (r['In.Model'] ?? '').toUpperCase() === 'Y',
    }));

    // `combo_coker_mut_background.csv` lists the mutational backgrounds the
    // paper screened; the file name in the repository is lower-cased for
    // portability while the canonical name keeps its Coker capitalisation.
    const backgroundRows = parseCsvObjects(read('combo_coker_mut_background.csv'));
    /** @type {Map<string, Array<{node:string, level:number}>>} */
    const byName = new Map();
    for (const r of backgroundRows) {
        const name = r.background;
        if (!byName.has(name)) byName.set(name, []);
        byName.get(name).push({ node: r.node, level: Number(r.activity) });
    }
    const publishedBackgrounds = [...byName.entries()].map(([name, mutations]) => ({ name, mutations }));

    cache = { nodesInfo, targetFunctions, edges, publishedBackgrounds };
    return cache;
}

/** Drop a memoised table load; used by tests. */
export function resetTables() { cache = undefined; }
