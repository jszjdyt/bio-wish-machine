/**
 * End-to-end validation against the manuscript's own published output.
 *
 * `results/combination_screen/processed_combination/additive_scores.csv` in
 * https://github.com/jfisher-lab/NSCLC is the paper's combination-benefit
 * matrix: rows are node pairs, columns are the seven screened backgrounds, and
 * the value is
 *
 *     additive = (Proliferation - Death | pair) - min over the two singles
 *
 * This script recomputes that matrix with this plugin's engine and prints the
 * agreement. Run it with
 *
 *     node test/validate-published.mjs /path/to/additive_scores.csv
 *
 * The reference CSV is not redistributed with the plugin (the paper's results
 * directory is ~200 MB); fetch the single file with
 *
 *     curl -L -o additive_scores.csv https://raw.githubusercontent.com/jfisher-lab/NSCLC/main/\
 *       results/combination_screen/processed_combination/additive_scores.csv
 *
 * @module dsh-nsclc-vt/test/validate-published
 */
import { readFileSync } from 'node:fs';
import { createAnalyst } from '../lib/analysis.js';
import { enumerateJobs, runJobs, scoreScreen } from '../lib/screen.js';
import { parseCsvObjects } from '../lib/csv.js';

const referencePath = process.argv[2];
if (referencePath === undefined) {
    console.error('usage: node test/validate-published.mjs <additive_scores.csv>');
    process.exit(2);
}

const analyst = createAnalyst();
const backgrounds = analyst.tables.publishedBackgrounds.map((b) => b.name);
console.log(`screening ${analyst.model.order.length} nodes over ${backgrounds.length} backgrounds: ${backgrounds.join(', ')}`);

const t0 = Date.now();
const { jobs } = enumerateJobs(analyst, { backgrounds, mode: 'both', levels: [0] });
console.log(`jobs: ${jobs.length}`);
const results = await runJobs(analyst, jobs, { inProcessLimit: 200 });
console.log(`solved in ${((Date.now() - t0) / 1000).toFixed(1)}s`);

const scored = scoreScreen(analyst, results);

// index: sorted "A+B" -> { background -> additive }
// The reference CSV labels rows by the two node names only, because the paper's
// screen only ever applies knockouts (level 0), so the level carries no information.
const mine = new Map();
for (const row of scored.rows) {
    if (row.additive === undefined) continue;
    const key = row.combination.slice().sort().join('+');
    let entry = mine.get(key);
    if (entry === undefined) { entry = new Map(); mine.set(key, entry); }
    entry.set(row.background.toLowerCase(), row.additive);
}

const ref = parseCsvObjects(readFileSync(referencePath, 'utf8'));
const header = Object.keys(ref[0] ?? {});
const labelKey = header[0];

let compared = 0, exact = 0, close = 0, missing = 0;
const failures = [];
const mismatches = [];
for (const r of ref) {
    const label = r[labelKey];
    if (label === undefined || label === '') continue;
    const key = label.split('+').sort().join('+');
    const entry = mine.get(key);
    if (entry === undefined) { missing += 1; continue; }
    for (const column of header.slice(1)) {
        const expected = Number(r[column]);
        if (!Number.isFinite(expected)) continue;
        const got = entry.get(column.toLowerCase());
        if (got === undefined) { missing += 1; continue; }
        compared += 1;
        if (Math.abs(got - expected) < 1e-9) exact += 1;
        else if (Math.abs(got - expected) < 1e-6) close += 1;
        else {
            failures.push({ label, column, expected, got });
            if (mismatches.length < 10) mismatches.push(`  ${label} / ${column}: got ${got}, expected ${expected}`);
        }
    }
}

console.log(`\ncompared ${compared} published cells: ${exact} exact, ${close} within 1e-6, ${failures.length} mismatched, ${missing} not covered`);
if (mismatches.length > 0) console.log(mismatches.join('\n'));
process.exit(failures.length === 0 ? 0 : 1);
