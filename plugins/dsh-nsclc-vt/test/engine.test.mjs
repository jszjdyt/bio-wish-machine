/**
 * Engine tests: the VMCAI port, the model importer, and the expression layer,
 * checked against results published with the manuscript.
 *
 *     node --test test/
 *
 * The golden fixture in `fixtures/published-cases.json` was sampled from
 * `results/combination_screen/node_results.csv` of
 * https://github.com/jfisher-lab/NSCLC — the authors' own BioCheckConsole
 * output — so a pass here means this port agrees with the engine that produced
 * the paper, including its refinement-step counter.
 *
 * @module dsh-nsclc-vt/test/engine.test
 */
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

import { loadModel, withKnockouts, summarizeModel, formulaOf } from '../lib/model.js';
import { analyzeStability, readout } from '../lib/vmcai.js';
import { parseFormula, evalExpr, isIncreasing, isDecreasing, freeVars } from '../lib/expr.js';
import { createAnalyst } from '../lib/analysis.js';
import { enumerateJobs, runJobsInline, scoreScreen } from '../lib/screen.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const fixture = JSON.parse(readFileSync(join(HERE, 'fixtures', 'published-cases.json'), 'utf8'));
const analyst = createAnalyst();
const model = analyst.model;

test('the importer reports the model the manuscript describes', () => {
    assert.equal(model.order.length, 162);
    assert.equal(model.name, 'NSCLC-master');
    const relationships = model.order.reduce((n, id) => n + model.activators.get(id).length + model.inhibitors.get(id).length, 0);
    assert.equal(relationships, 360);
    assert.match(summarizeModel(model), /162 nodes, 360 relationships/);
});

test('every published case reproduces exactly, including the Time counter', () => {
    const failures = [];
    for (const c of fixture.cases) {
        const { analysis } = analyst.solve({ background: c.background, interventions: c.interventions });
        const p = readout(model, analysis.bounds, model.byName.get('Proliferation'));
        const d = readout(model, analysis.bounds, model.byName.get('Death'));
        const got = { lo: p.lo, hi: p.hi };
        const gotD = { lo: d.lo, hi: d.hi };
        if (got.lo !== c.proliferation.lo || got.hi !== c.proliferation.hi
            || gotD.lo !== c.death.lo || gotD.hi !== c.death.hi
            || analysis.iterations !== c.iterations) {
            failures.push(`${c.background} [${c.interventions.join(', ')}]: got P${JSON.stringify(got)} D${JSON.stringify(gotD)} t=${analysis.iterations}, `
                + `expected P${JSON.stringify(c.proliferation)} D${JSON.stringify(c.death)} t=${c.iterations}`);
        }
    }
    assert.deepEqual(failures, [], `${failures.length} of ${fixture.cases.length} published cases diverged:\n${failures.slice(0, 8).join('\n')}`);
});

test('knockout replaces the target function with a constant, and the range is left intact', () => {
    const p53 = model.byName.get('p53');
    const knocked = withKnockouts(model, new Map([[p53, 0]]));
    assert.deepEqual(knocked.effectiveExprs.get(p53), ['const', 0]);
    assert.deepEqual(model.range.get(p53), knocked.range.get(p53), 'the declared range must survive a knockout');
    // A downstream reader still rescales the pinned constant by the declared range,
    // which is why `-ko` and a genuinely constant node are not interchangeable.
    const downstream = [...model.outputs.get(p53)][0];
    assert.ok(model.inputs.get(downstream).includes(p53));
});

test('interventions take precedence over the background, as in the NANSEN pipeline', () => {
    const ras = model.byName.get('RAS');
    const a = analyzeStability(withKnockouts(model, new Map([[ras, 3]])));
    const b = analyzeStability(withKnockouts(model, new Map([[ras, 3], [ras, 0]])));
    assert.notDeepEqual(a.bounds.get(ras), b.bounds.get(ras));
    assert.deepEqual(b.bounds.get(ras), [0, 0]);
});

test('default target functions are synthesised from the edges exactly as Marshal.fs does', () => {
    // Proliferation has one activator (Mitosis) and no inhibitors -> min(RangeTo, avg(activators)).
    const pro = model.byName.get('Proliferation');
    assert.equal(analyst.isDefaultFunction(pro), true);
    assert.deepEqual(model.exprs.get(pro), ['min', ['const', 3], ['ave', [['var', model.byName.get('Mitosis')]]]]);
    // A node with no edges at all becomes a constant at its lower bound.
    const orphan = model.order.find((id) => model.inputs.get(id).length === 0 && analyst.isDefaultFunction(id));
    assert.deepEqual(model.exprs.get(orphan), ['const', model.range.get(orphan)[0]]);
    // A node with both kinds of edge uses max(RangeFrom, avg(pos) - avg(neg)).
    const both = model.order.find((id) => analyst.isDefaultFunction(id)
        && model.activators.get(id).length > 0 && model.inhibitors.get(id).length > 0);
    assert.equal(model.exprs.get(both)[0], 'max');
});

test('evalExpr rescales inputs into the reading node range, rounds half-up and clamps', () => {
    const range = new Map([[1, [0, 3]], [2, [0, 1]], [3, [0, 3]]]);
    // var(2) lives on 0..1; reading it from a 0..3 node stretches it 3x.
    assert.equal(evalExpr(3, range, parseFormula('var(2)'), new Map([[2, 1]])), 3);
    assert.equal(evalExpr(3, range, parseFormula('var(2)'), new Map([[2, 0]])), 0);
    // Half-up rounding: 2.5 -> 3, 2.4 -> 2.
    assert.equal(evalExpr(3, range, ['div', ['const', 5], ['const', 2]], new Map()), 3);
    assert.equal(evalExpr(3, range, ['div', ['const', 12], ['const', 5]], new Map()), 2);
    // Clamping into the node's own range.
    assert.equal(evalExpr(2, range, ['const', 3], new Map()), 1);
    // A missing input evaluates to the node's lower bound.
    assert.equal(evalExpr(3, range, ['plus', ['const', 2], ['var', 2]], new Map()), 2);
});

test('monotonicity analysis matches the sign reasoning the engine depends on', () => {
    const e = parseFormula('3 - ceil(avg(var(184), var(387)))');
    assert.equal(isDecreasing(e, 184), true);
    assert.equal(isIncreasing(e, 184), false);
    const f = parseFormula('min(1, var(11))');
    assert.equal(isIncreasing(f, 11), true);
    assert.equal(isDecreasing(f, 11), false);
    assert.deepEqual([...freeVars(f)], [11]);
});

test('formulaOf substitutes node names for ids', () => {
    const death = model.byName.get('Death');
    const text = formulaOf(model, death);
    assert.match(text, /Apoptosis/);
    assert.match(text, /Necrosis/);
    assert.ok(!text.includes('var('), 'ids should be replaced by names');
});

test('a screen scores additivity as Difference(pair) - min(Difference of the singles)', () => {
    const { jobs } = enumerateJobs(analyst, {
        backgrounds: ['h1993'],
        mode: 'both',
        levels: [0],
        candidates: ['14_3_3', '53BP1', 'ATM'],
    });
    const scored = scoreScreen(analyst, runJobsInline(analyst, jobs));
    const row = scored.rows.find((r) => r.combination.slice().sort().join('+') === '14_3_3+53BP1');
    assert.ok(row, 'the pair should have been scored');
    assert.equal(row.difference, 1);
    assert.equal(row.singleADifference, 2);
    assert.equal(row.singleBDifference, 1);
    assert.equal(row.additive, 1 - Math.min(2, 1));
});

test('phenotype nodes are never perturbed by a screen', () => {
    const { candidateIds } = enumerateJobs(analyst, { backgrounds: ['healthy'], mode: 'mono' });
    for (const id of analyst.phenotypeIds) assert.ok(!candidateIds.includes(id));
});

test('an unstable model is reported as such rather than silently pinned', () => {
    const { analysis } = analyst.solve({ background: 'sw1573' });
    assert.equal(analysis.stabilizing, false, 'the sw1573 background leaves intervals in the manuscript data too');
    const p = readout(model, analysis.bounds, model.byName.get('Proliferation'));
    assert.deepEqual([p.lo, p.hi], [2, 3]);
});
