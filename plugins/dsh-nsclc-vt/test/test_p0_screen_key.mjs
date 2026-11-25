/**
 * P0 #10 回归：组合键必须含干预强度（(background|target@level)），内存与 CSV 同源。
 *
 * 修复前 `lib/screen.js:scoreScreen()` 的聚合键是
 * `row.combination.slice().sort().join('+')`，只保留基因名、丢掉 level，
 * 同一对基因在不同强度下的组合会被合并成一行；`lib/tools/nsclc-screen.js`
 * 的 CSV `additiveBy` 键同样丢 level，不同强度的 additive 互相覆盖。
 * 修复后统一用 `combinationKey()/pairKey()`（双向对称、含 level），
 * 并在写 CSV 后与内存结果做同源一致性校验。
 */
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, mkdtempSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

import { createAnalyst } from '../lib/analysis.js';
import {
    enumerateJobs, runJobsInline, scoreScreen, rankCombinations,
    combinationKey, pairKey,
} from '../lib/screen.js';
import { makeNsclcScreenTool } from '../lib/tools/nsclc-screen.js';

const analyst = createAnalyst();
const idOf = (name) => analyst.resolveNode(name).id;

const single = (background, name, level, difference) => ({
    key: `${background}|${name}:${level}`, background, kind: 'single',
    a: idOf(name), levelA: level, difference, proliferation: 1, death: 1,
});
const pair = (background, nameA, levelA, nameB, levelB, difference) => ({
    key: `${background}|${nameA}:${levelA}+${nameB}:${levelB}`, background, kind: 'pair',
    a: idOf(nameA), levelA, b: idOf(nameB), levelB, difference,
    proliferation: 1, death: 1, uncertainty: 0, stabilizing: true,
});
const baseline = (background, difference) => ({
    key: `${background}|baseline`, background, kind: 'baseline', difference,
    proliferation: 1, death: 1,
});

test('combination key keeps the level and is symmetric', () => {
    assert.equal(combinationKey(['MEK', 'PI3K'], [0, 3]), combinationKey(['PI3K', 'MEK'], [3, 0]));
    assert.notEqual(combinationKey(['MEK', 'PI3K'], [0, 3]), combinationKey(['MEK', 'PI3K'], [3, 0]));
    assert.match(combinationKey(['MEK', 'PI3K'], [0, 3]), /@/);
    assert.equal(pairKey('A549', ['MEK', 'PI3K'], [0, 3]), `A549|${combinationKey(['MEK', 'PI3K'], [0, 3])}`);
});

test('scoreScreen keeps the same gene pair at different levels apart', () => {
    const results = [
        baseline('A549', 10),
        single('A549', 'MEK', 0, 8), single('A549', 'PI3K', 0, 7),
        single('A549', 'MEK', 3, 6), single('A549', 'PI3K', 3, 5),
        pair('A549', 'MEK', 0, 'PI3K', 0, 4),   // additive = 4 - min(8,7) = -3
        pair('A549', 'MEK', 3, 'PI3K', 3, 2),   // additive = 2 - min(6,5) = -3
    ];
    const scored = scoreScreen(analyst, results);
    assert.equal(scored.combinations.length, 2, 'two different levels must stay two combinations');
    const keys = scored.combinations.map((c) => c.key);
    assert.equal(new Set(keys).size, 2);
    assert.ok(keys.includes(combinationKey(['MEK', 'PI3K'], [0, 0])));
    assert.ok(keys.includes(combinationKey(['MEK', 'PI3K'], [3, 3])));
    for (const combo of scored.combinations) {
        assert.equal(combo.levels.length, 2);
        assert.equal(combo.perBackground.length, 1);
    }
});

test('the same pair written in the opposite order aggregates across backgrounds', () => {
    const results = [
        single('A549', 'MEK', 0, 8), single('A549', 'PI3K', 0, 7),
        single('Healthy', 'MEK', 0, 8), single('Healthy', 'PI3K', 0, 7),
        pair('A549', 'MEK', 0, 'PI3K', 0, 4),
        pair('Healthy', 'PI3K', 0, 'MEK', 0, 3),   // reversed orientation
    ];
    const scored = scoreScreen(analyst, results);
    assert.equal(scored.combinations.length, 1);
    assert.equal(scored.combinations[0].perBackground.length, 2);
});

test('CSV additive column uses the level-aware key and matches in-memory results', async () => {
    const dir = mkdtempSync(join(tmpdir(), 'nsclc-p0-csv-'));
    const csvPath = join(dir, 'screen.csv');
    const tool = makeNsclcScreenTool(analyst, { workers: 1 });
    await tool.execute({
        backgrounds: ['healthy'], mode: 'both', candidates: ['MEK', 'PI3K'],
        levels: [0, 3], top: 5, output_csv: csvPath,
    }, {});

    // in-memory expectation, computed from the same enumeration/runner
    const { jobs } = enumerateJobs(analyst, {
        backgrounds: ['healthy'], mode: 'both', candidates: ['MEK', 'PI3K'], levels: [0, 3],
    });
    const expected = scoreScreen(analyst, runJobsInline(analyst, jobs));
    const additiveBy = new Map();
    for (const row of expected.rows) {
        additiveBy.set(pairKey(row.background, row.combination, row.levels), row.additive);
    }
    // 2 candidates × 2 levels => 2×2 cross pairs (same-node pairs are not screened)
    assert.equal(additiveBy.size, 4, 'four level-aware pair keys');

    const lines = readFileSync(csvPath, 'utf8').trim().split('\n');
    const header = lines[0].split(',');
    const col = (name) => header.indexOf(name);
    const pairRows = lines.slice(1).map((line) => line.split(',')).filter((cells) => cells[col('kind')] === 'pair');
    assert.equal(pairRows.length, 4);
    const seen = new Set();
    for (const cells of pairRows) {
        const background = cells[col('background')];
        const nameA = cells[col('a')];
        const nameB = cells[col('b')];
        const levelA = Number(cells[col('level_a')]);
        const levelB = Number(cells[col('level_b')]);
        const key = pairKey(background, [nameA, nameB], [levelA, levelB]);
        assert.ok(!seen.has(key), `duplicate CSV key ${key}`);
        seen.add(key);
        const csvAdditive = Number(cells[col('additive')]);
        const want = additiveBy.get(key);
        assert.equal(typeof want, 'number', `missing in-memory additive for ${key}`);
        assert.ok(Math.abs(csvAdditive - want) < 1e-9, `CSV vs memory mismatch for ${key}: ${csvAdditive} != ${want}`);
    }
    assert.equal(seen.size, 4);
});

test('rankCombinations still orders synergy most-negative-first', () => {
    const combos = [
        { combination: ['A'], levels: [0], additiveMean: 0.4, perBackground: [], backgrounds: 1 },
        { combination: ['B'], levels: [0], additiveMean: -0.9, perBackground: [], backgrounds: 1 },
    ];
    assert.equal(rankCombinations(combos, 'synergy')[0].combination[0], 'B');
});
