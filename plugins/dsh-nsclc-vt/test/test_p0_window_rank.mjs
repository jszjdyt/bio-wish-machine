/**
 * P0 #11 回归：window 排序方向 + 缺 Healthy 明确报错。
 *
 * 修复前 `lib/screen.js:rankCombinations()` 里
 *   benefit = -mean(benefit)   （benefit = -difference，越大越好）
 *   window  = -mean(difference)（= mean(benefit) 的相反数 → 排序方向与 benefit 相反，
 *                                且没有任何输入约束，缺 Healthy 组时仍在算）
 * window 没有定义。修复后 window 定义为
 *   「肿瘤背景平均 benefit − Healthy 背景 benefit」（越大 = 越肿瘤选择性，排前），
 * 缺 Healthy 背景时抛错说明为什么不能算。
 */
import { test } from 'node:test';
import assert from 'node:assert/strict';

import { rankCombinations } from '../lib/screen.js';

const combo = (id, rows) => ({
    combination: [id], levels: [0], perBackground: rows,
    additiveMean: 0, additiveSum: 0, improved: 0, worse: 0, equal: 0, backgrounds: rows.length,
});

const tumourSelective = combo('tumour-selective', [
    { background: 'A549', difference: -5, benefit: 5 },
    { background: 'Healthy', difference: 1, benefit: -1 },
]);
const systemic = combo('systemic', [
    { background: 'A549', difference: -1, benefit: 1 },
    { background: 'Healthy', difference: -2, benefit: 2 },
]);

test('window ranks the tumour-selective combination first', () => {
    // window(tumour-selective) = 5 - (-1) = 6; window(systemic) = 1 - 2 = -1
    const ranked = rankCombinations([systemic, tumourSelective], 'window');
    assert.equal(ranked[0].combination[0], 'tumour-selective');
    assert.equal(ranked[1].combination[0], 'systemic');
});

test('window refuses to guess when there is no Healthy background', () => {
    const noHealthy = combo('tumour-only', [
        { background: 'A549', difference: -5, benefit: 5 },
        { background: 'PC9', difference: -3, benefit: 3 },
    ]);
    assert.throws(() => rankCombinations([noHealthy], 'window'), /healthy/i);
});

test('window also refuses when only the Healthy background was screened', () => {
    const onlyHealthy = combo('healthy-only', [
        { background: 'Healthy', difference: 1, benefit: -1 },
    ]);
    assert.throws(() => rankCombinations([onlyHealthy], 'window'), /(tumour|tumor)/i);
});

test('benefit ranking is unchanged (best mean benefit first)', () => {
    const ranked = rankCombinations([systemic, tumourSelective], 'benefit');
    assert.equal(ranked[0].combination[0], 'tumour-selective');
});

test('death ranking is unchanged (most cell death first)', () => {
    const low = combo('low-death', [{ background: 'A549', difference: -1, benefit: 1, death: 0.1 }]);
    const high = combo('high-death', [{ background: 'A549', difference: -1, benefit: 1, death: 0.9 }]);
    assert.equal(rankCombinations([low, high], 'death')[0].combination[0], 'high-death');
});
