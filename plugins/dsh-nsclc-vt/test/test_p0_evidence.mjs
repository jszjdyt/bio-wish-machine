/**
 * P0 #12 回归：把 additive 改述为「较模型内最佳单干预的差值」，并给工具结果
 * 显式标 `evidence_kind: 'model_prediction'`（模型预测，不是实验结论）。
 */
import { test } from 'node:test';
import assert from 'node:assert/strict';

import { createAnalyst } from '../lib/analysis.js';
import { ADDITIVE_NOTE, MODEL_PREDICTION_NOTE } from '../lib/format.js';
import { makeNsclcModelTool } from '../lib/tools/nsclc-model.js';
import { makeNsclcNodeTool } from '../lib/tools/nsclc-node.js';
import { makeNsclcPredictTool } from '../lib/tools/nsclc-predict.js';
import { makeNsclcScreenTool } from '../lib/tools/nsclc-screen.js';

const analyst = createAnalyst();

test('ADDITIVE_NOTE is reframed as a predicted difference vs the best single intervention', () => {
    assert.match(ADDITIVE_NOTE, /best single/i);
    assert.match(ADDITIVE_NOTE, /model prediction/i);
    assert.doesNotMatch(ADDITIVE_NOTE, /^additive =/);
});

test('MODEL_PREDICTION_NOTE says it is not an experimental result', () => {
    assert.match(MODEL_PREDICTION_NOTE, /not (an )?experimental/i);
});

test('every tool result carries evidence_kind = model_prediction', async () => {
    const modelOut = await makeNsclcModelTool(analyst).execute({}, {});
    const nodeOut = await makeNsclcNodeTool(analyst).execute({ node: 'ATM' }, {});
    const predictOut = await makeNsclcPredictTool(analyst).execute({ background: 'pc9', interventions: ['ATM:0'] }, {});
    const screenOut = await makeNsclcScreenTool(analyst, { workers: 1 }).execute({
        backgrounds: ['healthy'], mode: 'mono', candidates: ['MEK'], levels: [0],
    }, {});
    for (const [name, out] of [
        ['nsclc_model', modelOut], ['nsclc_node', nodeOut],
        ['nsclc_predict', predictOut], ['nsclc_screen', screenOut],
    ]) {
        assert.equal(out.evidence_kind, 'model_prediction', `${name} must mark model_prediction`);
        assert.equal(typeof out.evidence_note, 'string', `${name} must carry evidence_note`);
        assert.match(out.evidence_note, /model prediction|not experimental|experiment/i, `${name} note`);
    }
});
