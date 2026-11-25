/**
 * dsh-nsclc-vt — the NSCLC "virtual tumor" as DeepSeek Harness tools.
 *
 * Wraps the interpretable mechanistic model of Clarke, Barker & Nicholls et al.,
 * *Virtual Tumors Enable Prediction of Personalized Therapeutic Combinations for
 * Non-Small Cell Lung Cancer* (Cancer Research 2026,
 * doi:10.1158/0008-5472.CAN-25-2476) and the BioModelAnalyzer VMCAI stability
 * engine it is analysed with, so the agent can run personalised combination
 * screens in-process without R, BioCheckConsole, Windows or a network call.
 *
 * The engine is a line-by-line port of `Stabilize.fs`, `GenLemmas.fs`,
 * `FNewLemmas.fs`, `Expr.fs` and `Marshal.fs` from
 * https://github.com/hallba/BioModelAnalyzer. It reproduces all 30,597 cells of
 * the manuscript's published combination-benefit matrix exactly; run
 * `node test/validate-published.mjs <additive_scores.csv>` to re-check that.
 *
 * @module dsh-nsclc-vt
 */
import { Schema, usingHarnessTools, hostResolution } from './host.js';
import { createAnalyst } from './analysis.js';
import { DEFAULT_MODEL_PATH } from './data.js';
import { defaultWorkers } from './screen.js';
import { makeNsclcModelTool } from './tools/nsclc-model.js';
import { makeNsclcNodeTool } from './tools/nsclc-node.js';
import { makeNsclcPredictTool } from './tools/nsclc-predict.js';
import { makeNsclcScreenTool } from './tools/nsclc-screen.js';
import { buildGuidance } from './prompt.js';

export { createAnalyst } from './analysis.js';
export { loadModel, withKnockouts, formulaOf, summarizeModel } from './model.js';
export { analyzeStability, readout, isConstantExpr } from './vmcai.js';
export { parseFormula, evalExpr, formatFormula } from './expr.js';
export { enumerateJobs, runJobs, runJobsInline, scoreScreen, rankCombinations } from './screen.js';
export { loadTables, DEFAULT_MODEL_PATH } from './data.js';
export { buildGuidance } from './prompt.js';
export { defineTool, Schema, usingHarnessTools, hostResolution } from './host.js';

/** Cordis plugin name used by loader diagnostics. */
export const name = 'nsclc-vt';

/** Services required before `apply` runs; `systemPrompt` is read lazily and stays optional. */
export const inject = ['tools'];

/** Configuration defaults, mirroring the schema below when schemastery is unavailable. */
const DEFAULTS = {
    modelPath: DEFAULT_MODEL_PATH,
    phenotypes: ['Proliferation', 'Death'],
    maxAnalyses: 60000,
    workers: 0,
    promptGuidance: true,
    promptOrder: 120,
};

/** Plugin configuration. */
export const Config = Schema === undefined ? undefined : Schema.object({
    modelPath: Schema.string()
        .default(DEFAULT_MODEL_PATH)
        .description('Path to a BMA model JSON. Defaults to the bundled NSCLC-master model.'),
    phenotypes: Schema.array(Schema.string())
        .default(['Proliferation', 'Death'])
        .description('Nodes reported as phenotype readouts. The manuscript screens on Proliferation and Death.'),
    maxAnalyses: Schema.number()
        .default(60000)
        .description('Refuse a screen whose enumeration exceeds this many model runs.'),
    workers: Schema.number()
        .default(0)
        .description('Worker threads for large screens; 0 picks min(8, cores - 1) automatically.'),
    promptGuidance: Schema.boolean()
        .default(true)
        .description('Register a system-prompt section explaining when to reach for these tools.'),
    promptOrder: Schema.number()
        .default(120)
        .description('Ordering weight of the system-prompt section.'),
});

/**
 * Register the virtual-tumor tools.
 * @param {import('@deepseek-ai/cordis').Context} ctx
 * @param {Partial<typeof DEFAULTS>} [config]
 */
export function apply(ctx, config) {
    const resolved = { ...DEFAULTS, ...(config ?? {}) };
    if (!Array.isArray(resolved.phenotypes) || resolved.phenotypes.length === 0) {
        throw new Error('nsclc-vt: `phenotypes` must list at least one model node');
    }
    if (typeof resolved.modelPath !== 'string' || resolved.modelPath.trim() === '') {
        throw new Error('nsclc-vt: `modelPath` must be a non-empty path');
    }
    if (Number.isFinite(resolved.workers) && resolved.workers > 0 && !Number.isInteger(resolved.workers)) {
        throw new Error('nsclc-vt: `workers` must be a whole number');
    }

    // Fail loudly at load time rather than on the first tool call: a bad model
    // path or an unknown phenotype is a configuration mistake, not a runtime one.
    const analyst = createAnalyst({
        modelPath: resolved.modelPath,
        phenotypes: resolved.phenotypes,
    });

    const workerCount = Number.isFinite(resolved.workers) && resolved.workers > 0
        ? resolved.workers
        : defaultWorkers();

    ctx.tools.register(makeNsclcModelTool(analyst, { source: resolved.modelPath }));
    ctx.tools.register(makeNsclcScreenTool(analyst, {
        maxAnalyses: Number.isFinite(resolved.maxAnalyses) ? resolved.maxAnalyses : 60000,
        workers: workerCount,
    }));
    ctx.tools.register(makeNsclcNodeTool(analyst));
    ctx.tools.register(makeNsclcPredictTool(analyst));

    if (resolved.promptGuidance) {
        ctx.inject(['systemPrompt'], (promptCtx) => {
            promptCtx.systemPrompt.section({
                name: 'tool:nsclc-vt',
                order: resolved.promptOrder,
                text: buildGuidance({
                    modelName: analyst.model.name,
                    nodeCount: analyst.model.order.length,
                    phenotypes: analyst.phenotypes,
                    backgrounds: analyst.tables.publishedBackgrounds.map((b) => b.name),
                }),
            });
        });
    }
}
