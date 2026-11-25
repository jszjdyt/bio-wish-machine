/**
 * `nsclc_predict` — run the virtual tumor for one genotype + intervention set
 * and report where the network settles.
 * @module dsh-nsclc-vt/tools/nsclc-predict
 */
import { defineTool } from '../host.js';
import { level, num, signed, table, DIFFERENCE_NOTE, MODEL_PREDICTION_NOTE } from '../format.js';
import { readout } from '../vmcai.js';

const READOUT_SCHEMA = {
    type: 'object',
    additionalProperties: false,
    properties: {
        node: { type: 'string', required: true },
        lo: { type: 'integer', required: true },
        hi: { type: 'integer', required: true },
        mean: { type: 'number', required: true, description: 'Interval midpoint — the value NANSEN screens on.' },
        uncertainty: { type: 'number', required: true, description: '(hi - lo) / (range span); 0 means the node is pinned to a single level.' },
    },
};

const OUTPUT_SCHEMA = {
    type: 'object',
    additionalProperties: false,
    properties: {
        evidence_kind: { type: 'string', required: true, enum: ['model_prediction'], description: 'Every number here is a simulation output, not a measurement.' },
        evidence_note: { type: 'string', required: true, description: MODEL_PREDICTION_NOTE },
        background: { type: 'string', required: true },
        interventions: { type: 'string', required: true },
        stabilizing: { type: 'boolean', required: true, description: 'True when every one of the model\'s nodes was pinned to a single level — the network proved a unique steady state.' },
        converged: { type: 'boolean', required: true, description: 'False only if the iteration cap was reached before the bounds stopped tightening.' },
        iterations: { type: 'integer', required: true, description: 'BMA\'s refinement-step count (the `Time` field in BioCheckConsole output).' },
        phenotypes: { type: 'array', required: true, items: READOUT_SCHEMA },
        difference: { type: 'number', description: DIFFERENCE_NOTE },
        changes: {
            type: 'array',
            description: 'Nodes that moved most relative to the comparison case, largest absolute change first.',
            items: {
                type: 'object',
                additionalProperties: false,
                properties: {
                    node: { type: 'string', required: true },
                    before: { type: 'number', required: true },
                    after: { type: 'number', required: true },
                    delta: { type: 'number', required: true },
                    afterRange: { type: 'string', required: true },
                },
            },
        },
        totalChanged: { type: 'integer', description: 'How many nodes moved, before the top_changes cut.' },
        state: { type: 'array', description: 'Full state, only when include_state is set.', items: READOUT_SCHEMA },
        degenerateNodes: { type: 'array', items: { type: 'string' }, description: 'Nodes left as an interval rather than a point, only present when the model did not fully stabilize.' },
        note: { type: 'string' },
    },
};

/**
 * @param {ReturnType<import('../analysis.js').createAnalyst>} analyst
 */
export function makeNsclcPredictTool(analyst) {
    const model = analyst.model;

    return defineTool({
        name: 'nsclc_predict',
        description: `Run the NSCLC virtual-tumor model forward for one genotype and a set of interventions (knockouts, activations or drug effects), and report the steady-state level of every phenotype plus the nodes that moved most. Use it to ask "in this mutational background, what does inhibiting X (and Y) do to proliferation and death?". ${DIFFERENCE_NOTE}`,
        parameters: {
            background: { type: 'string', description: 'Genotype: a published cell-line background ("a549", "pc9", "healthy") or explicit levels ("RAS:3,STK11:0,CDKN2A:0"). Omit for the unperturbed model.' },
            interventions: {
                type: 'array',
                items: { type: 'string' },
                description: 'Perturbations applied on top of the background, e.g. ["ATM:0", "MYC:3"]. A bare node name means a knockout at that node\'s minimum level. Interventions take precedence over the background on conflicting nodes, as in the manuscript pipeline.',
            },
            compare_to: { type: 'string', enum: ['background', 'none'], description: 'What to diff the changed-node list against: the background alone (default) or nothing.' },
            top_changes: { type: 'integer', description: 'How many moved nodes to return (default 15, 0 to skip).' },
            nodes: { type: 'array', items: { type: 'string' }, description: 'Optional whitelist of nodes to report individually, in addition to the phenotypes.' },
            include_state: { type: 'boolean', description: 'Return all model nodes, not just the phenotypes and the moved ones (default false).' },
        },
        output: {
            schema: OUTPUT_SCHEMA,
            render: (_args, value) => {
                const lines = [`Background: ${value.background}`, `Interventions: ${value.interventions}`];
                lines.push('');
                lines.push(table(['node', 'level', 'mean', 'uncertainty'], value.phenotypes.map((p) => [
                    p.node, level(p.lo, p.hi), num(p.mean), num(p.uncertainty),
                ]), { align: ['l', 'r', 'r', 'r'] }));
                if (value.difference !== undefined) lines.push(`Difference (Proliferation - Death): ${num(value.difference)}  — lower is better`);
                lines.push('');
                lines.push(value.stabilizing
                    ? `Model proved a unique steady state in ${value.iterations} refinement steps.`
                    : `Model did NOT stabilize to a point (${value.iterations} refinement steps).`
                        + (value.degenerateNodes && value.degenerateNodes.length > 0
                            ? ` Interval-valued nodes: ${value.degenerateNodes.join(', ')}.`
                            : ''));
                if (Array.isArray(value.changes) && value.changes.length > 0) {
                    lines.push('', `Largest changes vs ${value.background === 'unperturbed model' ? 'the unperturbed model' : 'the background alone'}:`);
                    lines.push(table(['node', 'before', 'after', 'delta', 'range'], value.changes.map((c) => [
                        c.node, num(c.before), num(c.after), signed(c.delta), c.afterRange,
                    ]), { align: ['l', 'r', 'r', 'r', 'l'] }));
                }
                if (value.note) lines.push('', value.note);
                return [{ type: 'text', text: lines.join('\n') }];
            },
        },
        timeoutMs: 180000,
        isConcurrencySafe: () => true,
        execute(args) {
            const compareTo = args.compare_to ?? 'background';
            const topChanges = Number.isInteger(args.top_changes) ? Math.max(0, Math.min(162, args.top_changes)) : 15;
            const interventions = Array.isArray(args.interventions) ? args.interventions : [];

            const primary = analyst.solve({ background: args.background, interventions });
            const pheno = analyst.phenotypeReadout(primary.analysis);

            const response = {
                evidence_kind: 'model_prediction',
                evidence_note: MODEL_PREDICTION_NOTE,
                background: primary.background.name
                    ?? (primary.background.mutations.map((m) => `${m.name}=${m.level}`).join(', ') || 'unperturbed model'),
                interventions: primary.interventions.map((p) => `${p.name}=${p.level}`).join(', ') || 'none',
                stabilizing: primary.analysis.stabilizing,
                converged: primary.analysis.converged,
                iterations: primary.analysis.iterations,
                phenotypes: pheno.phenotypes.map((p) => ({
                    node: p.node, lo: p.lo, hi: p.hi, mean: p.mean, uncertainty: p.uncertainty,
                })),
                ...(pheno.difference !== undefined ? { difference: pheno.difference } : {}),
            };

            const explicit = [];
            for (const token of args.nodes ?? []) {
                let id;
                try { ({ id } = analyst.resolveNode(token)); } catch { continue; }
                const r = readout(model, primary.analysis.bounds, id);
                explicit.push({ node: model.names.get(id), lo: r.lo, hi: r.hi, mean: r.mean, uncertainty: r.uncertainty });
            }

            if (compareTo === 'background') {
                const baseline = analyst.solve({ background: args.background });
                const changes = [];
                for (const id of model.order) {
                    const before = readout(model, baseline.analysis.bounds, id);
                    const after = readout(model, primary.analysis.bounds, id);
                    if (before.mean === after.mean && before.lo === after.lo && before.hi === after.hi) continue;
                    changes.push({
                        node: model.names.get(id),
                        before: before.mean,
                        after: after.mean,
                        delta: after.mean - before.mean,
                        afterRange: level(after.lo, after.hi),
                    });
                }
                changes.sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta) || a.node.localeCompare(b.node));
                response.changes = changes.slice(0, topChanges);
                response.totalChanged = changes.length;
            }

            if (args.include_state === true) {
                response.state = model.order.map((id) => {
                    const r = readout(model, primary.analysis.bounds, id);
                    return { node: model.names.get(id), lo: r.lo, hi: r.hi, mean: r.mean, uncertainty: r.uncertainty };
                });
            } else if (explicit.length > 0) {
                response.state = explicit;
            }

            const degenerate = model.order
                .filter((id) => primary.analysis.bounds.get(id)[0] !== primary.analysis.bounds.get(id)[1])
                .map((id) => model.names.get(id));
            if (!primary.analysis.stabilizing && degenerate.length > 0) {
                response.degenerateNodes = degenerate.slice(0, 40);
                const shown = degenerate.length > 40 ? ` (showing 40 of ${degenerate.length})` : '';
                response.note = `The interval analysis did not reduce every node to a single level, so the reported levels for `
                    + `${degenerate.join(', ')}${shown} are the fixed-point interval, not one trajectory. Phenotype readouts are still the manuscript's convention (interval midpoint).`;
            }
            return response;
        },
    });
}
