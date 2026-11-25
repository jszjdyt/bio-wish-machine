/**
 * `nsclc_node` — the mechanism card for one node of the virtual tumor:
 * its target function, its curated biology with citations, its wiring, and
 * optionally what it settles to in a given genotype.
 * @module dsh-nsclc-vt/tools/nsclc-node
 */
import { defineTool } from '../host.js';
import { level, signed, table, MODEL_PREDICTION_NOTE } from '../format.js';

const OUTPUT_SCHEMA = {
    type: 'object',
    additionalProperties: false,
    properties: {
        evidence_kind: { type: 'string', required: true, enum: ['model_prediction'], description: 'Every number here is a simulation output, not a measurement.' },
        evidence_note: { type: 'string', required: true, description: MODEL_PREDICTION_NOTE },
        node: { type: 'string', required: true },
        id: { type: 'integer', required: true },
        gene: { type: 'string' },
        protein: { type: 'string' },
        range: { type: 'string', required: true },
        targetFunction: { type: 'string', required: true, description: 'The node\'s effective target function with node names substituted, as the engine evaluates it.' },
        formulaIsSynthesised: { type: 'boolean', required: true, description: 'True when the model file leaves the formula empty and the importer synthesised it from the node\'s edges.' },
        curatedFunction: { type: 'string', description: 'The target function as curated in the model\'s companion table (may use pos/neg shorthand).' },
        curatedReasoning: { type: 'string' },
        citations: { type: 'array', items: { type: 'string' }, description: 'Sources cited for this node\'s target function.' },
        activators: { type: 'array', required: true, items: { type: 'string' } },
        inhibitors: { type: 'array', required: true, items: { type: 'string' } },
        activates: { type: 'array', required: true, items: { type: 'string' }, description: 'Nodes this one activates.' },
        inhibits: { type: 'array', required: true, items: { type: 'string' }, description: 'Nodes this one inhibits.' },
        edgeEvidence: {
            type: 'array',
            description: 'Curated edge evidence for this node\'s connections.',
            items: {
                type: 'object',
                additionalProperties: false,
                properties: {
                    from: { type: 'string', required: true },
                    to: { type: 'string', required: true },
                    type: { type: 'string', required: true },
                    citation: { type: 'string' },
                    doi: { type: 'string' },
                },
            },
        },
        prediction: {
            type: 'object',
            additionalProperties: false,
            description: 'Present only when a background or intervention was supplied.',
            properties: {
                background: { type: 'string', required: true },
                interventions: { type: 'string', required: true },
                steadyState: { type: 'string', required: true, description: 'The value, or interval, the node settles to.' },
                lo: { type: 'integer', required: true },
                hi: { type: 'integer', required: true },
                stabilizing: { type: 'boolean', required: true },
                interventionsApplied: { type: 'array', items: { type: 'string' } },
            },
        },
        note: { type: 'string' },
    },
};

/**
 * @param {ReturnType<import('../analysis.js').createAnalyst>} analyst
 */
export function makeNsclcNodeTool(analyst) {
    const model = analyst.model;

    return defineTool({
        name: 'nsclc_node',
        description: 'Explain one node of the NSCLC virtual-tumor model: its target function (the rule the simulation evaluates), its curated activators/inhibitors with the papers they came from, what it acts on downstream, and — if a genotype is given — the level or interval it settles to. Use it to justify why a predicted combination works, or to check how a gene is wired before proposing it as a target.',
        parameters: {
            node: { type: 'string', required: true, description: 'Node name, gene symbol or numeric variable id, e.g. "ATM", "TP53", "53BP1".' },
            background: { type: 'string', description: 'Optional genotype: a published background name ("a549", "healthy") or explicit levels ("RAS:3,STK11:0").' },
            interventions: { type: 'array', items: { type: 'string' }, description: 'Optional perturbations to apply on top of the background, e.g. ["ATM:0"]. A bare name means knockout at the node\'s minimum.' },
            neighbours: { type: 'integer', description: 'How many downstream/edge entries to list (default 12).' },
        },
        output: {
            schema: OUTPUT_SCHEMA,
            render: (_args, value) => {
                const lines = [`${value.node} (id ${value.id})${value.gene ? ` · ${value.gene}` : ''}${value.protein ? ` · ${value.protein}` : ''} — range ${value.range}`];
                lines.push(`Target function: ${value.targetFunction}${value.formulaIsSynthesised ? '  [synthesised from edges]' : ''}`);
                if (value.curatedFunction !== undefined && value.curatedFunction !== value.targetFunction) {
                    lines.push(`Curated form:    ${value.curatedFunction}`);
                }
                if (value.curatedReasoning) lines.push(`Rationale:       ${value.curatedReasoning}`);
                if (Array.isArray(value.citations) && value.citations.length > 0) lines.push(`Citations:       ${value.citations.join('; ')}`);
                lines.push(`Activated by:    ${value.activators.join(', ') || 'nothing'}`);
                lines.push(`Inhibited by:    ${value.inhibitors.join(', ') || 'nothing'}`);
                if (value.activates.length > 0) lines.push(`Activates:       ${value.activates.join(', ')}`);
                if (value.inhibits.length > 0) lines.push(`Inhibits:        ${value.inhibits.join(', ')}`);
                if (Array.isArray(value.edgeEvidence) && value.edgeEvidence.length > 0) {
                    lines.push('', 'Edge evidence:');
                    lines.push(table(['from', 'to', 'type', 'citation'], value.edgeEvidence.map((e) => [e.from, e.to, e.type, e.citation ?? ''])));
                }
                if (value.prediction !== undefined) {
                    lines.push('', `In ${value.prediction.background} with [${value.prediction.interventions}]: `
                        + `${value.node} settles to ${value.prediction.steadyState}`
                        + (value.prediction.stabilizing ? ' (model stabilized to a unique steady state)' : ' (model did not stabilize to a point for every node)'));
                }
                if (value.note) lines.push('', value.note);
                return [{ type: 'text', text: lines.join('\n') }];
            },
        },
        timeoutMs: 120000,
        isConcurrencySafe: () => true,
        execute(args) {
            const { id } = analyst.resolveNode(args.node);
            const name = model.names.get(id);
            const [lo, hi] = model.range.get(id);
            const neighbours = Number.isInteger(args.neighbours) ? Math.max(0, Math.min(200, args.neighbours)) : 12;

            const activates = [];
            const inhibits = [];
            for (const target of model.outputs.get(id)) {
                if (model.activators.get(target).includes(id)) activates.push(model.names.get(target));
                if (model.inhibitors.get(target).includes(id)) inhibits.push(model.names.get(target));
            }

            const edgeEvidence = [];
            for (const r of analyst.tables.edges) {
                if (!r.inModel) continue;
                if (r.from !== name && r.to !== name) continue;
                edgeEvidence.push({ from: r.from, to: r.to, type: r.type, citation: r.citation, doi: r.doi });
                if (edgeEvidence.length >= neighbours) break;
            }

            const curated = analyst.curatedFunctionOf(id);
            const citationList = (curated?.citation ?? '').split(/;\s*/).filter((s) => s !== '');

            /** @type {Record<string, unknown>} */
            const out = {
                node: name,
                evidence_kind: 'model_prediction',
                evidence_note: MODEL_PREDICTION_NOTE,
                id,
                gene: analyst.geneOf(id),
                protein: analyst.proteinOf(id),
                range: `${lo}..${hi}`,
                targetFunction: analyst.formulaOf(id, { effective: true }),
                formulaIsSynthesised: analyst.isDefaultFunction(id),
                ...(curated?.target ? { curatedFunction: curated.target } : {}),
                ...(curated?.reasoning ? { curatedReasoning: curated.reasoning } : {}),
                ...(citationList.length > 0 ? { citations: citationList } : {}),
                activators: model.activators.get(id).map((v) => model.names.get(v)),
                inhibitors: model.inhibitors.get(id).map((v) => model.names.get(v)),
                activates,
                inhibits,
                ...(edgeEvidence.length > 0 ? { edgeEvidence } : {}),
            };

            if (args.background !== undefined || (Array.isArray(args.interventions) && args.interventions.length > 0)) {
                const { background, interventions, analysis } = analyst.solve({
                    background: args.background,
                    interventions: args.interventions ?? [],
                });
                const [vl, vh] = analysis.bounds.get(id);
                out.prediction = {
                    background: background.name ?? (background.mutations.map((m) => `${m.name}=${m.level}`).join(', ') || 'unperturbed model'),
                    interventions: interventions.map((p) => `${p.name}=${p.level}`).join(', ') || 'none',
                    steadyState: level(vl, vh),
                    lo: vl,
                    hi: vh,
                    stabilizing: analysis.stabilizing,
                    interventionsApplied: interventions.map((p) => `${p.name}${p.level === lo ? '' : `=${signed(p.level, 0)}`}`),
                };
            }
            return out;
        },
    });
}
