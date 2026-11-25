/**
 * `nsclc_model` — describe the virtual-tumor model and search its nodes.
 * @module dsh-nsclc-vt/tools/nsclc-model
 */
import { defineTool } from '../host.js';
import { summarizeModel } from '../model.js';
import { table, MODEL_PREDICTION_NOTE } from '../format.js';

const OUTPUT_SCHEMA = {
    type: 'object',
    additionalProperties: false,
    properties: {
        evidence_kind: { type: 'string', required: true, enum: ['model_prediction'], description: 'Every number here is a simulation output, not a measurement.' },
        evidence_note: { type: 'string', required: true, description: MODEL_PREDICTION_NOTE },
        model: { type: 'string', required: true, description: 'One-line model summary.' },
        nodeCount: { type: 'integer', required: true },
        relationshipCount: { type: 'integer', required: true },
        defaultFunctionCount: { type: 'integer', required: true, description: 'Nodes whose target function is synthesised from their edges rather than written out.' },
        phenotypes: { type: 'array', required: true, items: { type: 'string' } },
        backgrounds: { type: 'array', required: true, items: { type: 'string' }, description: 'Cell-line mutational backgrounds the manuscript screened.' },
        nodes: {
            type: 'array',
            description: 'Matching nodes, when a query was supplied.',
            items: {
                type: 'object',
                additionalProperties: false,
                properties: {
                    node: { type: 'string', required: true },
                    id: { type: 'integer', required: true },
                    gene: { type: 'string' },
                    protein: { type: 'string' },
                    range: { type: 'string', required: true },
                    activators: { type: 'integer', required: true },
                    inhibitors: { type: 'integer', required: true },
                    downstream: { type: 'integer', required: true },
                },
            },
        },
        totalMatches: { type: 'integer' },
        source: { type: 'string', required: true, description: 'Where the model came from.' },
    },
};

/**
 * @param {import('../analysis.js').createAnalyst} analyst
 * @param {{limit?: number, source?: string}} ctx
 */
export function makeNsclcModelTool(analyst, ctx = {}) {
    const model = analyst.model;
    const relationships = model.order.reduce((n, id) => n + model.activators.get(id).length + model.inhibitors.get(id).length, 0);
    const defaultCount = model.order.filter((id) => analyst.isDefaultFunction(id)).length;

    return defineTool({
        name: 'nsclc_model',
        description: 'Describe the NSCLC "virtual tumor" qualitative network model bundled with this plugin: size, phenotype readouts, the seven published cell-line mutational backgrounds, and — given a query — matching nodes with their wiring. Use it first to discover exact node names before calling nsclc_predict or nsclc_screen.',
        parameters: {
            query: { type: 'string', description: 'Optional substring (case-insensitive) matched against node names, gene symbols and protein names, e.g. "ATM", "TP53", "rad51".' },
            limit: { type: 'integer', description: 'Maximum matching nodes to return (default 40).' },
            include_backgrounds: { type: 'boolean', description: 'List the published mutational backgrounds (default true).' },
        },
        output: {
            schema: OUTPUT_SCHEMA,
            render: (_args, value) => {
                const lines = [value.model, `Phenotype readouts: ${value.phenotypes.join(', ')}`, `Source: ${value.source}`];
                if (Array.isArray(value.nodes) && value.nodes.length > 0) {
                    lines.push('', `Matching nodes (${value.totalMatches} total):`);
                    lines.push(table(
                        ['node', 'gene', 'range', 'act', 'inh', 'out'],
                        value.nodes.map((n) => [n.node, n.gene ?? '', n.range, n.activators, n.inhibitors, n.downstream]),
                        { align: ['l', 'l', 'l', 'r', 'r', 'r'] },
                    ));
                }
                if (Array.isArray(value.backgrounds) && value.backgrounds.length > 0) {
                    lines.push('', `Published backgrounds: ${value.backgrounds.join(', ')}`);
                }
                return [{ type: 'text', text: lines.join('\n') }];
            },
        },
        timeoutMs: 30000,
        isConcurrencySafe: () => true,
        execute(args) {
            const limit = Number.isInteger(args.limit) ? Math.max(1, Math.min(400, args.limit)) : 40;
            const query = typeof args.query === 'string' ? args.query.trim().toLowerCase() : '';
            let nodes;
            let totalMatches;
            if (query !== '') {
                const matches = model.order.filter((id) => {
                    const n = model.names.get(id).toLowerCase();
                    const gene = (analyst.geneOf(id) ?? '').toLowerCase();
                    const protein = (analyst.proteinOf(id) ?? '').toLowerCase();
                    return n.includes(query) || gene.includes(query) || protein.includes(query);
                });
                totalMatches = matches.length;
                nodes = matches.slice(0, limit).map((id) => {
                    const [lo, hi] = model.range.get(id);
                    return {
                        node: model.names.get(id),
                        id,
                        gene: analyst.geneOf(id),
                        protein: analyst.proteinOf(id),
                        range: `${lo}..${hi}`,
                        activators: model.activators.get(id).length,
                        inhibitors: model.inhibitors.get(id).length,
                        downstream: model.outputs.get(id).size,
                    };
                });
            }
            return {
                model: summarizeModel(model),
                evidence_kind: 'model_prediction',
                evidence_note: MODEL_PREDICTION_NOTE,
                nodeCount: model.order.length,
                relationshipCount: relationships,
                defaultFunctionCount: defaultCount,
                phenotypes: analyst.phenotypes,
                backgrounds: args.include_backgrounds === false
                    ? []
                    : analyst.tables.publishedBackgrounds.map((b) => b.name),
                ...(nodes ? { nodes, totalMatches } : {}),
                source: ctx.source ?? analyst.modelPath,
            };
        },
    });
}
