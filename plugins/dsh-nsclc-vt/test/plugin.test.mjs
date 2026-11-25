/**
 * Plugin wiring tests: load the plugin through a stub Cordis context, check it
 * registers exactly the documented tools with valid JSON Schemas, and execute
 * every tool once so a schema/argument mismatch cannot ship.
 *
 * These run without the harness: `lib/host.js` falls back to a local
 * `defineTool` that produces the same definition shape.
 *
 * @module dsh-nsclc-vt/test/plugin.test
 */
import { test } from 'node:test';
import assert from 'node:assert/strict';

import { apply, name as pluginName, inject, hostResolution } from '../lib/index.js';

/** Build a stub context that records registrations and prompt sections. */
function makeContext() {
    const tools = new Map();
    const sections = [];
    const ctx = {
        tools: {
            register(definition) { tools.set(definition.name, definition); },
        },
        inject(deps, fn) {
            if (deps.includes('systemPrompt')) {
                fn({ systemPrompt: { section(s) { sections.push(s); } } });
            }
        },
    };
    return { ctx, tools, sections };
}

const { ctx, tools, sections } = makeContext();
apply(ctx, {});

test('the plugin registers under the documented name and requires only `tools`', () => {
    assert.equal(pluginName, 'nsclc-vt');
    assert.deepEqual(inject, ['tools']);
});

test('exactly the four documented tools are registered', () => {
    assert.deepEqual([...tools.keys()].sort(), ['nsclc_model', 'nsclc_node', 'nsclc_predict', 'nsclc_screen']);
});

test('every tool declares a JSON-Schema parameter object and an output schema', () => {
    for (const [toolName, definition] of tools) {
        assert.equal(typeof definition.description, 'string');
        assert.ok(definition.description.length > 40, `${toolName} needs a substantive description`);
        assert.equal(definition.parameters.type, 'object');
        assert.equal(typeof definition.parameters.properties, 'object');
        for (const [param, schema] of Object.entries(definition.parameters.properties)) {
            assert.equal(typeof schema.type, 'string', `${toolName}.${param} must declare a type`);
            assert.ok(!Object.hasOwn(schema, 'required') || schema.required === true,
                `${toolName}.${param}.required must be the literal true when present`);
        }
        assert.equal(definition.output.schema.type, 'object');
        assert.equal(typeof definition.output.render, 'function');
        assert.equal(typeof definition.execute, 'function');
        assert.ok(definition.timeoutMs > 0);
    }
});

test('the system-prompt section is registered and explains the sign conventions', () => {
    assert.equal(sections.length, 1);
    assert.equal(sections[0].name, 'tool:nsclc-vt');
    assert.match(sections[0].text, /additive/);
    assert.match(sections[0].text, /negative means the pair beats its better single component/);
});

test('nsclc_model lists the model, the backgrounds and a node query', async () => {
    const base = await tools.get('nsclc_model').execute({}, {});
    assert.equal(base.nodeCount, 162);
    assert.deepEqual(base.phenotypes, ['Proliferation', 'Death']);
    assert.equal(base.backgrounds.length, 7);

    const hit = await tools.get('nsclc_model').execute({ query: 'tp53' }, {});
    assert.ok(hit.nodes.some((n) => n.node === 'p53'), 'gene symbol TP53 should resolve to the p53 node');
});

test('nsclc_node explains a target function and predicts under a genotype', async () => {
    const plain = await tools.get('nsclc_node').execute({ node: 'ATM' }, {});
    assert.equal(plain.node, 'ATM');
    assert.ok(plain.activators.includes('MRN'), 'ATM should be listed as activated by the MRN complex');
    assert.ok(plain.inhibitors.includes('ATM_deficient'), 'ATM should be inhibited by ATM_deficient');
    assert.match(plain.targetFunction, /MRN/, 'target function should name its inputs');
    assert.ok(!plain.targetFunction.includes('var('), 'node ids must be substituted by names');

    const withPrediction = await tools.get('nsclc_node').execute({ node: 'ATM', background: 'a549' }, {});
    assert.equal(withPrediction.prediction.background.toLowerCase(), 'a549');
    assert.ok(Number.isInteger(withPrediction.prediction.lo));
});

test('nsclc_predict reports phenotypes, difference and the moved nodes', async () => {
    const out = await tools.get('nsclc_predict').execute({ background: 'pc9', interventions: ['ATM:0'] }, {});
    assert.equal(out.background.toLowerCase(), 'pc9', 'published background casing is preserved');
    assert.deepEqual(out.phenotypes.map((p) => p.node), ['Proliferation', 'Death']);
    assert.equal(typeof out.difference, 'number');
    assert.ok(Array.isArray(out.changes));
    assert.ok(out.iterations > 0);
});

test('nsclc_predict rejects an unknown node with a clear error', async () => {
    await assert.rejects(
        () => tools.get('nsclc_predict').execute({ interventions: ['NOT_A_GENE:0'] }, {}),
        /unknown node or gene/,
    );
});

test('nsclc_screen ranks a small pair screen and honours the additivity filter', async () => {
    const out = await tools.get('nsclc_screen').execute({
        backgrounds: ['a549', 'healthy'],
        mode: 'both',
        candidates: ['MEK', 'PI3K', 'Akt'],
        levels: [0],
        top: 5,
    }, {});
    assert.equal(out.mode, 'both');
    assert.equal(out.analysesRun, 2 * (1 + 3 + 3));
    assert.ok(out.singles.length > 0);
    assert.ok(out.combinations.length > 0);
    for (const c of out.combinations) {
        assert.equal(c.perBackground.length, 2);
        for (const row of c.perBackground) {
            assert.equal(typeof row.additive, 'number');
        }
    }
    // Ranked by synergy: the first entry must be the most negative mean.
    const means = out.combinations.map((c) => c.synergyMean);
    assert.deepEqual(means, [...means].sort((a, b) => a - b));
});

test('nsclc_screen refuses an enumeration above its cap instead of running forever', async () => {
    await assert.rejects(
        () => tools.get('nsclc_screen').execute({ mode: 'pair', max_analyses: 10 }, {}),
        /above the 10 cap/,
    );
});

test('the tools report which host factory they resolved', () => {
    assert.ok(['harness', 'local-fallback'].includes(hostResolution.tools));
    assert.ok(['harness', 'absent'].includes(hostResolution.schema));
});

/**
 * Recursively check a value against a declared output schema, including the
 * `additionalProperties: false` contract. The harness validates tool output, so
 * an undeclared key would be a runtime failure in production — this catches it
 * at test time instead.
 */
function schemaViolations(schema, value, path = 'result') {
    const out = [];
    if (schema === undefined || schema === null) return out;
    if (typeof schema.type === 'string') {
        const ok = schema.type === 'array' ? Array.isArray(value)
            : schema.type === 'object' ? value !== null && typeof value === 'object' && !Array.isArray(value)
                : schema.type === 'integer' ? Number.isInteger(value)
                    : schema.type === 'number' ? typeof value === 'number'
                        : schema.type === 'string' ? typeof value === 'string'
                            : schema.type === 'boolean' ? typeof value === 'boolean'
                                : true;
        if (!ok) { out.push(`${path} should be ${schema.type}, got ${Array.isArray(value) ? 'array' : typeof value}`); return out; }
    }
    if (Array.isArray(schema.enum) && !schema.enum.includes(value)) out.push(`${path} is not one of ${schema.enum.join(', ')}`);
    if (schema.type === 'object' && value !== null && typeof value === 'object') {
        for (const name of schema.required ?? []) {
            if (value[name] === undefined) out.push(`${path}.${name} is required but missing`);
        }
        if (schema.additionalProperties === false) {
            for (const key of Object.keys(value)) {
                if (!Object.hasOwn(schema.properties ?? {}, key)) out.push(`${path}.${key} is not declared in the schema`);
            }
        }
        for (const [name, child] of Object.entries(schema.properties ?? {})) {
            if (value[name] !== undefined) out.push(...schemaViolations(child, value[name], `${path}.${name}`));
        }
    }
    if (schema.type === 'array' && Array.isArray(value) && schema.items !== undefined) {
        value.forEach((entry, i) => out.push(...schemaViolations(schema.items, entry, `${path}[${i}]`)));
    }
    return out;
}

test('every tool result satisfies its own declared output schema', async () => {
    const problems = [];
    const invocations = [
        ['nsclc_model', {}],
        ['nsclc_model', { query: 'atm' }],
        ['nsclc_node', { node: 'p53', background: 'a549', interventions: ['ATM:0'] }],
        ['nsclc_predict', { background: 'sw1573', top_changes: 5 }],
        ['nsclc_predict', { background: 'healthY', interventions: ['MYC:3'], compare_to: 'none', nodes: ['p53', 'ATM'], include_state: true }],
        ['nsclc_screen', { backgrounds: ['healthy'], mode: 'both', candidates: ['MEK', 'PI3K'], top: 3 }],
    ];
    for (const [toolName, args] of invocations) {
        const definition = tools.get(toolName);
        const value = await definition.execute(args, {});
        for (const violation of schemaViolations(definition.output.schema, value, `${toolName}(${JSON.stringify(args)})`)) {
            problems.push(violation);
        }
    }
    assert.deepEqual(problems, []);
});
