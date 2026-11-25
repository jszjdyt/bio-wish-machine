/**
 * Host-service resolution.
 *
 * A DSH plugin normally writes `import { defineTool } from '@deepseek-ai/dsh-tools'`
 * and relies on Node finding the harness's own packages. That works when the
 * plugin sits inside the profile tree (`<DSH_HOME>/profiles/<name>/node_modules`),
 * because resolution walks up to `<DSH_HOME>/profiles/node_modules/@deepseek-ai`.
 * It does **not** work when the plugin is developed from a workspace and linked
 * in, because Node resolves symlinks to their real path first.
 *
 * This module therefore resolves the harness services by trying, in order: the
 * bare specifier, then the well-known profile locations. If none is reachable —
 * for example when running the plugin's own tests with plain `node` — it falls
 * back to a local `defineTool` that produces the identical definition shape, so
 * the engine and the tools stay runnable and testable everywhere.
 *
 * @module dsh-nsclc-vt/host
 */
import { existsSync, readdirSync } from 'node:fs';
import { homedir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const dshHome = () => process.env.DSH_HOME ?? join(homedir(), '.dsh');

/** Candidate absolute locations of a harness package, most specific first. */
function candidateUrls(pkg) {
    const home = dshHome();
    const roots = [
        join(home, 'profiles', 'node_modules'),
        join(home, 'node_modules'),
    ];
    // Also look inside every profile's own node_modules.
    try {
        const profiles = join(home, 'profiles');
        if (existsSync(profiles)) {
            for (const name of readdirSync(profiles)) {
                roots.push(join(profiles, name, 'node_modules'));
            }
        }
    } catch { /* unreadable profile directory is not fatal */ }
    const out = [];
    for (const root of roots) {
        const dir = join(root, ...pkg.split('/'));
        for (const entry of ['lib/index.js', 'index.js', 'lib/index.mjs', 'index.mjs']) {
            const file = join(dir, entry);
            if (existsSync(file)) out.push(pathToFileURL(file).href);
        }
    }
    return out;
}

/**
 * Import a harness package, trying the bare specifier then known locations.
 * @param {string} pkg
 * @returns {Promise<any|undefined>}
 */
async function importHost(pkg) {
    try {
        return await import(pkg);
    } catch { /* fall through to the filesystem candidates */ }
    for (const url of candidateUrls(pkg)) {
        try {
            return await import(url);
        } catch { /* try the next candidate */ }
    }
    return undefined;
}

/** Minimal JSON-Schema output for one value spec node. */
function compileValue(spec) {
    if (spec === null || typeof spec !== 'object') return spec;
    /** @type {Record<string, unknown>} */
    const out = {};
    for (const [key, value] of Object.entries(spec)) {
        if (key === 'required') continue;
        if (key === 'properties' && value !== null && typeof value === 'object') {
            const properties = {};
            const required = [];
            for (const [name, child] of Object.entries(value)) {
                properties[name] = compileValue(child);
                if (child !== null && typeof child === 'object' && child.required === true) required.push(name);
            }
            out.properties = properties;
            if (required.length > 0) out.required = required;
            continue;
        }
        if (key === 'items') { out.items = compileValue(value); continue; }
        out[key] = value;
    }
    return out;
}

/** Compile a parameters spec into an object-rooted JSON Schema. */
function compileParameters(spec) {
    const schema = compileValue({ type: 'object', properties: spec ?? {} });
    return { type: 'object', ...schema };
}

const typeMatches = (type, value) => {
    switch (type) {
        case 'string': return typeof value === 'string';
        case 'number': return typeof value === 'number' && Number.isFinite(value);
        case 'integer': return Number.isInteger(value);
        case 'boolean': return typeof value === 'boolean';
        case 'array': return Array.isArray(value);
        case 'object': return value !== null && typeof value === 'object' && !Array.isArray(value);
        case 'null': return value === null;
        default: return true;
    }
};

/** Report the first violation of `schema` by `value`, or undefined. */
function firstViolation(schema, value, path) {
    if (schema === null || typeof schema !== 'object') return undefined;
    if (typeof schema.type === 'string' && !typeMatches(schema.type, value)) {
        return `${path || 'arguments'} must be a ${schema.type}`;
    }
    if (Array.isArray(schema.enum) && !schema.enum.includes(value)) {
        return `${path || 'arguments'} must be one of ${schema.enum.join(', ')}`;
    }
    if (schema.type === 'object' && value !== null && typeof value === 'object' && !Array.isArray(value)) {
        for (const name of schema.required ?? []) {
            if (value[name] === undefined) return `missing required property "${path ? `${path}.` : ''}${name}"`;
        }
        for (const [name, child] of Object.entries(schema.properties ?? {})) {
            if (value[name] !== undefined) {
                const bad = firstViolation(child, value[name], path ? `${path}.${name}` : name);
                if (bad !== undefined) return bad;
            }
        }
    }
    if (schema.type === 'array' && Array.isArray(value) && schema.items !== undefined) {
        for (let i = 0; i < value.length; i += 1) {
            const bad = firstViolation(schema.items, value[i], `${path ? `${path}.` : ''}${i}`);
            if (bad !== undefined) return bad;
        }
    }
    return undefined;
}

/** The fallback definition factory; shape-identical to the harness `defineTool`. */
function localDefineTool(options) {
    const parameters = compileParameters(options.parameters);
    const schema = compileValue(options.output.schema);
    const userExecute = options.execute;
    const tool = {
        name: options.name,
        description: options.description,
        parameters,
        output: { schema, render: options.output.render },
        ...(options.timeoutMs !== undefined ? { timeoutMs: options.timeoutMs } : {}),
        ...(options.isConcurrencySafe !== undefined ? { isConcurrencySafe: options.isConcurrencySafe } : {}),
        ...(options.presentCall !== undefined ? { presentCall: options.presentCall } : {}),
        async execute(args, exec) {
            const bad = firstViolation(parameters, args ?? {}, '');
            if (bad !== undefined) throw new Error(`invalid arguments: ${bad}`);
            return userExecute(args ?? {}, exec);
        },
    };
    return tool;
}

const toolsModule = await importHost('@deepseek-ai/dsh-tools');
const schemaModule = await importHost('@deepseek-ai/schemastery');

/**
 * The harness `defineTool` when reachable, otherwise the local equivalent.
 * @type {typeof localDefineTool}
 */
export const defineTool = typeof toolsModule?.defineTool === 'function'
    ? toolsModule.defineTool
    : localDefineTool;

/** True when the real harness tool factory was resolved. */
export const usingHarnessTools = typeof toolsModule?.defineTool === 'function';

/**
 * Schemastery, when reachable. `undefined` means the plugin falls back to plain
 * defaults, which is only the case outside the harness.
 */
export const Schema = schemaModule?.default ?? schemaModule?.Schema;

/** Where the harness packages were found, for diagnostics. */
export const hostResolution = {
    tools: typeof toolsModule?.defineTool === 'function' ? 'harness' : 'local-fallback',
    schema: Schema === undefined ? 'absent' : 'harness',
};
