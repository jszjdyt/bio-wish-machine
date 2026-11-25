/**
 * P0-2 / P0-3 回归测试（SingleCell）：多步流程中切数据集不漂移、
 * 图按不可变 run_id 落盘、同参重画幂等、不同参数产生新 run_id、旧图可回查。
 *
 * 真实 lib/engine.mjs + 真实 lib/impl.mjs，Python 侧换成 JSON-lines 假 worker。
 *
 * 运行：node test/test_p0_state_freeze.mjs
 */

import { mkdtempSync, existsSync, readFileSync, readdirSync, rmSync, writeFileSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { createEngine } from '../lib/engine.mjs'
import * as impl from '../lib/impl.mjs'

const PLUGIN = 'singlecell'
const TMP = mkdtempSync(path.join(tmpdir(), `dsh-p0-state-${PLUGIN}-`))
const workerPath = path.join(TMP, 'fake_worker.py')

const FAKE_WORKER = `#!/usr/bin/env python3
import json, os, sys, time

def arg(name, default=None):
    try:
        return sys.argv[sys.argv.index(name) + 1]
    except Exception:
        return default

WORKDIR = arg('--workdir', '.')

def note(cmd, args):
    try:
        with open(os.path.join(WORKDIR, 'calls.jsonl'), 'a', encoding='utf-8') as handle:
            handle.write(json.dumps({'cmd': cmd, 'args': args}, ensure_ascii=False) + '\\n')
    except Exception:
        pass

def reply(msg_id, result=None, error=None):
    payload = {'id': msg_id, 'ok': error is None}
    if error is None:
        payload['result'] = result if result is not None else {}
    else:
        payload['error'] = {'code': 'internal', 'message': error}
    sys.stdout.write(json.dumps(payload) + '\\n')
    sys.stdout.flush()

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    try:
        msg = json.loads(line)
    except Exception:
        continue
    cmd = msg.get('cmd')
    args = msg.get('args') or {}
    note(cmd, args)
    if cmd == 'ping':
        reply(msg.get('id'), {'mode': 'builtin', 'packages': {}, 'layouts': ['umap'], 'platform': 'test'})
    elif cmd in ('list', 'state'):
        reply(msg.get('id'), {'datasets': [
            {'id': 'ds1', 'name': 'A', 'n_obs': 100, 'n_vars': 50, 'steps': []},
            {'id': 'ds2', 'name': 'B', 'n_obs': 200, 'n_vars': 60, 'steps': []},
        ], 'active': 'ds1'})
    elif cmd == 'recipe':
        rid = str(args.get('id') or '')
        steps = [{'step': 'qc'}, {'step': 'pca'}]
        reply(msg.get('id'), {'id': rid, 'title': rid, 'steps': steps} if rid else {'recipes': [{'id': 'standard', 'title': 'standard', 'steps': steps}]})
    elif cmd == 'run':
        reply(msg.get('id'), {'dataset': {'id': args.get('dataset')},
                              'summary': {'title': args.get('step'), 'bullets': [str(args.get('step')) + ' ok']},
                              'ms': 1})
    elif cmd == 'plot':
        kind = str(args.get('kind') or 'plot')
        name = 'f_%s_%s.png' % (kind, str(int(time.time() * 1000))[-6:])
        target = os.path.join(WORKDIR, 'figures', name)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with open(target, 'wb') as handle:
            handle.write(b'\\x89PNG\\r\\n\\x1a\\n fake')
        reply(msg.get('id'), {'id': name, 'file': target, 'filename': name, 'width': 100, 'height': 80,
                              'dataset': {'id': args.get('dataset')}})
    elif cmd == 'shutdown':
        reply(msg.get('id'), {'bye': True})
        break
    else:
        reply(msg.get('id'), {'echo': cmd})
`

writeFileSync(workerPath, FAKE_WORKER, { mode: 0o755 })

let passed = 0
let failed = 0
const failures = []

function check(name, condition, detail = '') {
  if (condition) {
    passed += 1
    console.log(`  PASS ${name}${detail ? `  [${detail}]` : ''}`)
  } else {
    failed += 1
    failures.push(`${name}${detail ? `：${detail}` : ''}`)
    console.log(`  FAIL ${name}${detail ? `  [${detail}]` : ''}`)
  }
}

function readCalls(workdir) {
  const file = path.join(workdir, 'calls.jsonl')
  if (!existsSync(file)) return []
  return readFileSync(file, 'utf8').split('\n').filter(Boolean).map((line) => {
    try { return JSON.parse(line) } catch { return { cmd: '?', args: {} } }
  })
}

function stagingLeftovers(workdir) {
  try { return readdirSync(path.join(workdir, 'staging')) } catch { return [] }
}

async function main() {
  const workdir = path.join(TMP, 'work')
  const engine = createEngine({
    workdir,
    engineScript: workerPath,
    python: process.env.DSH_TEST_PYTHON || 'python3',
  })
  const runtime = { engine, config: {}, ctx: null, activity: [], jobs: new Map() }
  // 单细胞工作台的运行时状态（resolveDatasetId 读这里）
  runtime.sc = {
    active: 'ds1',
    datasets: new Map([['ds1', { id: 'ds1' }], ['ds2', { id: 'ds2' }]]),
    activity: [],
    figures: new Map(),
    suggestCache: { key: '', value: [], at: 0 },
    lastRefresh: 0,
    bootAt: Date.now(),
  }
  check('0 引擎与 run store 就绪', (await engine.start()).mode === 'builtin' && Boolean(engine.store))

  // -------------------------------------------------------------------------
  // 1. 多步流程执行期间切活动数据集：后续步骤仍用开始时解析出的 ds1
  // -------------------------------------------------------------------------
  const recipePromise = impl.opRecipe(runtime, { id: 'standard' })
  await new Promise((resolve) => setTimeout(resolve, 80))
  runtime.sc.active = 'ds2' // 模拟用户在 UI 里切数据集
  const recipe = await recipePromise
  const runCalls = readCalls(workdir).filter((item) => item.cmd === 'run')
  check('1a 预设流程能跑完', recipe.completed === 2, `completed=${recipe.completed}`)
  check('1b 每一步都显式带开始时解析的 dataset（切数据集不漂移）',
    runCalls.length === 2 && runCalls.every((item) => item.args && item.args.dataset === 'ds1'),
    JSON.stringify(runCalls.map((item) => item.args && item.args.dataset)))

  // -------------------------------------------------------------------------
  // 2. 图：不可变 run 目录 + sha256 + 可回查
  // -------------------------------------------------------------------------
  const first = await impl.opPlot(runtime, { dataset: 'ds1', kind: 'embedding', params: { color_by: 'clusters' } })
  check('2a 出图返回 run_id', Boolean(first.run_id), `run=${first.run_id}`)
  const run = engine.getRun(first.run_id)
  const artifact = run && run.artifacts.find((item) => item.kind === 'figure')
  check('2b 图在 runs/<run_id>/figures 下且存在', Boolean(artifact) && artifact.file.startsWith(run.dir) &&
    existsSync(artifact.file), artifact ? artifact.file : 'no-artifact')
  check('2c 产物带 sha256', Boolean(artifact) && /^[0-9a-f]{64}$/.test(artifact.sha256 || ''),
    artifact ? String(artifact.sha256).slice(0, 12) : 'no-hash')
  check('2d 内存登记丢失后仍能按文件名回查', (() => {
    runtime.sc.figures.clear()
    const entry = impl.opFigure(runtime, artifact.name)
    return entry.file === artifact.file
  })(), artifact ? artifact.name : 'n/a')

  // -------------------------------------------------------------------------
  // 3. 同参重画幂等 / 不同参数新 run
  // -------------------------------------------------------------------------
  const plotCountBefore = readCalls(workdir).filter((item) => item.cmd === 'plot').length
  const second = await impl.opPlot(runtime, { dataset: 'ds1', kind: 'embedding', params: { color_by: 'clusters' } })
  const plotCountAfter = readCalls(workdir).filter((item) => item.cmd === 'plot').length
  check('3a 同参数重画复用同一 run_id', second.run_id === first.run_id && second.run_reused === true,
    `first=${first.run_id} second=${second.run_id}`)
  check('3b 幂等复用时不再渲染', plotCountAfter === plotCountBefore, `before=${plotCountBefore} after=${plotCountAfter}`)
  const third = await impl.opPlot(runtime, { dataset: 'ds1', kind: 'embedding', params: { color_by: 'cell_type' } })
  check('3c 不同参数产生新 run_id', third.run_id !== first.run_id, `a=${first.run_id} b=${third.run_id}`)

  // -------------------------------------------------------------------------
  // 4. 失败不新增已提交 run、staging 无残留
  // -------------------------------------------------------------------------
  const before = engine.listRuns({ limit: 100 }).filter((item) => item.status === 'committed').length
  let threw = false
  try {
    await impl.opPlot(runtime, { dataset: 'ds-never', kind: 'embedding' })
  } catch (error) {
    threw = Boolean(error && ['not_found', 'bad_request'].includes(error.code))
  }
  const after = engine.listRuns({ limit: 100 }).filter((item) => item.status === 'committed').length
  check('4a 不存在的数据集明确失败', threw === true)
  check('4b 失败请求不新增已提交 run', after === before, `before=${before} after=${after}`)
  check('4c staging 无残留', stagingLeftovers(workdir).length === 0, JSON.stringify(stagingLeftovers(workdir)))

  await engine.dispose()
  try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }

  console.log(`\n== P0-2/P0-3 状态回归（${PLUGIN}）：${passed} 通过 / ${failed} 失败 ==`)
  if (failed) {
    for (const item of failures) console.log(`  ✗ ${item}`)
    process.exit(1)
  }
}

main().catch((error) => {
  console.error('测试本身异常：', error)
  try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }
  process.exit(1)
})
