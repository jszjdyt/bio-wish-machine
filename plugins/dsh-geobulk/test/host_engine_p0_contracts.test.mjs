/**
 * dsh-* — P0-2 / P0-3 契约回归（任务单 §7 第 2、3 条）
 *
 *   ② 冻结每个任务的数据版本；逐步状态；下载按哈希复用；参考对象用 artifact ID；
 *   ③ 图/表/日志按不可变 run_id 落独立目录；缓存键含输入与参数版本；
 *
 * 用真实 engine.mjs + 假 Python worker 驱动，覆盖：
 *   - freezeDataset / resolveDatasetRef：活动数据集切换、数据 revision 变化 → stale_artifact；
 *   - withTask：同参数幂等复用 / 不同参数新 run / 同键修订 -v2 / 单步失败不得 ok:true；
 *   - staging → runs/<run_id> 原子提交，失败/超时 run 不产生任何已提交产物；
 *   - 下载文件按内容 sha256 复用、期望哈希不符报错、缓存自校验；
 *   - kill worker 后已提交 run 不被污染；进程重启后仍能按 run_id 回查并幂等复用。
 *
 * 运行：node test/host_engine_p0_contracts.test.mjs
 * 每个测试都有外置超时保护（watchdog + 单请求 withTimeout + 测试脚本最终 process.exit）。
 */

import { mkdtempSync, writeFileSync, readFileSync, existsSync, rmSync, mkdirSync, readdirSync, statSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { createEngine } from '../lib/engine.mjs'

const TMP = mkdtempSync(path.join(tmpdir(), 'dsh-engine-p0-'))
const WORK = path.join(TMP, 'work')
const WORKER_LOG = path.join(TMP, 'worker.log')

const FAKE_WORKER = `#!/usr/bin/env python3
import json, os, sys, time

LOG = ${JSON.stringify(WORKER_LOG)}

datasets = [
    {"id": "ds1", "name": "A", "n_obs": 10, "n_vars": 5, "steps": ["qc"]},
    {"id": "ds2", "name": "B", "n_obs": 8, "n_vars": 5, "steps": []},
]
active = "ds1"

def reply(msg_id, result=None, error=None):
    payload = {"id": msg_id, "ok": error is None}
    if error is None:
        payload["result"] = result if result is not None else {}
    else:
        payload["error"] = {"code": "internal", "message": error}
    sys.stdout.write(json.dumps(payload) + "\\n")
    sys.stdout.flush()

def record(cmd, args):
    with open(LOG, "a", encoding="utf-8") as handle:
        handle.write(cmd + "|" + json.dumps(args) + "\\n")
        handle.flush()

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    try:
        msg = json.loads(line)
    except Exception:
        continue
    cmd = msg.get("cmd")
    args = msg.get("args") or {}
    record(cmd, args)
    if cmd == "ping":
        reply(msg.get("id"), {"mode": "builtin", "packages": {}, "layouts": ["umap"], "platform": "test"})
    elif cmd == "list":
        reply(msg.get("id"), {"datasets": datasets, "active": active})
    elif cmd == "set_datasets":
        if "datasets" in args:
            datasets = args.get("datasets") or []
        if "active" in args:
            active = args.get("active")
        reply(msg.get("id"), {"datasets": datasets, "active": active})
    elif cmd == "echo":
        reply(msg.get("id"), {"echo": args.get("value")})
    elif cmd == "block":
        time.sleep(float(args.get("seconds") or 1.0))
        reply(msg.get("id"), {"blocked": True})
    elif cmd == "kill":
        reply(msg.get("id"), {"bye": True})
        sys.stdout.flush()
        os._exit(9)
    elif cmd == "selfdestruct":
        # 睡指定秒数后直接退出、不回复：模拟任务执行中途 worker 被杀
        time.sleep(float(args.get("seconds") or 1.0))
        os._exit(9)
    elif cmd == "shutdown":
        reply(msg.get("id"), {"bye": True})
        break
    else:
        reply(msg.get("id"), error="unknown cmd " + str(cmd))
`

const workerPath = path.join(TMP, 'fake_worker.py')
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

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms))

function withTimeout(promise, ms, label) {
  return Promise.race([
    promise.then((value) => ({ settled: 'ok', value })).catch((error) => ({ settled: 'error', error })),
    new Promise((resolve) => setTimeout(() => resolve({ settled: 'hung', label }), ms)),
  ])
}

const DATASETS_R0 = [
  { id: 'ds1', name: 'A', n_obs: 10, n_vars: 5, steps: ['qc'] },
  { id: 'ds2', name: 'B', n_obs: 8, n_vars: 5, steps: [] },
]
const DATASETS_R1 = [
  { id: 'ds1', name: 'A', n_obs: 10, n_vars: 5, steps: ['qc', 'de'] },
  { id: 'ds2', name: 'B', n_obs: 8, n_vars: 5, steps: [] },
]

async function setDatasets(engine, datasets, active = 'ds1') {
  const reply = await engine.request('set_datasets', { datasets, active }, { timeout: 15_000 })
  if (!reply || !Array.isArray(reply.datasets)) throw new Error('set_datasets 没有返回数据集')
  return reply
}

// 每个任务自己产生一份内容不同的图（用来验证不同 run 的产物目录互不覆盖）
function makeTask(engine, params, extra = {}) {
  return engine.withTask(
    { op: 'make_figure', dataset: 'ds1', params, timeoutMs: 30_000, ...extra },
    async (ctx) => {
      const src = path.join(TMP, `src-${params.method}-${params.tag || 'x'}.png`)
      writeFileSync(src, `PNG-${params.method}-${params.tag || 'x'}`)
      ctx.step('compute', 'ok', { ms: 1 })
      const artifact = ctx.artifact({ kind: 'figure', source: src, name: 'fig.png' })
      ctx.step('write', 'ok', { artifact: artifact && artifact.relPath })
      return { result: { method: params.method, artifact: artifact && artifact.relPath } }
    },
  )
}

async function main() {
  const engine = createEngine({ workdir: WORK, engineScript: workerPath, python: 'python3', queueGrace: 400 })
  const started = await withTimeout(engine.start(), 25_000, 'start')
  check('P0-0a 引擎握手在时限内结束（上一轮 deferred bug 会让它永久挂住）', started.settled === 'ok',
    started.settled === 'ok' ? 'ok' : String(started.error && started.error.message))
  await setDatasets(engine, DATASETS_R0, 'ds1')

  // -------------------------------------------------------------------------
  // ② 数据版本冻结 / 参考对象 artifact ID
  // -------------------------------------------------------------------------
  const frozen = await engine.freezeDataset()
  const artifactPattern = /^dataset:ds1:[0-9a-f]{16}$/
  check('P0-2a 冻结活动数据集，产出 artifact ID', frozen.id === 'ds1' && artifactPattern.test(frozen.artifactId),
    `${frozen.id} ${frozen.artifactId}`)
  const frozenAgain = await engine.freezeDataset()
  check('P0-2b 数据未变时 revision 与 artifact ID 稳定',
    frozenAgain.revision === frozen.revision && frozenAgain.artifactId === frozen.artifactId,
    frozen.revision.slice(0, 12))
  const resolved = await engine.resolveDatasetRef(frozen.artifactId)
  check('P0-2c artifact ID 可解析回同一 revision',
    resolved.id === 'ds1' && resolved.revision === frozen.revision, artifactPattern.test(frozen.artifactId) ? 'ok' : 'bad id')
  const explicit = await engine.freezeDataset('ds2')
  check('P0-2d 显式 id 冻结不跟随活动数据集', explicit.id === 'ds2' && explicit.artifactId.startsWith('dataset:ds2:'))
  const stillActive = await engine.freezeDataset()
  check('P0-2e 活动数据集仍是 ds1（没有被上一步改掉）', stillActive.id === 'ds1')
  const artifactForTask = frozen.artifactId

  // 数据 revision 变化（模拟用户改了分组/步骤/重载）：旧 artifact ID 必须失效，不能漂移
  await setDatasets(engine, DATASETS_R1, 'ds2')
  let stale = null
  try { await engine.resolveDatasetRef(artifactForTask) } catch (error) { stale = error }
  check('P0-2f 数据 revision 变化后旧 artifact ID → stale_artifact（绝不静默复用）',
    Boolean(stale) && stale.code === 'stale_artifact', stale ? stale.code : 'no-throw')
  check('P0-2g stale_artifact 带上 expected/actual 供诊断',
    Boolean(stale && stale.details && stale.details.expected && stale.details.actual),
    stale && stale.details ? `${String(stale.details.expected).slice(0, 8)}≠${String(stale.details.actual).slice(0, 8)}` : 'no details')
  await setDatasets(engine, DATASETS_R0, 'ds1')

  // -------------------------------------------------------------------------
  // ③ 不可变 run 目录 + 幂等键（输入/参数版本）
  // -------------------------------------------------------------------------
  const run1 = await makeTask(engine, { method: 'a' })
  check('P0-3a 首次任务提交为 runs/<run_id> 独立目录', run1.reused === false && existsSync(run1.run.dir) &&
    readdirSync(run1.run.dir).includes('manifest.json'), run1.runId)
  const manifest1 = JSON.parse(readFileSync(path.join(run1.run.dir, 'manifest.json'), 'utf8'))
  check('P0-3b manifest 固化 run_id / 数据 revision / 参数 / 逐步状态 / 产物哈希',
    manifest1.run_id === run1.runId &&
      manifest1.dataset && manifest1.dataset.revision === frozen.revision &&
      Array.isArray(manifest1.steps) && manifest1.steps.length === 2 &&
      manifest1.artifacts.length === 1 && /^[0-9a-f]{64}$/.test(manifest1.artifacts[0].sha256),
    `steps=${manifest1.steps.length}`)
  check('P0-3c 产物按内容 hash 校验通过', engine.verifyRun(run1.runId).ok === true)
  const record1 = engine.getRun(run1.runId)
  check('P0-3c2 逐步状态持久化（runs/attempts/steps 落 SQLite）',
    engine.store.kind === 'sqlite' && Array.isArray(record1.steps) && record1.steps.length === 2 &&
      record1.steps.every((item) => item.status === 'ok') && record1.steps[0].step === 'compute',
    `store=${engine.store.kind} steps=${(record1.steps || []).map((item) => `${item.step}:${item.status}`).join(',')}`)
  check('P0-3c3 status() 暴露队列与 run store / 下载缓存',
    Boolean(engine.status().runs && engine.status().runs.kind === 'sqlite' && engine.status().queue) &&
      Boolean(engine.status().downloads), JSON.stringify(engine.status().runs))
  const run1File = path.join(run1.run.dir, manifest1.artifacts[0].relPath)
  check('P0-3d 产物落在 run 目录的 figures/ 下且内容正确',
    manifest1.artifacts[0].relPath.startsWith('figures/') && readFileSync(run1File, 'utf8') === 'PNG-a-x',
    manifest1.artifacts[0].relPath)

  const run1again = await makeTask(engine, { method: 'a' })
  check('P0-3e 同参数重复提交幂等（复用同一 run_id，不重跑）',
    run1again.reused === true && run1again.runId === run1.runId,
    run1again.reused ? 'reused' : 'rerun')

  const run2 = await makeTask(engine, { method: 'b' })
  const manifest2 = JSON.parse(readFileSync(path.join(run2.run.dir, 'manifest.json'), 'utf8'))
  check('P0-3f 不同参数 → 新 run_id / 新目录',
    run2.runId !== run1.runId && run2.run.dir !== run1.run.dir, run2.runId)
  check('P0-3g 新 run 不能覆盖旧 run 的同名产物（fig.png 各自独立）',
    readFileSync(run1File, 'utf8') === 'PNG-a-x' &&
      readFileSync(path.join(run2.run.dir, manifest2.artifacts[0].relPath), 'utf8') === 'PNG-b-x' &&
      engine.verifyRun(run1.runId).ok === true)

  // 图/表/日志三类产物分别落 runs/<run_id>/{figures,tables,logs}/
  const mixed = await engine.withTask(
    { op: 'mixed_outputs', dataset: 'ds1', params: { m: 1 }, timeoutMs: 30_000 },
    async (ctx) => {
      const figSrc = path.join(TMP, 'mixed-fig.png')
      const tabSrc = path.join(TMP, 'mixed.tsv')
      writeFileSync(figSrc, 'FIG')
      writeFileSync(tabSrc, 'gene\teffect\nA\t1\n')
      ctx.artifact({ kind: 'figure', source: figSrc, name: 'mixed.png' })
      ctx.artifact({ kind: 'table', source: tabSrc, name: 'de.tsv' })
      ctx.log('run.log', 'hello log')
      ctx.step('done', 'ok')
      return { result: { mixed: true } }
    },
  )
  const mixedManifest = JSON.parse(readFileSync(path.join(mixed.run.dir, 'manifest.json'), 'utf8'))
  const byKind = Object.fromEntries(mixedManifest.artifacts.map((item) => [item.kind, item]))
  check('P0-3q 图/表/日志分别落 runs/<run_id>/{figures,tables,logs}/',
    Boolean(byKind.figure && byKind.figure.relPath.startsWith('figures/') &&
      byKind.table && byKind.table.relPath.startsWith('tables/') &&
      byKind.log && byKind.log.relPath.startsWith('logs/')),
    mixedManifest.artifacts.map((item) => item.relPath).join(' | '))
  check('P0-3r 图/表/日志都登记 sha256 且按内容校验通过',
    ['figure', 'table', 'log'].every((kind) => byKind[kind] && /^[0-9a-f]{64}$/.test(byKind[kind].sha256)) &&
      engine.verifyRun(mixed.runId).ok === true,
    `${Object.keys(byKind).join(',')}`)

  const run1v2 = await makeTask(engine, { method: 'a' }, { revise: true })
  check('P0-3h 同键显式修订 → -v2 新版本，旧版本仍完整可回查',
    run1v2.runId !== run1.runId && /-v2$/.test(run1v2.runId) && engine.verifyRun(run1.runId).ok === true,
    run1v2.runId)

  // 数据 revision 变化 → 同参数也必须产生新 run（输入版本进缓存键）
  await setDatasets(engine, DATASETS_R1, 'ds1')
  const run3 = await makeTask(engine, { method: 'c' })
  check('P0-3i 数据 revision 变化 → 产生新 run（不沿用旧 revision 的缓存）',
    run3.runId !== run1.runId && run3.reused === false, run3.runId)
  const record3 = engine.getRun(run3.runId)
  check('P0-3j 数据 revision 变化后 run 记录的是新 revision（不是活动 ds1 的旧值）',
    record3 && record3.dataset && record3.dataset.revision !== frozen.revision,
    record3 && record3.dataset ? `${String(record3.dataset.revision).slice(0, 8)}…` : 'no dataset')
  await setDatasets(engine, DATASETS_R0, 'ds1')

  // -------------------------------------------------------------------------
  // ① 单步失败不得整体 ok:true；失败/超时任务不产生已提交产物
  // -------------------------------------------------------------------------
  let stepError = null
  try {
    await engine.withTask(
      { op: 'fail_step', dataset: 'ds1', params: { x: 1 }, timeoutMs: 30_000 },
      async (ctx) => {
        ctx.step('qc', 'ok')
        ctx.step('de', 'failed', { error: 'boom' })
        return { result: { ok: true } } // 故意返回“成功”，必须被拒绝
      },
    )
  } catch (error) { stepError = error }
  check('P0-1a 单步失败不得整体 ok:true', Boolean(stepError) && stepError.code === 'step_failed',
    stepError ? stepError.code : 'no-throw')
  const failedRecord = stepError && stepError.runId ? engine.getRun(stepError.runId) : null
  check('P0-1b 失败 run 记为 failed、不产生已提交产物目录',
    Boolean(failedRecord) && failedRecord.status === 'failed' && !failedRecord.dir,
    failedRecord ? failedRecord.status : 'no record')
  const stagingDirs = (() => { try { return readdirSync(path.join(WORK, 'staging')) } catch { return [] } })()
  check('P0-1c 失败 run 的 staging 已清理（不残留半成品）',
    !stagingDirs.some((name) => stepError && name.startsWith(stepError.runId)), stagingDirs.length ? stagingDirs.join(',') : 'empty')
  check('P0-1c2 失败 run 的逐步状态仍可回查（qc:ok / de:failed）',
    Boolean(failedRecord) && Array.isArray(failedRecord.steps) &&
      failedRecord.steps.some((item) => item.step === 'de' && item.status === 'failed') &&
      failedRecord.steps.some((item) => item.step === 'qc' && item.status === 'ok'),
    (failedRecord && failedRecord.steps || []).map((item) => `${item.step}:${item.status}`).join(','))

  let timeoutError = null
  try {
    await engine.withTask(
      { op: 'timeout_task', dataset: 'ds1', params: { x: 2 }, timeoutMs: 1500 },
      async (ctx) => {
        // 任务内请求超出剩余预算：worker 要睡 3s，deadline 只给 ~1.5s
        await ctx.request('block', { seconds: 3 }, { timeout: ctx.remaining() })
        return { result: { done: true } }
      },
    )
  } catch (error) { timeoutError = error }
  check('P0-1d 任务内请求超时 → 整个任务失败（不提交）',
    Boolean(timeoutError) && ['timeout', 'cancelled'].includes(timeoutError.code),
    timeoutError ? timeoutError.code : 'no-throw')
  const timeoutRecord = timeoutError && timeoutError.runId ? engine.getRun(timeoutError.runId) : null
  check('P0-1e 超时 run 状态为 cancelled、无已提交目录',
    Boolean(timeoutRecord) && timeoutRecord.status === 'cancelled' && !timeoutRecord.dir,
    timeoutRecord ? timeoutRecord.status : 'no record')
  // 等 worker 睡完，避免污染后续请求
  await engine.request('list', {}, { timeout: 20_000 })

  // -------------------------------------------------------------------------
  // ② 下载按哈希复用
  // -------------------------------------------------------------------------
  const downloadSrc = path.join(TMP, 'download-src.tsv')
  writeFileSync(downloadSrc, 'gene\tsample\nA\t1\n')
  const dl1 = engine.acquireFile(downloadSrc, { name: 'matrix.tsv' })
  const dl2 = engine.acquireFile(downloadSrc, { name: 'matrix.tsv' })
  check('P0-2h 同一文件第二次进入按内容哈希复用（不重复复制）',
    dl1.reused === false && dl2.reused === true && dl1.path === dl2.path && dl1.sha256 === dl2.sha256,
    `${dl1.sha256.slice(0, 12)}…`)
  writeFileSync(downloadSrc, 'gene\tsample\nA\t1\n')
  const dl3 = engine.acquireFile(downloadSrc, { name: 'another-name.tsv' })
  check('P0-2i 内容相同（文件名不同）→ 复用同一缓存对象', dl3.reused === true && dl3.path === dl1.path)
  let mismatch = null
  try { engine.acquireFile(downloadSrc, { name: 'x.tsv', sha256: 'deadbeef'.repeat(8) }) } catch (error) { mismatch = error }
  check('P0-2j 期望哈希不符 → hash_mismatch（防止下错文件）',
    Boolean(mismatch) && mismatch.code === 'hash_mismatch', mismatch ? mismatch.code : 'no-throw')
  check('P0-2k 下载缓存自校验通过', engine.downloadCache.verify().ok === true,
    `entries=${engine.downloadCache.stats().count}`)

  // 输入清单进缓存键：同参数不同输入 → 新 run；再提交同输入 → 复用
  const makeInputTask = (sha) => engine.withTask(
    { op: 'input_task', dataset: 'ds1', params: { same: true }, inputs: [{ role: 'matrix', sha256: sha }], timeoutMs: 20_000 },
    async (ctx) => { ctx.step('consume', 'ok'); return { result: { sha } } },
  )
  const inputRun1 = await makeInputTask(dl1.sha256)
  const inputRun2 = await makeInputTask(dl1.sha256)
  check('P0-3k 同参数同输入重复提交 → 复用', inputRun2.runId === inputRun1.runId && inputRun2.reused === true)
  const otherSha = 'f'.repeat(64)
  const inputRun3 = await makeInputTask(otherSha)
  check('P0-3l 输入版本变化 → 同参数也产生新 run（输入进缓存键）',
    inputRun3.runId !== inputRun1.runId && inputRun3.reused === false)

  // -------------------------------------------------------------------------
  // ①/③ kill 不污染已提交结果；重启后仍可回查 + 幂等复用
  // -------------------------------------------------------------------------
  const keeper = await makeTask(engine, { method: 'keeper' })
  const keeperManifest = JSON.parse(readFileSync(path.join(keeper.run.dir, 'manifest.json'), 'utf8'))
  const keeperFile = path.join(keeper.run.dir, keeperManifest.artifacts[0].relPath)
  const keeperSha = keeperManifest.artifacts[0].sha256
  const keeperBytes = statSync(keeperFile).size

  // 在一个任务执行中途让 worker 直接死掉：这个任务必须失败，已提交的 keeper 不能受影响。
  // （队列是串行的，排在后面的 kill 命令追不上正在运行的请求，所以用 worker 自毁来模拟。）
  const doomed = withTimeout((async () => {
    try {
      const run = await engine.withTask(
        { op: 'doomed_task', dataset: 'ds1', params: { x: 9 }, timeoutMs: 30_000 },
        async (ctx) => {
          ctx.step('compute', 'ok')
          await ctx.request('selfdestruct', { seconds: 1 }, { timeout: 25_000 })
          ctx.step('write', 'failed', { error: 'worker killed' })
          return { result: { shouldNotCommit: true } }
        },
      )
      return { committed: true, run }
    } catch (error) { return { error } }
  })(), 30_000, 'doomed')
  const doomedOutcome = await doomed
  check('P0-1f worker 崩溃时在途任务在时限内明确结束（不会永久挂住）',
    doomedOutcome.settled !== 'hung', doomedOutcome.settled)
  const doomedValue = doomedOutcome.value || {}
  check('P0-1g 被 kill 打断的任务失败，且不会被提交',
    Boolean(doomedValue.error) && doomedValue.committed !== true,
    doomedValue.error ? String(doomedValue.error.code) : 'committed?!')
  const doomedRunId = doomedValue.error && doomedValue.error.runId
  const doomedRecord = doomedRunId ? engine.getRun(doomedRunId) : null
  check('P0-1h 被 kill 的任务记录为失败、无已提交目录',
    Boolean(doomedRecord) && doomedRecord.status === 'failed' && !doomedRecord.dir,
    doomedRecord ? doomedRecord.status : 'no record')
  check('P0-1i kill 后已提交 run 仍 committed，产物内容与 hash 未变',
    engine.getRun(keeper.runId).status === 'committed' &&
      engine.verifyRun(keeper.runId).ok === true &&
      existsSync(keeperFile) && statSync(keeperFile).size === keeperBytes,
    `${keeperSha.slice(0, 12)}…`)

  await engine.dispose()

  // 重启（新 engine / 新 worker）：SQLite 与 runs/ 目录都还在
  const engine2 = createEngine({ workdir: WORK, engineScript: workerPath, python: 'python3', queueGrace: 400 })
  const restarted = await withTimeout(engine2.start(), 25_000, 'restart')
  check('P0-3m 重启后引擎可用', restarted.settled === 'ok')
  const reloaded = engine2.getRun(run1.runId)
  check('P0-3n 重启后按 run_id 仍能回查历史 run 与产物',
    Boolean(reloaded) && reloaded.status === 'committed' && reloaded.artifacts.length === 1,
    reloaded ? reloaded.runId : 'missing')
  check('P0-3o 重启后产物内容 hash 仍校验通过（不可变）',
    engine2.verifyRun(keeper.runId).ok === true && engine2.verifyRun(run1.runId).ok === true)
  let reuseAfterRestart = null
  try {
    reuseAfterRestart = await makeTask(engine2, { method: 'a' })
  } catch (error) { reuseAfterRestart = { error } }
  const reuseKeyRunIds = [run1.runId, run1v2.runId]
  check('P0-3p 重启后同参数同数据 → 复用已提交 run，不重跑',
    Boolean(reuseAfterRestart) && reuseAfterRestart.reused === true && reuseKeyRunIds.includes(reuseAfterRestart.runId),
    reuseAfterRestart && reuseAfterRestart.runId ? reuseAfterRestart.runId : String(reuseAfterRestart && reuseAfterRestart.error))
  await engine2.dispose()

  console.log(`\n== 结果：${passed} 通过 / ${failed} 失败 ==`)
  if (failed) {
    for (const item of failures) console.log(`  ✗ ${item}`)
  }
  return failed === 0
}

const WATCHDOG_MS = 180_000
const watchdog = setTimeout(() => {
  console.error(`\n测试挂起超过 ${WATCHDOG_MS / 1000}s，强制退出（防止出现不退出测试）`)
  try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }
  process.exit(2)
}, WATCHDOG_MS)

main()
  .then((ok) => {
    clearTimeout(watchdog)
    try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }
    process.exit(ok ? 0 : 1)
  })
  .catch((error) => {
    clearTimeout(watchdog)
    console.error('测试本身异常：', error)
    try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }
    process.exit(1)
  })
