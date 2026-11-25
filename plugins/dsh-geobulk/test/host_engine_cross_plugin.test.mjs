/**
 * dsh-* ×6 — 跨插件引擎冒烟 + 共享模块一致性（任务单 §7 第 1 条，六份实现同源）
 *
 * 上一轮改写把六个插件的 lib/engine.mjs 都换成了「TaskQueue + RunStore」版本，但
 * sendRequest 的错误实现让**每一个请求（包括 start() 的 ping 握手）都永久挂住**。
 * 本测试对六个插件逐一复验：
 *   1. 握手 + 业务请求在时限内结束（死锁回归）；
 *   2. 排队期间超时 → 该请求明确失败，且 worker 从未收到它（零执行）；
 *   3. 共享模块 q_engine_queue.mjs / q_run_store.mjs 六份逐字相同；
 *   4. 每份 engine.mjs 都带本轮修复标记（deferred 修复、queueGrace、artifact ID、步骤失败）。
 *
 * 运行：node test/host_engine_cross_plugin.test.mjs
 */

import { mkdtempSync, writeFileSync, readFileSync, rmSync, existsSync, mkdirSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

const HERE = path.dirname(fileURLToPath(import.meta.url))
const PLUGINS_ROOT = path.resolve(HERE, '..', '..')
const TMP = mkdtempSync(path.join(tmpdir(), 'dsh-engine-cross-'))
const PLUGINS = ['geobulk', 'singlecell', 'spatial', 'metabolomics', 'plotlab', 'virtualcell']

const FAKE_WORKER = `#!/usr/bin/env python3
import json, os, sys

LOG = os.environ.get("DSH_TEST_WORKER_LOG")

def reply(msg_id, result=None, error=None):
    payload = {"id": msg_id, "ok": error is None}
    if error is None:
        payload["result"] = result if result is not None else {}
    else:
        payload["error"] = {"code": "internal", "message": error}
    sys.stdout.write(json.dumps(payload) + "\\n")
    sys.stdout.flush()

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
    if LOG:
        with open(LOG, "a", encoding="utf-8") as handle:
            handle.write(cmd + "|" + json.dumps(args) + "\\n")
            handle.flush()
    if cmd == "ping":
        reply(msg.get("id"), {"mode": "builtin", "packages": {}, "layouts": ["umap"], "platform": "test"})
    elif cmd == "echo":
        reply(msg.get("id"), {"echo": args.get("value"), "plugin": args.get("plugin")})
    elif cmd == "block":
        import time
        time.sleep(float(args.get("seconds") or 1.0))
        reply(msg.get("id"), {"blocked": True})
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

function logCommands(logFile) {
  try {
    return readFileSync(logFile, 'utf8').split('\n').filter(Boolean).map((line) => line.split('|')[0])
  } catch { return [] }
}

const FIX_MARKERS = [
  ['deferred 死锁修复', 'if (ok) resolve(value)'],
  ['queueGrace 可配置', 'const queueGrace ='],
  ['taskId 唯一', '++state.taskSeq'],
  ['artifact ID 参考对象', 'resolveDatasetRef'],
  ['下载哈希复用', 'openDownloadCache'],
  ['步骤失败不得成功', "code: 'step_failed'"],
  ['不可变 run 校验', 'verifyRun'],
]

async function checkOnePlugin(plugin) {
  const enginePath = path.join(PLUGINS_ROOT, `dsh-${plugin}`, 'lib', 'engine.mjs')
  if (!existsSync(enginePath)) { check(`X ${plugin}: engine.mjs 存在`, false, enginePath); return }
  const mod = await import(`file://${enginePath}`)
  check(`X-${plugin}-a createEngine 可导入`, typeof mod.createEngine === 'function')

  const logFile = path.join(TMP, `${plugin}.log`)
  const workdir = path.join(TMP, `work-${plugin}`)
  mkdirSync(workdir, { recursive: true })
  process.env.DSH_TEST_WORKER_LOG = logFile
  const engine = mod.createEngine({
    workdir,
    engineScript: workerPath,
    python: process.env.DSH_TEST_PYTHON || 'python3',
    queueGrace: 400,
  })

  const started = await withTimeout(engine.start(), 30_000, `${plugin}-start`)
  check(`X-${plugin}-b 握手在时限内结束（deferred 死锁回归）`, started.settled === 'ok',
    started.settled === 'ok' ? 'ok' : String(started.error && started.error.message).slice(0, 80))

  const roundtrip = await withTimeout(engine.request('echo', { value: 42, plugin }, { timeout: 15_000 }), 20_000, `${plugin}-echo`)
  check(`X-${plugin}-c 业务请求正常往返`, roundtrip.settled === 'ok' && roundtrip.value.echo === 42,
    roundtrip.settled === 'ok' ? `echo=${roundtrip.value.echo}` : String(roundtrip.error && roundtrip.error.code))

  // 排队超时 → 不派发
  const blocker = withTimeout(engine.request('block', { seconds: 0.8 }, { timeout: 15_000 }), 20_000, `${plugin}-block`)
  await sleep(120)
  const mark = `TIMEOUT-${plugin}`
  const impatient = await withTimeout(engine.request('echo', { value: mark, plugin }, { timeout: 100 }), 15_000, `${plugin}-impatient`)
  check(`X-${plugin}-d 排队超时的请求按 timeout 失败`,
    impatient.settled === 'error' && impatient.error && impatient.error.code === 'timeout',
    impatient.settled === 'error' ? String(impatient.error.code) : impatient.settled)
  await blocker
  await sleep(400)
  const cmds = logCommands(logFile)
  const dispatched = [] // 从 worker 日志里找带 mark 的 echo
  try {
    const lines = readFileSync(logFile, 'utf8').split('\n').filter(Boolean)
    for (const line of lines) if (line.startsWith('echo|') && line.includes(mark)) dispatched.push(line)
  } catch { /* ignore */ }
  check(`X-${plugin}-e 超时的排队项从未被执行（worker 零收到）`, dispatched.length === 0,
    `echo 命令 ${cmds.filter((cmd) => cmd === 'echo').length} 次`)

  await engine.dispose()
  delete process.env.DSH_TEST_WORKER_LOG
  check(`X-${plugin}-f dispose 后 running=false`, engine.running === false)
}

async function main() {
  // 1) 共享模块一致性（六份必须逐字相同）
  const queueHashes = new Set()
  const storeHashes = new Set()
  for (const plugin of PLUGINS) {
    const q = path.join(PLUGINS_ROOT, `dsh-${plugin}`, 'lib', 'q_engine_queue.mjs')
    const s = path.join(PLUGINS_ROOT, `dsh-${plugin}`, 'lib', 'q_run_store.mjs')
    queueHashes.add(readFileSync(q, 'utf8'))
    storeHashes.add(readFileSync(s, 'utf8'))
  }
  check('X-share-1 q_engine_queue.mjs 六份逐字相同', queueHashes.size === 1, `${queueHashes.size} 个版本`)
  check('X-share-2 q_run_store.mjs 六份逐字相同', storeHashes.size === 1, `${storeHashes.size} 个版本`)

  // 2) 每份 engine.mjs 都有本轮修复标记
  for (const plugin of PLUGINS) {
    const text = readFileSync(path.join(PLUGINS_ROOT, `dsh-${plugin}`, 'lib', 'engine.mjs'), 'utf8')
    const missing = FIX_MARKERS.filter(([, marker]) => !text.includes(marker)).map(([name]) => name)
    check(`X-${plugin}-marker 修复标记齐备`, missing.length === 0, missing.length ? `缺：${missing.join('、')}` : 'all')
  }

  // 3) 六个引擎逐个冒烟（真实 engine.mjs + 假 worker）
  for (const plugin of PLUGINS) {
    await checkOnePlugin(plugin)
  }

  console.log(`\n== 结果：${passed} 通过 / ${failed} 失败 ==`)
  if (failed) for (const item of failures) console.log(`  ✗ ${item}`)
  return failed === 0
}

const WATCHDOG_MS = 240_000
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
