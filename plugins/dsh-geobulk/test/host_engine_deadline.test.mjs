/**
 * dsh-geobulk / dsh-singlecell / dsh-spatial / dsh-metabolomics / dsh-plotlab / dsh-virtualcell
 * — 队列项超时/取消后不得执行（任务单 7.1，HOST-06）
 *
 * 复现的问题：request() 的 guard 只 reject 前端 Promise；排队中的任务在前序请求释放后
 * 仍会照常被 sendRequest 派发，原写入照跑（实测 GeoBulk：第二个请求已经对用户超时返回，
 * 几秒后 worker 仍然收到并执行了它）。
 *
 * 修复后必须满足：
 *   1. 排队期间 deadline 到期 → 该项永不派发（零写入）；
 *   2. 排队期间被取消（AbortSignal）→ 永不派发；
 *   3. 运行中被取消 → 只影响该任务，不杀 worker、不影响后续请求（任务隔离）；
 *   4. 和 host_engine_queue.test.mjs 的 HOST-02 语义兼容（ping 握手不排队、崩溃后队列不互相等待）。
 *
 * 运行：node test/host_engine_deadline.test.mjs
 */

import { mkdtempSync, writeFileSync, readFileSync, existsSync, statSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { createEngine } from '../lib/engine.mjs'

const TMP = mkdtempSync(path.join(tmpdir(), 'dsh-engine-deadline-'))
const LOG = path.join(TMP, 'received.log')
const WRITES = path.join(TMP, 'writes.log')

const FAKE_WORKER = `#!/usr/bin/env python3
import json, os, sys, time

LOG = ${JSON.stringify(LOG)}
WRITES = ${JSON.stringify(WRITES)}

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
    elif cmd == "block":
        time.sleep(float(args.get("seconds") or 1.0))
        reply(msg.get("id"), {"blocked": True})
    elif cmd == "write":
        # 这是「原写入」：只有真的被派发才会发生
        with open(WRITES, "a", encoding="utf-8") as handle:
            handle.write(str(args.get("mark") or "?") + "\\n")
            handle.flush()
        reply(msg.get("id"), {"written": True, "mark": args.get("mark")})
    elif cmd == "list":
        reply(msg.get("id"), {"datasets": [], "active": None})
    elif cmd == "shutdown":
        reply(msg.get("id"), {"bye": True})
        break
    else:
        reply(msg.get("id"), {})
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

function receivedCommands() {
  try {
    return readFileSync(LOG, 'utf8').split('\n').filter(Boolean).map((line) => line.split('|')[0])
  } catch {
    return []
  }
}

function writes() {
  try {
    return readFileSync(WRITES, 'utf8').split('\n').filter(Boolean)
  } catch {
    return []
  }
}

function writesSize() {
  try { return statSync(WRITES).size } catch { return 0 }
}

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms))

async function settled(promise) {
  return promise.then((value) => ({ ok: true, value })).catch((error) => ({ ok: false, error }))
}

async function main() {
  // queueGrace 调小，让「排队超时」在测试里可观测（真实默认 30s）
  const engine = createEngine({
    workdir: path.join(TMP, 'work'),
    engineScript: workerPath,
    python: process.env.DSH_TEST_PYTHON || 'python3',
    queueGrace: 400,
  })
  await engine.start()
  check('HOST-06a 握手 ping 不经过业务队列', receivedCommands().filter((cmd) => cmd === 'ping').length === 1,
    `ping=${receivedCommands().filter((cmd) => cmd === 'ping').length}`)

  // ---- 1. 排队超时 → 永不派发（核心修复点） ----
  const blocker = settled(engine.request('block', { seconds: 1.5 }, { timeout: 10_000 }))
  await sleep(150)  // 让 block 真正执行、占住队列
  const beforeSize = writesSize()
  const impatient = await settled(engine.request('write', { mark: 'TIMEOUT-ITEM' }, { timeout: 100 }))
  check('HOST-06b 排队超时的请求按 timeout 失败', impatient.ok === false && impatient.error && impatient.error.code === 'timeout',
    impatient.ok ? 'unexpected ok' : String(impatient.error && impatient.error.code))
  await blocker
  await sleep(400)  // 给旧实现「前序释放后执行原写入」的机会
  const cmdsAfter = receivedCommands()
  const writesAfterTimeout = writes().filter((mark) => mark === 'TIMEOUT-ITEM')
  check('HOST-06c 超时的排队项没有被派发（worker 从未收到 write）',
    !cmdsAfter.includes('write'), `write 出现 ${cmdsAfter.filter((cmd) => cmd === 'write').length} 次`)
  check('HOST-06d 超时的排队项零写入（原写入不执行）',
    writesAfterTimeout.length === 0 && writesSize() === beforeSize,
    `writes=${JSON.stringify(writes())}`)

  // ---- 2. 排队期间取消 → 永不派发 ----
  const blocker2 = settled(engine.request('block', { seconds: 1.2 }, { timeout: 10_000 }))
  await sleep(150)
  const controller = new AbortController()
  const cancellable = settled(engine.request('write', { mark: 'ABORTED-QUEUED' }, { timeout: 10_000, signal: controller.signal }))
  await sleep(120)
  controller.abort()
  const cancelled = await cancellable
  check('HOST-06e 排队中被取消的请求以 cancelled/timeout 明确失败',
    cancelled.ok === false && ['cancelled', 'aborted', 'timeout'].includes(cancelled.error && cancelled.error.code),
    cancelled.ok ? 'unexpected ok' : String(cancelled.error && cancelled.error.code))
  await blocker2
  await sleep(400)
  check('HOST-06f 被取消的排队项没有被派发、零写入',
    !receivedCommands().includes('write') && writes().filter((mark) => mark === 'ABORTED-QUEUED').length === 0,
    `writes=${JSON.stringify(writes())}`)

  // ---- 3. 运行中取消：任务隔离，不杀 worker ----
  const inflight = new AbortController()
  const running = settled(engine.request('block', { seconds: 0.9 }, { timeout: 10_000, signal: inflight.signal }))
  await sleep(200)
  inflight.abort()
  const runCancel = await running
  check('HOST-06g 运行中取消只让该任务失败',
    runCancel.ok === false && ['cancelled', 'aborted'].includes(runCancel.error && runCancel.error.code),
    runCancel.ok ? 'unexpected ok' : String(runCancel.error && runCancel.error.code))
  check('HOST-06h 运行中取消不杀 worker', engine.running === true, `running=${engine.running}`)
  const afterCancel = await settled(engine.request('list', {}, { timeout: 20_000 }))
  check('HOST-06i 取消之后后续请求仍可执行', afterCancel.ok === true,
    afterCancel.ok ? 'ok' : String(afterCancel.error && afterCancel.error.code))
  const pendingAfter = engine._internal ? engine._internal.pending.size : -1
  check('HOST-06j 取消/超时后没有残留 pending', pendingAfter === 0, `pending=${pendingAfter}`)

  await engine.dispose()
  try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }

  console.log(`\n== 结果：${passed} 通过 / ${failed} 失败 ==`)
  if (failed) {
    for (const item of failures) console.log(`  ✗ ${item}`)
    process.exit(1)
  }
}

main().catch((error) => {
  console.error('测试本身异常：', error)
  process.exit(1)
})