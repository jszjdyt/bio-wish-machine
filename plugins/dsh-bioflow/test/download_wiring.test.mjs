/**
 * 下载编排回归（Windows 实战踩坑修复）：
 *   1. geo 下载参数必须按引擎契约传 dest_dir（曾因传 destDir 驼峰被引擎忽略，
 *      文件下进引擎默认目录，bioflow 扫 destDir 一无所获）；
 *   2. files 参数（数组或 CLI 逗号串）必须透传；
 *   3. 插件引擎没落任何文件 → 必须 ok:false 报 no_files，绝不静默成功；
 *   4. 引擎自身 ok:false → 必须透传错误码。
 * 全部离线：pluginDownload 用桩，不联网、不 spawn python。
 */
import { test } from 'node:test'
import assert from 'node:assert/strict'
import { mkdtempSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'

import * as S from '../lib/search.mjs'
import { normalizeFileList } from '../lib/impl.mjs'

test('geo 下载：dest_dir（蛇形）与 files 按引擎契约透传', async () => {
  const resolved = await S.resolveDownloadTasks({
    source: 'geo', id: 'gse246662', sourcesClient: null,
    pluginDownload: async () => ({}),
    destDir: '/tmp/dl', files: ['GSE246662_RAW.tar'],
  })
  assert.equal(resolved.mode, 'plugin')
  assert.equal(resolved.params.gse, 'GSE246662')
  assert.equal(resolved.params.dest_dir, '/tmp/dl', '必须传 dest_dir 而不是 destDir')
  assert.equal(resolved.params.destDir, undefined)
  assert.deepEqual(resolved.params.files, ['GSE246662_RAW.tar'])
})

test('files 归一化：数组原样、逗号串拆数组、空值 → null', () => {
  assert.deepEqual(normalizeFileList(['a.txt', ' b.txt ']), ['a.txt', 'b.txt'])
  assert.deepEqual(normalizeFileList('a.txt, b.txt'), ['a.txt', 'b.txt'])
  assert.equal(normalizeFileList(undefined), null)
  assert.equal(normalizeFileList(null), null)
})

test('插件引擎没落任何文件 → ok:false no_files（不许静默成功）', async () => {
  const root = mkdtempSync(path.join(tmpdir(), 'bf-dlwire-'))
  try {
    const result = await S.executeDownload({
      store: null, downloader: null, sourcesClient: null,
      source: 'geo', id: 'GSE246662', destDir: root,
      pluginDownload: async () => ({ ok: true, files: [] }),
    })
    assert.equal(result.ok, false)
    assert.equal(result.error.code, 'no_files')
    assert.match(result.error.message, /没有任何文件/)
  } finally { rmSync(root, { recursive: true, force: true, maxRetries: 5, retryDelay: 100 }) }
})

test('插件引擎 ok:false → 错误码透传', async () => {
  const root = mkdtempSync(path.join(tmpdir(), 'bf-dlwire2-'))
  try {
    const result = await S.executeDownload({
      store: null, downloader: null, sourcesClient: null,
      source: 'geo', id: 'GSE246662', destDir: root,
      pluginDownload: async () => ({ ok: false, error: { code: 'not_found', message: '模拟引擎失败' } }),
    })
    assert.equal(result.ok, false)
    assert.equal(result.error.code, 'not_found')
  } finally { rmSync(root, { recursive: true, force: true, maxRetries: 5, retryDelay: 100 }) }
})
