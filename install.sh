#!/usr/bin/env bash
# 万能生信许愿机 · macOS 安装脚本
# 作用：把本包里的 8 个 DSH 插件登记到你的 DSH 桌面版 profile，并建立符号链接。
# 幂等：重复运行会覆盖同名链接；每次都会先备份 package.json。

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLUGINS="$HERE/plugins"

if [ ! -d "$PLUGINS" ]; then
  echo "✗ 找不到 plugins/ 目录（请把整个文件夹一起拷贝，不要只拷 install.sh）" >&2
  exit 1
fi

# ---- 1. 定位 DSH profile ----
PROFILE=""
if [ -f "$HOME/.dsh/profiles/desktop/package.json" ]; then
  PROFILE="$HOME/.dsh/profiles/desktop"
else
  for cand in "$HOME"/.dsh/profiles/*/; do
    if [ -f "$cand/package.json" ]; then PROFILE="${cand%/}"; break; fi
  done
fi
if [ -z "$PROFILE" ]; then
  echo "✗ 没找到 DSH profile（$HOME/.dsh/profiles/<名字>/package.json）" >&2
  echo "  请先安装并启动一次 DeepSeek Harness 桌面版，再运行本脚本。" >&2
  exit 1
fi
echo "→ 目标 profile：$PROFILE"

# ---- 2. 备份 ----
BACKUP="$PROFILE/package.json.bak-wish-$(date +%Y%m%d-%H%M%S)"
cp "$PROFILE/package.json" "$BACKUP"
echo "→ 已备份：$BACKUP"

# ---- 3. 登记依赖 + 建符号链接 ----
node - "$PROFILE" "$PLUGINS" <<'NODE'
const fs = require('fs')
const path = require('path')
const [profile, plugins] = process.argv.slice(2)
const pkgPath = path.join(profile, 'package.json')
const pkg = JSON.parse(fs.readFileSync(pkgPath, 'utf8'))
pkg.dependencies = pkg.dependencies || {}
const nm = path.join(profile, 'node_modules')
fs.mkdirSync(nm, { recursive: true })

let n = 0
for (const name of fs.readdirSync(plugins).sort()) {
  const dir = path.join(plugins, name)
  const manifest = path.join(dir, 'package.json')
  if (!fs.existsSync(manifest)) continue
  const meta = JSON.parse(fs.readFileSync(manifest, 'utf8'))
  const id = meta.name || name
  pkg.dependencies[id] = `link:${dir}`
  const link = path.join(nm, id)
  try { fs.rmSync(link, { recursive: true, force: true }) } catch { /* 旧链接可能不存在 */ }
  fs.symlinkSync(dir, link, 'dir')
  console.log(`  ✓ ${id}  →  ${dir}`)
  n += 1
}
fs.writeFileSync(pkgPath, JSON.stringify(pkg, null, 2) + '\n')
console.log(`\n一共登记 ${n} 个插件。`)
NODE

echo
echo "✅ 安装完成。"
echo "   请「⌘Q 完全退出 DSH」，再重新打开。"
echo "   重开后工具表里应能看到：bioflow（许愿机）/ geobulk / singlecell / spatial / metabolomics / plotlab / virtualcell"
echo
echo "   验证：在对话里说一句「用 bioflow 看看状态」，或运行"
echo "         node \"$PLUGINS/dsh-bioflow/bin/bioflow.mjs\" status --json"
