#!/usr/bin/env bash
# 万能生信许愿机 · 卸载脚本
# 只清理本包登记的依赖与符号链接，不动你的数据、不动其它插件。

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLUGINS="$HERE/plugins"

PROFILE=""
if [ -f "$HOME/.dsh/profiles/desktop/package.json" ]; then
  PROFILE="$HOME/.dsh/profiles/desktop"
else
  for cand in "$HOME"/.dsh/profiles/*/; do
    if [ -f "$cand/package.json" ]; then PROFILE="${cand%/}"; break; fi
  done
fi
if [ -z "$PROFILE" ]; then
  echo "✗ 没找到 DSH profile" >&2
  exit 1
fi

cp "$PROFILE/package.json" "$PROFILE/package.json.bak-unwish-$(date +%Y%m%d-%H%M%S)"

node - "$PROFILE" "$PLUGINS" <<'NODE'
const fs = require('fs')
const path = require('path')
const [profile, plugins] = process.argv.slice(2)
const pkgPath = path.join(profile, 'package.json')
const pkg = JSON.parse(fs.readFileSync(pkgPath, 'utf8'))
pkg.dependencies = pkg.dependencies || {}
const nm = path.join(profile, 'node_modules')
let n = 0
for (const name of fs.readdirSync(plugins)) {
  const dir = path.join(plugins, name)
  if (!fs.existsSync(path.join(dir, 'package.json'))) continue
  const meta = JSON.parse(fs.readFileSync(path.join(dir, 'package.json'), 'utf8'))
  const id = meta.name || name
  if (pkg.dependencies[id]) { delete pkg.dependencies[id]; n += 1 }
  const link = path.join(nm, id)
  try { fs.rmSync(link, { recursive: true, force: true }) } catch { /* ignore */ }
}
fs.writeFileSync(pkgPath, JSON.stringify(pkg, null, 2) + '\n')
console.log(`已移除 ${n} 个插件登记。`)
NODE

echo "✅ 卸载完成，⌘Q 退出 DSH 再重开即可。"
