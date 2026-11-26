# 仓库与更新说明

本仓库（`jszjdyt/bio-wish-machine`）是「万能生信许愿机」的**可安装合集**：
包含统一入口 `dsh-bioflow` 与 7 个生信工作台的源码快照、`install.sh` / `uninstall.sh`。

- 主开发仓库（最新代码 / issue / 版本）：<https://github.com/jszjdyt/dsh-bioflow>
- 本合集仓库：<https://github.com/jszjdyt/bio-wish-machine>
- 当前版本：**0.2.0**（2026-11-26，评审 P0 修复与许愿交互重做）

## 从主仓库同步到合集

```bash
# 1) 拉主仓库最新代码
git clone https://github.com/jszjdyt/dsh-bioflow.git /tmp/dsh-bioflow
# 2) 同步到合集内的插件目录（保留合集自己的 README/INSTALL ）
rsync -a --exclude='.git' --exclude='node_modules' --exclude='.bioflow-engine' \
  /tmp/dsh-bioflow/ /path/to/bio-wish-machine/plugins/dsh-bioflow/
# 3) 提交并推送
cd /path/to/bio-wish-machine
git add -A && git commit -m "chore: 同步 dsh-bioflow 最新版本"
git push
```

## 安装

见 `README.md`：把整个文件夹放到固定位置 → `bash install.sh` → 完全退出 DSH 再打开。
更新时重复上面「同步 → 安装」即可；数据目录 `~/.dsh/bioflow` 不会被覆盖。
