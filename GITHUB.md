# 推到 GitHub（三步）

本仓库已经初始化好并完成首次提交（`main` 分支、479 个文件、commit `468f11f`）。
这台机器上**没有配置 GitHub 凭据**（`gh` CLI 未安装、没有全局 git 身份），所以「最后一脚」需要你用自己的账号完成，二选一即可。

---

## 路线 A：用网页建仓 + 命令行推送（推荐，最简单）

1. 浏览器打开 <https://github.com/new>
2. 仓库名填 **万能生信许愿机**
   - 如果觉得中文 URL 不方便，可以改用 `bio-wish-machine`，Description 里写「万能生信许愿机 · 许愿式生物信息全自动分析平台」即可。
3. **不要**勾选 "Add a README / .gitignore / license"（否则会有冲突）。选好 Public / Private，点 Create。
4. 回到终端，把 `<你的用户名>` 和仓库名替换掉：

```bash
cd ~/Desktop/万能生信许愿机
git remote add origin https://github.com/<你的用户名>/万能生信许愿机.git
git push -u origin main
```

第一次推送会要求认证：GitHub 从 2021 年起**不支持账号密码**，密码框要填 **Personal Access Token**：

- 生成入口：<https://github.com/settings/tokens> → Generate new token (classic) → 勾 `repo` → 生成后复制
- 用户名填你的 GitHub 用户名，密码位置粘贴 token

---

## 路线 B：装 gh CLI 一条命令搞定

```bash
brew install gh
gh auth login                     # 按提示浏览器授权
cd ~/Desktop/万能生信许愿机
gh repo create 万能生信许愿机 --public --source=. --remote=origin --push
```

---

## 推完之后建议做两件事

1. **打 tag 和 Release**（方便师弟直接下 zip）：
   ```bash
   git tag -a v0.1.0 -m "万能生信许愿机 v0.1.0"
   git push origin v0.1.0
   ```
   然后在 GitHub 仓库页 → Releases → Draft a new release → 选 `v0.1.0`，把桌面上的
   `万能生信许愿机-mac.zip`（14 MB）作为附件传上去。

2. **写仓库简介**（右侧 About 齿轮）：
   `许愿式生物信息全自动分析平台：一句话 → 检索/下载/立项/分析/出图/留证，包含 GEO Bulk、单细胞、空间、代谢组、绘图、虚拟细胞 7 个工作台`

---

## 备注

- 桌面上的 `万能生信许愿机-mac.zip` 是**给师弟的安装包**（14 MB，sha256 前 24 位 `af8b1add33cedb12608e206e`），里面不含 `.git`，解压后 `bash install.sh` 即可。
- 仓库里 `.gitignore` 已排除 `runs/`、`node_modules/`、`__pycache__/`、`.zip`，不会把运行产物或安装包重复提交。
- 这个仓库是**独立仓库**；原工作区里每个插件各自也有自己的 `.git`，两者互不影响。
