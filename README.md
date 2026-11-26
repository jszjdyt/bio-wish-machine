# 万能生信许愿机

一句话说人话，让 agent 把「找数据 → 下载 → 立项 → 跑分析 → 出图 → 留证」整条链跑完。
包装了 7 个生信工作台 + 1 个统一入口（许愿机本体），全部跑在 DSH 桌面版里。

> 版本：0.2.1 · 平台：macOS（Apple Silicon / Intel 均可）· 需要 DeepSeek Harness 桌面版
> 最新代码：<https://github.com/jszjdyt/dsh-bioflow> · 本仓库为可安装合集（含插件源码快照）

---

## 里面有什么

| 插件 | 干什么 |
|---|---|
| **dsh-bioflow**（许愿式生物信息全自动分析平台） | 统一入口。自然语言许愿 / 疾病引导 / 问题入手 → 自动挑数据（中文摘要弹窗）→ draft 计划（先给提案，人只检查）→ 确认 → 不可变 run → 产物与证据；GEO/PubMed/GDC-TCGA/MetaboLights/CELLxGENE 检索、统一下载、多数据集去批次合并（`ds_harm_*`）、联合分析、删除与导出 |
| dsh-geobulk | GEO Bulk：检索/下载/导入、QC、差异表达（DESeq2/edgeR/limma-voom）、富集、免疫浸润、生存与预后、circRNA/miRNA/lncRNA/融合/剪接 |
| dsh-singlecell | 单细胞：QC、聚类、UMAP、标记基因、注释、供体级 pseudobulk（edgeR QL）、组成分析、轨迹、细胞通讯 |
| dsh-spatial | 空间转录组：Visium/Stereo-seq 导入、空间域、SVG（Moran's I）、邻域富集、反卷积 |
| dsh-metabolomics | 代谢组：LC-MS/GC-MS 特征表、PCA/PLS-DA、差异、通路富集、标志物筛选、MS2 鉴定 |
| dsh-plotlab | 绘图工作台：22 种图、一句话出图、统计检验、300dpi PNG + 可复现代码 |
| dsh-virtualcell | 虚拟细胞：状态表征、GRN、虚拟敲除/过表达、扰动响应预测、状态转移 |
| dsh-nsclc-vt | NSCLC 虚拟试验：组合干预排名（产物标为模型预测，不是实验结论） |

---

## 本轮更新（0.2.1 · 2026-11-26）

依据两份评审（UX review / 产品升级方案）重做了许愿交互，并修掉 4 个会改变结果的高危缺陷：

- **四种起手式**：我有数据 / 我只有方向 / 帮我找数据 / 照着图来（6 张出图样式卡，选图反推所需数据与分组）。
- **活体愿望卡**：绿/黄/红体检灯 + 待定字段 + 每条闸门的中文补救动作；缺字段时**一次只问一个问题**（2–4 个可点候选）。
- **数据选择**：许愿后自动弹「候选数据集（中文摘要）」，写明物种/样本量/模态/分组/矩阵/临床/适合做什么/风险；疾病引导里本地与在线（GEO/TCGA-GDC/MetaboLights）都可直接选，在线选中**自动下载+登记**；本地/在线追加多个数据时自动区分「同类型→先去批次合并」「不同类型→一脉相承」。
- **接着做**：run 成功后自动给「顺着结果继续许愿」（QC/批次、趋势 Top 基因、独立队列验证、TCGA 生存、单细胞表达）；零显著明确写「当前阈值下 0 项显著」。
- **高危修复**：点名未登记的 GSE/TCGA 不再静默用本地数据顶替；矩阵方向与样本身份核对共用同一轴；负整数/小数不再自动标 `counts/verified`；队列×分组完全混杂时 harmonize 前置拦截；合并样本按 `dataset_id:原样本名` 逐队列回填并做样本数守恒校验；单细胞无供体列不再生成伪重复 pseudobulk。
- 中文输入法选词按 Enter 不会误提交；下载按钮语义改为「停止等待」。
- **视觉减负**：修掉 `.bf-progress` 类名冲突导致的巨型圆角块；「进度与产物」改成扁平细线（无 run 时只有一行提示，有 run 才出现内联统计与 4px 细进度条）；顶栏健康灯改成暗色背景上的低调文字「数据源 x/5 · 工作台 y/6」。
- **回收站**：删除的 run 先移入 `trash/runs/`，顶栏/历史 run 面板的「回收站」可查看、**恢复**（重建索引，可继续看结果）或**彻底删除/清空**；对应 `bioflow action=trash trash_op=list|restore|purge|empty` 与 CLI `bioflow trash --trash-op=...`。
- **页面减负**：顶栏健康灯收成一行「数据源 5/5 · 工作台 6/6」；示例/常用套路默认折叠；没有计划/run 时不再显示空白流水线；历史 run 空态压成一行。

回归：`node --test plugins/dsh-bioflow/test/*.test.mjs` → **64 项 / 63 通过 / 1 skip / 0 失败**（新增 `review_p0` 与 `harmonize_merge`）。

---

## 安装

1. 先装好 **DeepSeek Harness 桌面版**，并启动过至少一次（这样才会有 `~/.dsh/profiles/desktop`）。
2. 把整个文件夹放到一个**不会被挪走**的位置，例如 `~/Applications/万能生信许愿机`（脚本用绝对路径登记，之后移动文件夹需要重跑安装）。
3. 打开终端，执行：

```bash
cd ~/Applications/万能生信许愿机
bash install.sh
```

4. **⌘Q 完全退出** DSH，再重新打开。

卸载：`bash uninstall.sh`（只移除插件登记，不动你的数据）。

---

## 环境依赖（Mac 自带的不够，需要补）

- **Node ≥ 18**：DSH 自带运行时即可，无需另装。
- **Python 3.10+**，且装齐：
  ```bash
  pip3 install numpy pandas scipy scikit-learn matplotlib seaborn statsmodels lifelines openpyxl anndata h5py networkx
  ```
- **R 4.4+（可选，但强烈建议）**：供体级 edgeR QL、limma-voom、Seurat 6、timeROC、fgsea 都会用到，缺了会自动降级并**如实记录** `method_actual / fallback`，不会伪造结果。
- 空间插件的 squidpy / 单细胞的 scanpy 是**可选**的：装了用原生实现，没装走内置实现（结果口径一样，速度慢些）。

联网仅用于 GEO / PubMed / GDC / MetaboLights 检索与下载。CELLxGENE 接口在大陆网络通常 60 秒仍超时，插件会明确报 `source_unavailable`，不会卡住整条流程。

---

## 三分钟上手

在 DSH 对话里直接说人话：

```
帮我把 GSE84433 的胃癌肝转移，按 NK 细胞比例跟原发灶比一下，出图
```

许愿机会：
1. 回一份 **draft 计划**（数据源、队列、比较方向、设计、产物），**不会直接开跑**；
2. 缺关键字段（分组、配对、终点）时只问你最少的问题，并给候选；
3. 你确认后它给 `run_id`，按步执行，任一步失败不会整体报成功；
4. 产物在 `~/.dsh/bioflow/runs/<run_id>/`：`manifest.json`、`design.tsv`、`qc.tsv`、`results/*.tsv`、`figures/*.png`、`evidence.jsonl`、`methods.md`。

也可以用命令行（和 agent 工具是同一套实现）：

```bash
node plugins/dsh-bioflow/bin/bioflow.mjs status --json
node plugins/dsh-bioflow/bin/bioflow.mjs search --query="gastric cancer liver metastasis" --sources=pubmed,gdc
node plugins/dsh-bioflow/bin/bioflow.mjs health          # 五个数据源健康灯
node plugins/dsh-bioflow/bin/bioflow.mjs templates       # 六个分析模板
node plugins/dsh-bioflow/bin/bioflow.mjs advise --need="胃癌肝转移 单细胞 NK 细胞"   # AI 挑数据（打分+中文摘要+下一步）
node plugins/dsh-bioflow/bin/bioflow.mjs topics --disease="胃癌肝转移"               # 候选临床问题 + 近两年综述热点
node plugins/dsh-bioflow/bin/bioflow.mjs harmonize --dataset-ids='["ds_a","ds_b"]' --method=quantile
node plugins/dsh-bioflow/bin/bioflow.mjs suggest-next --run-id=<run_id>              # 顺着结果继续许愿
```

---

## 关键设计（为什么它不会瞎报结论）

- **闸门真拦**：样本表未知、counts 尺度未验证、每组独立 n<2、设计完全混杂、contrast 方向未定义 → 一律 `blocked`，并写明卡在哪一条。
- **不重算分析**：许愿机只做编排与留证，统计全部交给各工作台引擎，原样落盘再抽表。
- **失败不长绿**：任何一步失败，run 状态是 `failed`，`ok:false`；失败 run 仍写 manifest/design/qc/methods 供诊断，但不产图。
- **模型推演 ≠ 实验结论**：virtualcell / 签名映射 / NSCLC 组合的产物统一标 `evidence_kind: model_prediction`。
- **可回查**：每条结论能绑到 `run_id / table_id / row_id / figure_id`。

---

## 已知限制（先看这几条，省得踩坑）

1. **出图要显式给角色**：`plotlab-figure` 模板必须传 `roles`（如 `{"x":"group","y":"score"}`），不给会报「缺少角色 x 的列」。
2. **矩阵宽表不能直接出图**：矩阵会被读成「基因 + 样本列」，需要先转成 tidy（一行一个观测）再出图。
3. **数据集身份包含 data_kind**：同一份文件按 matrix 与 tidy 注册会得到两个不同数据集，这是有意的。
4. **单细胞组成分析目前只有单一分母**（全部合格细胞）；「NK / 全部合格免疫细胞」双分母尚未实现，表注已固定写明当前口径。
5. **队列×分组完全混杂现在会前置拦截**：`harmonize` 若发现每个队列只有一组且分属不同组，直接 `status=blocked`（`dataset_confounded`），不再生成可跨队列推断的合并集；历史上用旧逻辑生成的合并集需从原始队列重建。
6. **旧虚拟细胞状态需重跑**：早期结果没有保存 HVG 的 mean/std，投影量纲会降级并给 warning。
7. R 侧依赖（edgeR/limma/Seurat）缺失时功能降级而非报错，报告里一定写明实际方法。

---

## 验收自查

```bash
# 1) 插件是否都加载
node plugins/dsh-bioflow/bin/bioflow.mjs status --json     # 期望 "插件可用 6/6"

# 2) 数据源
node plugins/dsh-bioflow/bin/bioflow.mjs health            # geo/pubmed/gdc/metabolights 应为 ok

# 3) 自检
node plugins/dsh-bioflow/bin/bioflow.mjs audit --json
```

本次交付前，我们跑过：**43 个 P0 反例测试全绿**、新增 5 项评审回归（尺度/方向/混杂/串组/闸门补救）全绿、六个插件既有回归（GEO Bulk 150/31/41/94/54/14/40、SingleCell 36/10893 断言、Spatial pipeline、Metabolomics 138/89、PlotLab 204、VirtualCell 8/11/23）全绿，以及许愿机自身的 **64 项测试（63 通过 / 1 skip）**与真实端到端（两个 bulk 合并 → `ds_harm_*` → 7/7 步 succeeded、31 个产物、41 条证据）。

---

## 目录结构

```
万能生信许愿机/
├── install.sh            # 安装（登记到 DSH profile）
├── uninstall.sh          # 卸载
├── README.md             # 本文件
└── plugins/
    ├── dsh-bioflow/      # 许愿机本体（工具 + CLI）
    ├── dsh-geobulk/
    ├── dsh-singlecell/
    ├── dsh-spatial/
    ├── dsh-metabolomics/
    ├── dsh-plotlab/
    ├── dsh-virtualcell/
    └── dsh-nsclc-vt/
```
