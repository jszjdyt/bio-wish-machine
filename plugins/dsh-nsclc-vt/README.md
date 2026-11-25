# dsh-nsclc-vt — NSCLC「虚拟肿瘤」插件

把 [jfisher-lab/NSCLC](https://github.com/jfisher-lab/NSCLC) 那套**虚拟肿瘤**算法搬进 DeepSeek Harness：
用纯 JavaScript 重写了 BioModelAnalyzer 的 **VMCAI 稳定性引擎**，在本地进程里直接跑模型，
不需要 R、不需要 Windows、不需要 `BioCheckConsole.exe`、不联网。

> Clarke, Barker & Nicholls et al., *Virtual Tumors Enable Prediction of Personalized
> Therapeutic Combinations for Non-Small Cell Lung Cancer*, **Cancer Research 2026**,
> doi:[10.1158/0008-5472.CAN-25-2476](https://doi.org/10.1158/0008-5472.CAN-25-2476)

---

## 一句话：它能干什么

给它一个基因型（比如「RAS 突变 + p53 缺失」）和几个干预（敲低某个基因 / 用药），
它告诉你这个虚拟肿瘤会**稳定在什么状态**、增殖和死亡各是多少，以及**哪些节点变化最大**；
再往下，它能在 7 个已发表的细胞系背景上做**单药 / 双药全筛**，按论文自己的
additivity 打分给组合排序。

---

## 四个工具

| 工具 | 用途 |
|---|---|
| `nsclc_model` | 模型概览：162 个节点、360 条关系、7 个已发表背景；给 `query` 可按名字/基因符号检索节点 |
| `nsclc_node` | 单个节点的「机制卡」：target function、上下游、策展的通路证据与文献出处；给定基因型还能给出它的稳态水平 |
| `nsclc_predict` | 跑一次模拟：输出表型（Proliferation / Death）、`Difference`、变化最大的节点、收敛步数 |
| `nsclc_screen` | 单药 / 双药筛选，跨背景排序；支持论文的 additivity 打分与「至少在 N 个背景里优于最优单药」过滤；可导出 CSV |

### 命名约定

- **敲除 = 把节点钉到该节点的最低水平**（这个模型里就是 0）。`interventions: ["ATM"]` 表示敲除，
  `["MYC:3"]` 表示强制激活到 3。范围是整数等级 `0..3`，不是浓度。
- `Difference = mean(Proliferation) - mean(Death)`，**越低越好**。
- `additive = Difference(双药) - min(两个单药的 Difference)`，**负 = 协同**（比最优单药还好），
  正 = 拮抗。这就是论文 Figure 5 的定义。
- `uncertainty = (hi - lo) / 范围跨度`。等于 0 才代表引擎把该节点钉死在唯一取值上；
  非 0 说明该节点在一个区间里，网络没有证明唯一稳态。

---

## 安装

插件位于本仓库的 `plugins/dsh-nsclc-vt/`。挂载到某个 DSH profile：

```bash
# 1) 在 profile 里登记为依赖（link 指向本目录）与 bundle
cd ~/.dsh/profiles/desktop
# 编辑 package.json：
#   "dependencies": { "dsh-nsclc-vt": "link:/绝对路径/plugins/dsh-nsclc-vt" }
#   "dsh": { "profile": { "bundles": [ ..., "dsh-nsclc-vt" ] } }

# 2) 安装依赖
pnpm install

# 3) 重启 DeepSeek Harness（长驻进程在启动时解析 bundle）
```

装好后对话里直接说需求即可，例如：

- 「A549 背景下敲掉 MEK 和 PI3K，虚拟肿瘤怎么说？」
- 「在 pc9 和 hcc827 两个背景里筛一下 PI3K/AKT 通路的双靶点组合，按协同排序。」
- 「ATM 这个节点在模型里是怎么被调控的？有文献支持吗？」

也可以让它导出完整结果表做图：`nsclc_screen` 的 `output_csv` 参数会写出逐背景的 CSV。

---

## 算法到底是什么

这不是「神经网络预测」，是一个**可解释的定性网络 + 抽象解释**：

1. **模型**是一组带整数等级（0–3）的变量和它们的目标函数（target function），比如
   `ATM = MRN * max(0, min(1, 1 - ATM_deficient)) + min(1, ATM_deficient * MRN)`。
   公式里没写 target function 的节点，导入时按 BioModelAnalyzer 的规则从它的
   activator / inhibitor 边**合成**一个（`avg(pos) - avg(neg)` 那一套）。
2. **更新是同步的**：一个 tick 里所有变量同时按各自的目标函数更新。
3. **VMCAI 引擎**（Cook, Fisher, Krepska & Piterman, *Proving Stabilization of
   Biological Systems*, VMCAI 2011）做的是**区间混沌迭代**：每个变量从自己的完整取值范围出发，
   反复用「输入区间 → 目标函数取值上下界」收紧自己，被收紧的变量把它的下游重新推入 frontier，
   直到 frontier 排空。排空时如果所有区间都塌缩成一个点，就证明了**唯一稳态**——
   `Ticks[0]` 就是那个稳态；否则返回的是**稳态集合的区间包络**，不是某条轨迹。
4. `(Lo, Hi)` 界定的是**稳态（不动点）集合**，不是可达状态集，也不是时间轨迹。
   `Time` 字段是**收敛步数**（frontier 弹出次数），不是生物学时间。

### 实现是逐行对齐源码的

移植自 `hallba/BioModelAnalyzer` 的 `Stabilize.fs` / `GenLemmas.fs` / `FNewLemmas.fs` /
`Expr.fs` / `Marshal.fs`，包括这些容易踩的细节：

- `var(id)` 不是直接取值，而是**按读取方的范围线性重标定**，再四舍五入（`x + 0.5` 后截断）并夹取到节点自身范围；
- 公式里出现、但不是该节点关系边的输入，在求值时按节点下界处理；
- frontier 按**变量 id 从大到小**弹栈，所以 `Time` 可以逐位复现；
- `-ko id level` 是把目标函数**替换成常量**（永久钉住），不是初值约束；同一节点多次 `-ko` 时**命令行里靠前的那个生效**，
  所以 NANSEN 把干预放在背景之前——本插件同样「干预优先于背景」。

### 验证

| 验证 | 结果 |
|---|---|
| 论文 `results/combination_screen/node_results.csv` 全量逐例比对（123,711 个 case、247,422 个表型值，含 `Time`） | **0 处不一致** |
| 论文 `processed_combination/additive_scores.csv`（Figure 5 的协同矩阵，4,371 组合 × 7 细胞系 = 30,597 格） | **30,597 / 30,597 完全一致** |
| 单元测试 | `node --test test/*.test.mjs` → 23 项全通过 |

自己复核：

```bash
curl -L -o /tmp/additive_scores.csv \
  https://raw.githubusercontent.com/jfisher-lab/NSCLC/main/results/combination_screen/processed_combination/additive_scores.csv
node test/validate-published.mjs /tmp/additive_scores.csv
# → compared 30597 published cells: 30597 exact, 0 mismatched
```

---

## 数据来源与引用

`data/` 里的文件逐字复制自 [jfisher-lab/NSCLC](https://github.com/jfisher-lab/NSCLC) `model/`：

- `NSCLC-master.json` — 模型本体（162 节点 / 360 关系）
- `nodes_info.csv`、`edges_master.csv`、`target_functions_master.csv` — 节点与边的策展注释、文献出处
- `combo_coker_mut_background.csv` — 论文筛选用的 7 个细胞系背景
- `spec_master.csv` — 916 行模型-实验一致性规格（随插件分发，供查询/延伸使用）

模型与数据版权归原作者。本插件只做算法移植，若用于发表请引用上面的论文。

---

## 配置

```yaml
- id: nsclc-vt
  config:
    modelPath: /path/to/another-bma-model.json   # 默认用内置的 NSCLC-master.json
    phenotypes: [Proliferation, Death]           # 表型读数节点
    maxAnalyses: 60000                           # 超过这个模拟次数的筛选直接拒绝执行
    workers: 0                                   # 0 = 自动 min(8, 核数-1)
    promptGuidance: true                         # 是否注入系统提示说明
    promptOrder: 120
```

引擎与模型解耦：把 `modelPath` 指向任何 BMA 定性的 `.json`，工具就换一套网络跑。

---

## 开发

```bash
node --test test/*.test.mjs      # 单元 + 黄金样例 + 工具契约测试
node test/validate-published.mjs <additive_scores.csv>   # 全量复现论文矩阵
```

模块划分：

```
lib/expr.js       公式词法/语法、求值、单调性分析（Expr.fs 移植）
lib/model.js      BMA JSON 导入、默认目标函数合成、-ko 变换（Marshal.fs / QN.fs 移植）
lib/vmcai.js      stabilize_lazy + tighten_fast（GenLemmas.fs / FNewLemmas.fs 移植）
lib/analysis.js   基因型解析、单次预测
lib/screen.js     筛选枚举、worker 并行、additivity 打分与排序
lib/tools/        四个 DSH 工具的 schema + 渲染 + 执行
lib/host.js       宿主服务解析（带本地 fallback，脱离 DSH 也能跑测试）
data/             模型与策展注释
test/             测试与黄金样例
```

---

## 边界

- 结果是**定性模型上的预测**，用来排优先级、设计实验，**不是临床建议**。
- 等级 0–3 是抽象水平，不可解读为剂量、IC50 或表达量。
- 论文原文的 19 基因放疗响应 signature 属于另一条分析链路，本插件没有包含。
- 模型只有 162 个节点、360 条边，覆盖的是肿瘤内在信号通路；微环境、免疫、药代动力学都不在里面。

## English summary

`dsh-nsclc-vt` brings the NSCLC virtual-tumor model (Clarke et al., Cancer Research 2026)
and the BioModelAnalyzer VMCAI stability engine into DeepSeek Harness as four tools —
`nsclc_model`, `nsclc_node`, `nsclc_predict`, `nsclc_screen` — implemented as a line-by-line
JavaScript port of `Stabilize.fs`, `GenLemmas.fs`, `FNewLemmas.fs`, `Expr.fs` and `Marshal.fs`,
with no R, .NET, Windows or network dependency. It reproduces all 30,597 cells of the
manuscript's published combination-benefit matrix and all 247,422 phenotype values of the
published screen exactly, including BMA's refinement-step counter.

MIT licensed. Model and data © the original authors.
