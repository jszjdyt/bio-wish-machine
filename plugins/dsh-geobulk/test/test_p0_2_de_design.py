#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — P0-2 回归测试：差异设计与回退

反例来源（BIO_P0_BACKLOG §1 GEO Bulk 第 2 条）：
  * ``fail_on_fallback`` 解析后从未使用：显式要求「失败不回退」时仍在缺依赖时静默切 Welch；
  * ``_welch_paired`` 丢弃协变量；``_ols_fallback`` 把数值协变量 ``str()`` 成哑变量、
    设计退化时静默丢协变量（声称「协变量校正」但实际没校正）；
  * 三组取二时 ``_counts_matrix(ctx, X)`` 用切片后的 X 与全量 counts 比形状，
    DESeq2 误报「没有原始 counts」；
  * ``bg_de.R`` 的 ``round(as.matrix(expr))`` 静默圆整计数；
  * 真 counts 走 limma 时被 log2(x+1) 后当普通 limma（应为 limma-voom）；
  * 结果与方法里没有记录 actual method / contrast / 有效 n / pairs / 实际公式。

本文件是**反例矩阵**：每条断言在修复前 FAIL、修复后 PASS。为了让「修复前」也能一次
跑完并留下 before 证据，每个 section 都被包在 try/except 里，单条失败不中断。

用法（在 plugins/dsh-geobulk 下）：
    export MPLBACKEND=Agg MPLCONFIGDIR=$(mktemp -d)
    python3 test/test_p0_2_de_design.py
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="bg-mpl-"))

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bg_de  # noqa: E402
import rbridge  # noqa: E402
from bg_common import EngineError  # noqa: E402

OUTDIR = tempfile.mkdtemp(prefix="bg-p0-2-")
PASSED: list[str] = []
FAILED: list[str] = []
BEFORE: dict[str, str] = {}


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


def build_ctx(X, genes, samples, info, *, kind="lognorm", counts=None, groups=None) -> dict:
    return {
        "dataset_id": "ds_p0_2",
        "name": "synthetic",
        "genes": list(genes),
        "samples": list(samples),
        "X": np.asarray(X, dtype=np.float64),
        "X_kind": kind,
        "counts": None if counts is None else np.asarray(counts, dtype=np.float64),
        "sample_info": info,
        "group_key": "group",
        "groups": groups or [],
        "species": "human",
        "results": {},
        "modules": {},
    }


def run(ctx, params, tag):
    out = os.path.join(OUTDIR, tag)
    os.makedirs(out, exist_ok=True)
    return bg_de.run_de(ctx, params, out)


def run_catch(ctx, params, tag) -> tuple[dict | None, EngineError | None]:
    """跑一次，失败时返回 (None, EngineError)，供「修复前行为」取证。"""
    try:
        return run(ctx, params, tag), None
    except EngineError as error:
        return None, error
    except Exception as error:  # noqa: BLE001
        return None, EngineError("internal", f"{type(error).__name__}: {error}")


def expect_raise(fn) -> tuple[bool, str]:
    try:
        fn()
    except EngineError as error:
        return True, f"{error.code}: {error.message}"
    except Exception as error:  # noqa: BLE001
        return True, f"其他异常 {type(error).__name__}: {error}"
    return False, "没有报错"


def _counts_3group(n_per=4, n_genes=50, seed=17):
    """三组整数 counts：tumor 前 10 基因上调（相对 normal），met 是第三组。"""
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    base = rng.uniform(80, 500, n_genes)
    samples, labels, X = [], [], []
    for label, scale in (("tumor", 3.0), ("normal", 1.0), ("met", 1.0)):
        for j in range(n_per):
            mu = base * rng.uniform(0.85, 1.15) * scale
            mu[:10] *= 2.2 if label == "tumor" else 1.0
            X.append(rng.negative_binomial(3.0, 3.0 / (3.0 + mu)))
            samples.append(f"{label[0].upper()}{j}")
            labels.append(label)
    return np.vstack(X), genes, samples, labels


def _direct_deseq2(X, genes, samples, labels, case="tumor", control="normal") -> pd.DataFrame:
    """独立参照：直接 Rscript 跑 DESeq2（与工作台实现无关）。"""
    tmp = Path(tempfile.mkdtemp(prefix="p0-2-deseq2-"))
    pd.DataFrame(np.asarray(X).T, index=genes, columns=samples).to_csv(tmp / "counts.csv")
    pd.DataFrame({"group": labels}, index=samples).to_csv(tmp / "info.csv")
    script = tmp / "ref.R"
    script.write_text(
        """
args <- commandArgs(trailingOnly=TRUE)
d <- args[1]
counts <- as.matrix(read.csv(file.path(d, "counts.csv"), row.names=1, check.names=FALSE))
info <- read.csv(file.path(d, "info.csv"), row.names=1, check.names=FALSE)
info$group <- factor(info$group, levels=c("%s", "%s"))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=info, design=~group)
dds <- DESeq(dds, quiet=TRUE)
res <- results(dds, contrast=c("group", "%s", "%s"))
out <- data.frame(gene=rownames(res), logFC=res$log2FoldChange,
                  P.Value=res$pvalue, adj.P.Val=res$padj)
write.csv(out, file.path(d, "ref.csv"), row.names=FALSE)
""" % (control, case, case, control),
        encoding="utf-8",
    )
    proc = subprocess.run(["Rscript", "--vanilla", str(script), str(tmp)],
                          capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        raise RuntimeError(f"参考 DESeq2 调用失败：{proc.stderr[-1200:]}")
    return pd.read_csv(tmp / "ref.csv")


# ---------------------------------------------------------------------------
# A. 三组取二：矩阵 / 元数据 / counts 必须共享同一份已验证样本索引
# ---------------------------------------------------------------------------

def test_three_group_subset():
    print("\n[A] 三组取二：counts 按同一样本索引切片")
    counts, genes, samples, labels = _counts_3group()
    info = pd.DataFrame({"group": labels}, index=samples)
    # 真实场景：矩阵先标准化成 log 尺度（X_kind=lognorm），原始 counts 层保留；
    # 三组里只比较 tumor vs normal（4 vs 4），第三组 met 不参与。
    ctx = build_ctx(np.log2(counts + 1.0), genes, samples, info, kind="lognorm", counts=counts)
    params = {"method": "deseq2", "group_col": "group", "case": "tumor", "control": "normal"}
    out, error = run_catch(ctx, params, "p0-2-A")
    if error is not None:
        BEFORE["A"] = f"运行失败 {error.code}: {error.message[:90]}"
        check("A1 三组取二 + DESeq2 能跑（修复前误报缺 counts）", False, BEFORE["A"])
        for name in ("A2", "A3", "A4", "A5"):
            check(name, False, "上游未产出结果")
        return
    metrics = out["summary"]["metrics"]
    check("A1 三组取二 + DESeq2 能跑（修复前误报缺 counts）",
          metrics.get("actual_method") == "deseq2", f"actual={metrics.get('actual_method')}")
    check("A2 有效 n = 4 vs 4（第三组不参与）",
          int(metrics.get("n_case") or 0) == 4 and int(metrics.get("n_control") or 0) == 4,
          f"{metrics.get('n_case')} vs {metrics.get('n_control')}")
    keep = [i for i, lab in enumerate(labels) if lab in ("tumor", "normal")]
    ref = _direct_deseq2(counts[keep], genes, [samples[i] for i in keep],
                         [labels[i] for i in keep])
    merged = out["de_store"]["table"].set_index("gene").join(ref.set_index("gene"), rsuffix="_ref")
    lfc_diff = float(np.max(np.abs(merged["logFC"] - merged["logFC_ref"])))
    check("A3 切片后的 logFC 与直接 DESeq2 子集一致", lfc_diff < 1e-6, f"max|diff|={lfc_diff:.2e}")
    check("A4 记录 contrast / 实际公式",
          metrics.get("contrast") == "tumor-normal" and bool(metrics.get("design")),
          f"contrast={metrics.get('contrast')} design={metrics.get('design')}")
    check("A5 未发生静默回退（lite=False）", not metrics.get("lite"),
          f"lite={metrics.get('lite')}")


# ---------------------------------------------------------------------------
# B. fail_on_fallback：显式要求不回退时必须失败
# ---------------------------------------------------------------------------

def test_fail_on_fallback():
    print("\n[B] fail_on_fallback 真正生效")
    rng = np.random.default_rng(5)
    n, g = 8, 30
    genes = [f"G{i:02d}" for i in range(g)]
    samples = [f"S{i}" for i in range(n)]
    info = pd.DataFrame({"group": ["tumor"] * 4 + ["normal"] * 4}, index=samples)
    X = rng.normal(7, 0.7, (n, g))
    X[:4, :5] += 2.0
    ctx = build_ctx(X, genes, samples, info)
    original = bg_de._r_diff
    try:
        bg_de._r_diff = lambda *a, **k: {"frame": None, "failed": True,
                                         "error": "本机没有可用的 Rscript", "code": "no_r"}
        out, error = run_catch(ctx, {"method": "limma", "group_col": "group", "fail_on_fallback": True},
                               "p0-2-B")
        BEFORE["B"] = (f"运行成功（静默回退 Welch）" if out is not None
                       else f"运行失败 {error.code}: {error.message[:80]}")
        check("B1 fail_on_fallback=true + R 不可用 → 明确失败", out is None and error is not None,
              BEFORE["B"])
        msg = error.message if error is not None else ""
        check("B2 报错点名 fail_on_fallback / 不回退",
              ("fail_on_fallback" in msg) or ("回退" in msg), msg[:110])

        bg_de._r_diff = lambda *a, **k: {"frame": None, "failed": True,
                                         "error": "limma 没有安装", "code": "prerequisite"}
        out3, error3 = run_catch(ctx, {"method": "limma", "group_col": "group",
                                       "fail_on_fallback": "true"}, "p0-2-B2")
        check("B3 缺包 + fail_on_fallback（字符串 true）→ 明确失败",
              out3 is None and error3 is not None,
              f"{(error3.code + ': ' + error3.message[:80]) if error3 else '未报错'}")
        check("B3b 缺包失败保留 prerequisite 语义",
              error3 is not None and error3.code == "prerequisite",
              f"code={getattr(error3, 'code', None)}")

        bg_de._r_diff = lambda *a, **k: {"frame": None, "failed": True,
                                         "error": "本机没有可用的 Rscript", "code": "no_r"}
        out4, error4 = run_catch(ctx, {"method": "limma", "group_col": "group"}, "p0-2-B3")
        m = (out4 or {}).get("summary", {}).get("metrics", {})
        BEFORE["B4"] = (f"{m.get('method_requested')}->{m.get('actual_method')} lite={m.get('lite')}"
                        if out4 is not None else f"运行失败 {error4.message[:80]}")
        check("B4 未设 fail_on_fallback 时仍可回退并记录 requested/actual",
              out4 is not None and m.get("method_requested") == "limma"
              and m.get("actual_method") == "welch" and bool(m.get("lite")),
              BEFORE["B4"])
    finally:
        bg_de._r_diff = original


# ---------------------------------------------------------------------------
# C. 配对回退不得静默丢协变量
# ---------------------------------------------------------------------------

def _paired_fixture(n_pairs=8, n_genes=40, seed=3):
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    patient = rng.normal(0, 2.0, (n_pairs, 1))
    effect = np.zeros(n_genes)
    effect[:8] = 1.5
    X, samples, labels, pairs, purity = [], [], [], [], []
    for p in range(n_pairs):
        for arm, lab in ((0, "control"), (1, "case")):
            X.append(7.0 + patient[p] + arm * effect + rng.normal(0, 0.3, n_genes))
            samples.append(f"P{p:02d}{'T' if lab == 'case' else 'C'}")
            labels.append(lab)
            pairs.append(f"P{p:02d}")
            purity.append(float(rng.uniform(0.2, 0.9)))
    info = pd.DataFrame({"group": labels, "patient": pairs, "purity": purity}, index=samples)
    return np.vstack(X), genes, samples, info


def test_paired_covariate():
    print("\n[C] 配对回退不能静默丢协变量")
    X, genes, samples, info = _paired_fixture()
    ctx = build_ctx(X, genes, samples, info)
    original = bg_de._r_diff
    base = {"method": "limma", "group_col": "group", "paired": True, "pair_col": "patient"}
    try:
        bg_de._r_diff = lambda *a, **k: {"frame": None, "failed": True,
                                         "error": "本机没有可用的 Rscript", "code": "no_r"}
        out, error = run_catch(ctx, {**base, "covars": ["purity"]}, "p0-2-C")
        BEFORE["C"] = ("运行成功（协变量被静默丢弃）" if out is not None
                       else f"运行失败 {error.code}: {error.message[:80]}")
        check("C1 配对回退 + 协变量 → 明确失败（不静默丢弃）",
              out is None and error is not None, BEFORE["C"])
        msg = error.message if error is not None else ""
        check("C2 报错点名协变量", ("协变量" in msg) or ("purity" in msg), msg[:120])

        out2, error2 = run_catch(ctx, dict(base), "p0-2-C2")
        store = (out2 or {}).get("de_store", {})
        m = (out2 or {}).get("summary", {}).get("metrics", {})
        check("C3 无协变量时配对回退照常工作",
              out2 is not None and len(store.get("table", [])) == len(genes),
              f"rows={len(store.get('table', []))}"
              + ("" if out2 is not None else f" err={error2.message[:60]}"))
        check("C4 配对回退记录实际方法/公式/配对数/有效 n",
              str(m.get("actual_method") or "").startswith("welch")
              and "pair" in str(m.get("design") or "")
              and int(m.get("n_pairs") or 0) == 8
              and int(m.get("n_effective") or 0) == 16,
              f"actual={m.get('actual_method')} design={m.get('design')} "
              f"n_pairs={m.get('n_pairs')} n_effective={m.get('n_effective')}")
    finally:
        bg_de._r_diff = original


# ---------------------------------------------------------------------------
# D. _ols_fallback：数值协变量保持 numeric；不可识别设计必须阻止
# ---------------------------------------------------------------------------

def test_ols_numeric_covariate():
    print("\n[D] 协变量 OLS 回退（数值协变量 + 不可识别设计）")
    # case 年龄 30/40/50/60（均值 45），control 35/45/55/65（均值 50）。
    # 生成：expr = 6 + 0.5*(age-45) + 1.0*(case) + 噪声。
    # 未校正组间差 ≈ 1 - 0.5*5 = -1.5；校正 age 后应为 +1.0。
    ages = [30.0, 40.0, 50.0, 60.0, 35.0, 45.0, 55.0, 65.0]
    labels = ["case"] * 4 + ["control"] * 4
    rng = np.random.default_rng(23)
    n, g = 8, 20
    genes = [f"G{i:02d}" for i in range(g)]
    samples = [f"S{i}" for i in range(n)]
    base = 6.0 + np.outer(np.array(ages) - 45.0, np.full(g, 0.5))
    group_effect = np.zeros((n, g))
    group_effect[:4, :] = 1.0
    X = base + group_effect + rng.normal(0, 0.05, (n, g))
    info = pd.DataFrame({"group": labels, "age": ages}, index=samples)
    ctx = build_ctx(X, genes, samples, info)
    raw, raw_err = run_catch(ctx, {"method": "welch", "group_col": "group"}, "p0-2-D-raw")
    adj, adj_err = run_catch(ctx, {"method": "welch", "group_col": "group", "covars": ["age"]},
                             "p0-2-D-adj")
    raw_lfc = float(raw["de_store"]["table"].set_index("gene").loc["G00", "logFC"]) if raw else float("nan")
    adj_lfc = float(adj["de_store"]["table"].set_index("gene").loc["G00", "logFC"]) if adj else float("nan")
    BEFORE["D"] = f"未校正 logFC={raw_lfc:.3f}；加 age 后 logFC={adj_lfc:.3f}"
    check("D1 未校正的 group 差 ≈ -1.5（年龄混杂）", abs(raw_lfc - (-1.5)) < 0.12, f"{raw_lfc:.3f}")
    check("D2 校正 age 后 logFC ≈ +1.0（数值协变量按连续变量进入模型）",
          abs(adj_lfc - 1.0) < 0.12, f"{adj_lfc:.3f}")
    metrics = (adj or {}).get("summary", {}).get("metrics", {})
    design = str(metrics.get("design") or "")
    check("D3 设计公式记录 cov_age 为连续项（不是每个年龄一个哑变量）",
          "cov_age" in design and "cov_age=" not in design, f"design={design}")
    check("D4 记录有效 n", int(metrics.get("n_effective") or 0) == 8,
          f"n_effective={metrics.get('n_effective')}")

    # 完全混杂：case 全男、control 全女 → 分组与性别不可同时识别
    conf_info = info.copy()
    conf_info["sex"] = ["M"] * 4 + ["F"] * 4
    conf_ctx = build_ctx(X, genes, samples, conf_info)
    out_c, err_c = run_catch(conf_ctx, {"method": "welch", "group_col": "group", "covars": ["sex"]},
                             "p0-2-D-conf")
    BEFORE["D5"] = ("运行成功（静默丢弃共线协变量）" if out_c is not None
                    else f"运行失败 {err_c.code}: {err_c.message[:80]}")
    check("D5 完全混杂（sex=group）→ 明确失败，不静默丢协变量",
          out_c is None and err_c is not None, BEFORE["D5"])
    msg = err_c.message if err_c is not None else ""
    check("D6 报错点名不可识别的协变量", "sex" in msg or "协变量" in msg, msg[:130])

    # 请求不存在协变量 → 阻止，而不是 warning 后忽略
    out_m, err_m = run_catch(ctx, {"method": "welch", "group_col": "group",
                                   "covars": ["not_a_column"]}, "p0-2-D-missing")
    check("D7 协变量不在样本表 → 阻止并点名",
          out_m is None and err_m is not None and "not_a_column" in err_m.message,
          f"{(err_m.code + ': ' + err_m.message[:90]) if err_m else '未报错'}")

    # 恒定协变量（只有一个取值）→ 不可估计，必须点名阻止
    const_info = info.copy()
    const_info["plate"] = ["P1"] * n
    const_ctx = build_ctx(X, genes, samples, const_info)
    out_cst, err_cst = run_catch(const_ctx, {"method": "welch", "group_col": "group",
                                             "covars": ["plate"]}, "p0-2-D-const")
    check("D8 恒定协变量 → prerequisite 并点名",
          out_cst is None and err_cst is not None and "plate" in err_cst.message,
          f"{(err_cst.code + ': ' + err_cst.message[:90]) if err_cst else '未报错'}")


# ---------------------------------------------------------------------------
# E. 非法 counts：负数/分数必须阻止（Python 与 R 两侧）；整数 counts 能跑通
# ---------------------------------------------------------------------------

def _basic_counts_ctx(frac=False):
    X, genes, samples, labels = _counts_3group(n_per=4, n_genes=20)
    keep = [i for i, lab in enumerate(labels) if lab in ("tumor", "normal")]
    X = X[keep].astype(float)
    samples = [samples[i] for i in keep]
    labels = [labels[i] for i in keep]
    if frac:
        X[0, 0] = X[0, 0] + 0.5
    info = pd.DataFrame({"group": labels}, index=samples)
    return X, genes, samples, info


def _r_payload(expr, genes, samples, info):
    return {
        "expr": expr, "genes": genes, "samples": samples, "info": info,
        "params": {"method": "deseq2", "group_col": "group", "case": "tumor",
                   "control": "normal", "covars": [], "paired": False, "pair_col": ""},
    }


def test_illegal_counts():
    print("\n[E] 非法 counts（负数 / 分数）与合法整数 counts")
    X, genes, samples, info = _basic_counts_ctx()
    neg = X.copy()
    neg[1, 2] = -5.0
    neg_ctx = build_ctx(neg, genes, samples, info, kind="counts", counts=neg)
    out_n, err_n = run_catch(neg_ctx, {"method": "deseq2", "group_col": "group"}, "p0-2-E-neg")
    BEFORE["E1"] = ("运行成功（负值被 clip 成 0）" if out_n is not None
                    else f"运行失败 {err_n.code}: {err_n.message[:80]}")
    check("E1 DESeq2 遇负 counts → 阻止（不再 clip 成 0）", out_n is None and err_n is not None,
          BEFORE["E1"])
    check("E2 报错说明负值", err_n is not None and ("负" in err_n.message or "negative" in err_n.message.lower()),
          (err_n.message[:120] if err_n else "未报错"))

    Xf, genes, samples, info = _basic_counts_ctx(frac=True)
    frac_ctx = build_ctx(Xf, genes, samples, info, kind="counts", counts=Xf)
    out_f, err_f = run_catch(frac_ctx, {"method": "deseq2", "group_col": "group"}, "p0-2-E-frac")
    BEFORE["E3"] = ("运行成功（分数被四舍五入）" if out_f is not None
                    else f"运行失败 {err_f.code}: {err_f.message[:80]}")
    check("E3 DESeq2 遇分数 counts → 阻止（不再静默四舍五入）", out_f is None and err_f is not None,
          BEFORE["E3"])
    check("E4 报错说明非整数", err_f is not None and ("整数" in err_f.message or "integer" in err_f.message.lower()),
          (err_f.message[:120] if err_f else "未报错"))

    ok_ctx = build_ctx(X, genes, samples, info, kind="counts", counts=X)
    ok_out, ok_err = run_catch(ok_ctx, {"method": "deseq2", "group_col": "group"}, "p0-2-E-ok")
    check("E4b 合法整数 counts → DESeq2 成功（不是一律拒绝）",
          ok_out is not None and ok_out["summary"]["metrics"].get("actual_method") == "deseq2",
          f"{(ok_out or {}).get('summary', {}).get('metrics', {}).get('actual_method')}"
          + (f" err={ok_err.message[:60]}" if ok_err else ""))

    # R 侧兜底：绕过 Python 直接调 bg_de.R，分数/负数 counts 也不得静默 round/clip
    if rbridge.r_available() and rbridge.has_pkg("DESeq2"):
        outdir = os.path.join(OUTDIR, "p0-2-E-r")
        os.makedirs(outdir, exist_ok=True)
        out_r, err_r = None, None
        try:
            rbridge.run_r("bg_de.R", _r_payload(Xf, genes, samples, info), outdir, timeout=600)
        except EngineError as error:
            err_r = error
        BEFORE["E5"] = ("R 成功（round 后照常出结果）" if err_r is None
                        else f"R 失败 {err_r.code}: {err_r.message[:80]}")
        check("E5 bg_de.R 不再 round 分数 counts（R 侧明确失败）", err_r is not None, BEFORE["E5"])
        check("E5b R 侧分数报错说明整数", err_r is not None and "整数" in err_r.message,
              (err_r.message[:120] if err_r else "未报错"))

        out_r2, err_r2 = None, None
        try:
            rbridge.run_r("bg_de.R", _r_payload(neg, genes, samples, info), outdir, timeout=600)
        except EngineError as error:
            err_r2 = error
        check("E6 bg_de.R 对负 counts 明确失败",
              err_r2 is not None and "负" in err_r2.message,
              f"{(err_r2.code + ': ' + err_r2.message[:90]) if err_r2 else 'R 未报错'}")
    else:
        print("  SKIP E5/E6 本机没有 R + DESeq2")


# ---------------------------------------------------------------------------
# F. 矩阵尺度 → 方法路由（真 counts 走 limma-voom）
# ---------------------------------------------------------------------------

def _counts_2group_clean(n_per=5, n_genes=40, seed=29, effect_genes=8, effect=4.0):
    """两组整数 counts：只有第一组前 effect_genes 个基因真实上调，文库大小相当。

    F 段只验证「尺度 → 方法路由」，所以不制造整体文库大小偏移，让方向判据稳定。
    """
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    base = rng.uniform(120, 600, n_genes)
    samples, labels, X = [], [], []
    for label in ("tumor", "normal"):
        for j in range(n_per):
            mu = base * rng.uniform(0.9, 1.1, n_genes)
            if label == "tumor":
                mu[:effect_genes] *= effect
            X.append(rng.negative_binomial(5.0, 5.0 / (5.0 + mu)))
            samples.append(f"{label[0].upper()}{j}")
            labels.append(label)
    return np.vstack(X), genes, samples, labels


def test_scale_routing():
    print("\n[F] 矩阵尺度 → 方法路由")
    if not (rbridge.r_available() and rbridge.has_pkg("limma")):
        print("  SKIP 本机没有 R + limma")
        for name in ("F1", "F2", "F3", "F4", "F5"):
            check(name, False, "本机没有 R + limma")
        return
    Xk, genes, samples, labels = _counts_2group_clean()
    info = pd.DataFrame({"group": labels}, index=samples)
    counts_ctx = build_ctx(Xk, genes, samples, info, kind="counts", counts=Xk)
    out, error = run_catch(counts_ctx, {"method": "limma", "group_col": "group"}, "p0-2-F-counts")
    m = (out or {}).get("summary", {}).get("metrics", {})
    BEFORE["F"] = f"actual_method={m.get('actual_method')} value_kind={m.get('value_kind')}"
    check("F1 真 counts + limma → limma-voom（不是 log2(x+1) 普通 limma）",
          m.get("actual_method") == "limma-voom",
          f"actual={m.get('actual_method')}" + (f" err={error.message[:60]}" if error else ""))
    table = (out or {}).get("de_store", {}).get("table")
    clean = table.set_index("gene") if table is not None and len(table) else pd.DataFrame()
    up_genes = [f"G{i:03d}" for i in range(8)]
    if not clean.empty:
        up_lfc = clean.loc[[g for g in up_genes if g in clean.index], "logFC"].astype(float)
        mean_up = float(up_lfc.mean())
        n_pos = int((up_lfc > 0).sum())
    else:
        mean_up, n_pos = float("nan"), 0
    check("F2 voom 方向正确（真上调基因集平均 logFC>0.5 且 ≥7/8 为正）",
          mean_up > 0.5 and n_pos >= 7, f"mean logFC={mean_up:.3f}, 正号 {n_pos}/8")

    log_ctx = build_ctx(np.log2(Xk + 1.0), genes, samples, info, kind="lognorm", counts=Xk)
    out2, err2 = run_catch(log_ctx, {"method": "limma", "group_col": "group"}, "p0-2-F-log")
    m2 = (out2 or {}).get("summary", {}).get("metrics", {})
    check("F3 log 矩阵 + limma → limma（不误用 voom）", m2.get("actual_method") == "limma",
          f"actual={m2.get('actual_method')}" + (f" err={err2.message[:60]}" if err2 else ""))
    check("F4 log 矩阵结果记录输入尺度", m2.get("value_kind") == "lognorm",
          f"value_kind={m2.get('value_kind')}")

    plain_ctx = build_ctx(np.log2(Xk + 1.0), genes, samples, info, kind="lognorm")
    out5, err5 = run_catch(plain_ctx, {"method": "deseq2", "group_col": "group"}, "p0-2-F-nocounts")
    check("F5 log 矩阵 + DESeq2 → 阻止（要求原始 counts）",
          out5 is None and err5 is not None,
          f"{(err5.code + ': ' + err5.message[:80]) if err5 else '未报错'}")


# ---------------------------------------------------------------------------
# G. 记录 actual method / contrast / 有效 n / pairs / 实际公式
# ---------------------------------------------------------------------------

def test_record_provenance():
    print("\n[G] 结果与方法记录完整性")
    counts, genes, samples, labels = _counts_3group(n_per=4, n_genes=24)
    keep = [i for i, lab in enumerate(labels) if lab in ("tumor", "normal")]
    counts = counts[keep]
    samples = [samples[i] for i in keep]
    labels = [labels[i] for i in keep]
    info = pd.DataFrame({"group": labels}, index=samples)
    ctx = build_ctx(counts, genes, samples, info, kind="counts", counts=counts)
    if not (rbridge.r_available() and rbridge.has_pkg("DESeq2")):
        print("  SKIP G 本机没有 R + DESeq2")
        check("G0 记录完整性（需要 R + DESeq2）", False, "SKIPPED")
        return
    out, error = run_catch(ctx, {"method": "deseq2", "group_col": "group"}, "p0-2-G")
    if out is None:
        BEFORE["G"] = f"运行失败 {error.code}: {error.message[:80]}"
        for name in ("G1", "G2", "G3"):
            check(name, False, BEFORE["G"])
        return
    metrics = out["summary"]["metrics"]
    store = out["de_store"]
    BEFORE["G"] = (f"actual={metrics.get('actual_method')} contrast={metrics.get('contrast')} "
                   f"design={metrics.get('design')} n_eff={metrics.get('n_effective')}")
    check("G1 metrics 记录 actual method / contrast / 有效 n / 公式",
          metrics.get("actual_method") == "deseq2"
          and metrics.get("contrast") == "tumor-normal"
          and int(metrics.get("n_effective") or 0) == 8
          and bool(metrics.get("design")),
          BEFORE["G"])
    check("G2 de_store 保存 method/design/contrast/n_effective/n_pairs/value_kind",
          store.get("method") == "deseq2"
          and bool(store.get("design"))
          and store.get("contrast") == "tumor-normal"
          and int(store.get("n_effective") or 0) == 8
          and int(store.get("n_pairs") or 0) == 0
          and store.get("value_kind") == "counts",
          f"store keys: method={store.get('method')} design={store.get('design')} "
          f"n_eff={store.get('n_effective')} n_pairs={store.get('n_pairs')} "
          f"value_kind={store.get('value_kind')}")
    notes = " ".join(str(item) for item in (out.get("notes") or []))
    bullets = " ".join(str(item) for item in ((out.get("summary") or {}).get("bullets") or []))
    check("G3 结果文字写明实际方法与对比方向",
          ("DESeq2" in bullets or "DESeq2" in notes) and "tumor" in bullets,
          f"bullets={bullets[:120]}")


# ---------------------------------------------------------------------------
# 入口
# ---------------------------------------------------------------------------

SECTIONS = [
    ("A 三组取二", test_three_group_subset),
    ("B fail_on_fallback", test_fail_on_fallback),
    ("C 配对+协变量", test_paired_covariate),
    ("D 数值协变量/不可识别设计", test_ols_numeric_covariate),
    ("E 非法 counts", test_illegal_counts),
    ("F 尺度路由 voom", test_scale_routing),
    ("G 记录完整性", test_record_provenance),
]


def main() -> int:
    print("== dsh-geobulk P0-2 差异设计与回退 ==")
    print(f"Python {sys.version.split()[0]} | R {rbridge.probe().get('version') or '不可用'}")
    for name, fn in SECTIONS:
        try:
            fn()
        except Exception as error:  # noqa: BLE001
            FAILED.append(f"{name} 崩溃")
            print(f"  FAIL {name} 崩溃：{type(error).__name__}: {error}")
            traceback.print_exc(limit=3)
    print("\n== 修复前行为快照 ==")
    for key, value in BEFORE.items():
        print(f"  [{key}] {value}")
    print()
    print(f"== 结果：{len(PASSED)} 通过 / {len(FAILED)} 失败 ==")
    if FAILED:
        for item in FAILED:
            print(f"  - {item}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
