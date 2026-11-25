#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P0-2 回归：供体级差异表达不能再拿「细胞」当生物学重复。

反例（修复前实测）：
  1. unit=donor 时实现把每个供体的 log 表达取均值后跑独立 Welch t ——
     既不是 counts 层，也没有库大小 offset，method 却写成「pseudobulk」；
  2. 供体列里的 NaN 被 ``str(value)`` 变成字符串 ``"nan"``，凭空多出一个供体；
  3. ``donor_min_cells`` 没进 DE 缓存键 → 改这个参数仍命中旧结果；
  4. 只含单一细胞类型的对象被 ``_group_labels`` 以「只有一个取值」拒绝。

本测试直接驱动 ``sc_analyses_sets.run("de_group", ...)``：
  * [1] pseudobulk 必须对**已验证原始 counts**按 (供体 × 细胞类型) 求和；
  * [2] 主路径是 edgeR quasi-likelihood（库大小 offset），并记录
        method_requested / method_actual / fallback / fallback_reason；
  * [3] n 是供体数，不是细胞数；FDR 校正后真差异基因可检出；
  * [4] donor_min_cells 进缓存键（不同取值不得复用同一结果）；
  * [5] 缺失 participant 不得变成字符串供体，必须排除并记录；
  * [6] 单一细胞类型对象不得被拒；
  * [7] counts 未验证时只能降级为「均值比较」，不得冒充 counts pseudobulk；
  * [8] 配对设计（同一供体跨两组）必须走配对模型并记录；
  * [9] 全量整数校验能识别「整数前缀 + 小数尾部」的假 counts。

运行：cd plugins/dsh-singlecell && python3 test/test_p0_donor_stats.py
退出码 0 = 全部 PASS。
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
LIB = HERE.parent / "lib"
TMP_ROOT = Path(tempfile.mkdtemp(prefix="sc-p0-donor-"))
OUT = TMP_ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(TMP_ROOT / "mpl"))
sys.dont_write_bytecode = True

PASSED = []
FAILED = []


def check(name, cond, detail=""):
    if cond:
        PASSED.append(name)
        print("  PASS %s%s" % (name, ("  [" + str(detail)[:200] + "]") if detail else ""))
        return True
    FAILED.append(name)
    print("  FAIL %s%s" % (name, ("  [" + str(detail)[:200] + "]") if detail else ""))
    return False


def _load(module_name, filename):
    spec = importlib.util.spec_from_file_location(module_name, str(LIB / filename))
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


MODULE = _load("_p0_donor_analyses", "sc_analyses_sets.py")


# ---------------------------------------------------------------------------
# 合成上下文：每个供体 × 细胞类型 5 个细胞；group B 的 G0000–G0019 上调
# ---------------------------------------------------------------------------


def _log_norm(counts):
    counts = sp.csr_matrix(np.asarray(counts, dtype=np.float64))
    row_sum = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.zeros_like(row_sum)
    np.divide(1e4, row_sum, out=scale, where=row_sum > 0)
    matrix = (sp.diags(scale) @ counts).tocsr()
    matrix.data = np.log1p(matrix.data)
    matrix.eliminate_zeros()
    return matrix.astype(np.float32)


def build_ctx(n_genes=240, per_donor_cells=5, seed=3, cell_types=("T cell", "B cell"),
              donors_per_group=4, paired=False, missing_donor_for=0):
    rng = np.random.default_rng(seed)
    genes = ["G%04d" % i for i in range(n_genes)]
    obs_rows = []
    count_rows = []
    if paired:
        donor_groups = {"d%d" % (i + 1): [0, 1] for i in range(donors_per_group)}
        names = sorted(donor_groups)
    else:
        donor_groups = {}
        names = []
        for i in range(donors_per_group):
            name = "dA%d" % (i + 1)
            donor_groups[name] = [0]
            names.append(name)
        for i in range(donors_per_group):
            name = "dB%d" % (i + 1)
            donor_groups[name] = [1]
            names.append(name)
    group_labels = ["ctrl", "case"]
    for donor in names:
        for group_index in donor_groups[donor]:
            for cell_type in cell_types:
                baseline = rng.poisson(4.0, size=n_genes).astype(np.float64)
                block = np.zeros((per_donor_cells, n_genes), dtype=np.float64)
                for row in range(per_donor_cells):
                    block[row] = rng.poisson(baseline + 1.0).astype(np.float64)
                if group_index == 1:
                    block[:, 0:20] += rng.poisson(30.0, size=(per_donor_cells, 20))
                for row in range(per_donor_cells):
                    counts = block[row]
                    if counts.sum() <= 0:
                        counts[0] = 1.0
                    count_rows.append(counts)
                    obs_rows.append({
                        "donor": donor,
                        "group": group_labels[group_index],
                        "cell_type": cell_type,
                    })
    counts = np.vstack(count_rows)
    obs = pd.DataFrame(obs_rows, index=["cell_%04d" % i for i in range(counts.shape[0])])
    if missing_donor_for:
        # 把 case 组某几个供体的 participant 置空：绝不能变成字符串 "nan" 供体
        bad = ["dB%d" % (i + 1) for i in range(missing_donor_for)]
        mask = obs["donor"].isin(bad).to_numpy()
        obs.loc[mask, "donor"] = np.nan
        obs.loc[mask & (obs["group"] == "ctrl").to_numpy(), "donor"] = ""
    matrix = sp.csr_matrix(counts.astype(np.float32))
    return {
        "dataset_id": "ds-donor",
        "dataset_version": 0,
        "dataset_revision": 0,
        "content_revision": "donor-fixture",
        "name": "donor-fixture",
        "n_obs": int(counts.shape[0]),
        "n_vars": int(n_genes),
        "X": _log_norm(counts),
        "X_kind": "lognorm",
        "counts": matrix,
        "counts_verified": True,
        "counts_source": "fixture:layers.counts",
        "obs": obs,
        "var_names": genes,
        "obsm": {},
        "clusters": None,
        "group_keys": ["group", "cell_type", "donor"],
        "engine_mode": "builtin",
        "rng": np.random.default_rng(0),
    }


def run_de(ctx, params, clear_cache=True):
    if clear_cache and hasattr(MODULE, "_DE_CACHE"):
        MODULE._DE_CACHE.clear()
    return MODULE.run("de_group", ctx, params, str(OUT / ("de-%d" % len(os.listdir(OUT)))))


# group1=case / group2=ctrl：logFC>0 表示 case 组上调（与差异方向约定一致）
DONOR_PARAMS = {"group": "group", "cell_type": "cell_type", "unit": "donor",
                "donor": "donor", "min_cells": 1, "donor_min_cells": 1, "top_n": 20,
                "group1": "case", "group2": "ctrl"}


# ---------------------------------------------------------------------------
# [1] pseudobulk = 按 (供体 × 细胞类型) 求和原始 counts
# ---------------------------------------------------------------------------


def test_pseudobulk_sums_counts():
    print("\n[1] pseudobulk 必须是原始 counts 求和（不是 log 表达均值）")
    ctx = build_ctx()
    fn = getattr(MODULE, "_pseudobulk_counts", None)
    if not callable(fn):
        check("1a 存在 _pseudobulk_counts(counts, rows, donors, levels)", False, "函数缺失")
        return
    counts = sp.csr_matrix(ctx["counts"])
    donors = np.asarray(ctx["obs"]["donor"].astype(str).tolist(), dtype=object)
    levels = sorted(set(donors.tolist()))
    rows = np.flatnonzero(np.asarray(ctx["obs"]["cell_type"] == "T cell"))
    sums, sizes, used = fn(counts, rows, donors, levels)
    manual = np.zeros((len(levels), counts.shape[1]), dtype=np.float64)
    for i, donor in enumerate(levels):
        index = rows[np.asarray([donors[row] == donor for row in rows])]
        if index.size:
            manual[i] = np.asarray(counts[index].sum(axis=0)).ravel()
    check("1a pseudobulk 行=供体数、列=基因数", np.asarray(sums).shape == manual.shape,
          "%s vs %s" % (np.asarray(sums).shape, manual.shape))
    check("1b pseudobulk 逐值等于该供体该细胞类型 counts 之和",
          np.allclose(np.asarray(sums), manual), "max diff=%.3g" % float(np.max(np.abs(np.asarray(sums) - manual))))
    expected_cells = int(np.asarray(ctx["obs"]["donor"] == levels[0]).sum() * 0 )
    first_donor_cells = int(sum(1 for row in rows if donors[row] == levels[0]))
    check("1c 每个供体的细胞数单独记录（不是把细胞数当 n）",
          [int(item) for item in sizes][:1] == [first_donor_cells],
          "%s vs %s" % (list(sizes)[:2], first_donor_cells))
    check("1d 求和结果是计数尺度（最大值远大于 log 值）", float(np.max(manual)) > 5.0,
          "max=%s" % float(np.max(manual)))


def test_full_integer_scan():
    print("\n[2] 全量整数校验（不能只抽样前 2000 个非零值）")
    fn = getattr(MODULE, "_all_integer_counts", None)
    if not callable(fn):
        check("2a 存在 _all_integer_counts", False, "函数缺失")
        return
    n_genes = 2600
    data = np.ones((2, n_genes), dtype=np.float64)
    data[:, 2000:] = 0.375
    ok, reason = fn(sp.csr_matrix(data))
    check("2a 整数前缀 + 小数尾部 → 判为非法 counts", ok is False, str(reason))
    ok2, _ = fn(sp.csr_matrix(np.array([[3.0, 0.0], [1.0, 2.0]])))
    check("2b 全整数矩阵 → 判为合法 counts", ok2 is True)
    negative = np.array([[3.0, -1.0]], dtype=np.float64)
    ok3, _ = fn(sp.csr_matrix(negative))
    check("2c 负数 counts → 判为非法", ok3 is False)
    infinite = np.array([[3.0, np.inf]], dtype=np.float64)
    ok4, _ = fn(sp.csr_matrix(infinite))
    check("2d 非有限值 counts → 判为非法", ok4 is False)


# ---------------------------------------------------------------------------
# [3] 主路径：edgeR QL / 供体级 n / 真差异可检出
# ---------------------------------------------------------------------------


def test_donor_ql_primary():
    print("\n[3] unit=donor 走 counts pseudobulk + edgeR QL")
    ctx = build_ctx()
    result = run_de(ctx, dict(DONOR_PARAMS))
    metrics = result["summary"]["metrics"]
    check("3a statistical_unit=donor", metrics.get("statistical_unit") == "donor",
          str(metrics.get("statistical_unit")))
    requested = str(metrics.get("method_requested") or "")
    actual = str(metrics.get("method_actual") or "")
    available = bool(getattr(MODULE, "_edger_available", lambda: False)())
    if available:
        check("3b 请求与实走方法都是 edgeR QL（库大小 offset）",
              requested == "edger_ql" and actual == "edger_ql",
              "requested=%s actual=%s" % (requested, actual))
    else:
        check("3b 无 edgeR 时明确降级并记录原因",
              actual in ("donor_logcpm_welch", "donor_mean_welch")
              and bool(metrics.get("fallback")) and bool(metrics.get("fallback_reason")),
              "requested=%s actual=%s" % (requested, actual))
    check("3c 方法标签不再写「Welch t on 细胞」",
          "Welch t + BH" not in str(result.get("method") or ""),
          str(result.get("method"))[:80])
    check("3d n_in/n_out 是供体数（4/4），不是细胞数",
          int(metrics.get("n_donors_in") or 0) == 4 and int(metrics.get("n_donors_out") or 0) == 4,
          "%s/%s" % (metrics.get("n_donors_in"), metrics.get("n_donors_out")))
    summary_table = next((t for t in result["tables"] if t.get("id") == "de_summary"), None)
    if summary_table is None:
        check("3e de_summary 存在", False, str([t.get("id") for t in result["tables"]]))
        return
    rows = summary_table["rows"]
    n_in = [int(row.get("n_in") or 0) for row in rows]
    check("3e 表里的 n 也是供体数（不是细胞数 20/细胞）",
          n_in and max(n_in) == 4, str(n_in[:3]))
    de_table = next((t for t in result["tables"] if t.get("id") == "de_genes"), None)
    if de_table is None:
        check("3f de_genes 存在", False, "")
        return
    hits = [row for row in de_table["rows"]
            if str(row.get("gene")) in {"G%04d" % i for i in range(20)}
            and float(row.get("logfc") or 0) > 0]
    sig = [row for row in hits if row.get("pval_adj") is not None and float(row["pval_adj"]) < 0.05]
    check("3f 20 个真上调基因里至少 10 个 FDR<0.05 且 logFC>0",
          len(sig) >= 10, "%d/%d 显著" % (len(sig), len(hits)))
    check("3g 每个细胞类型都产出结果（2 个细胞类型）",
          int(metrics.get("n_cell_types") or 0) == 2, str(metrics.get("n_cell_types")))


# ---------------------------------------------------------------------------
# [4] donor_min_cells 进缓存键
# ---------------------------------------------------------------------------


def test_donor_min_cells_cache_key():
    print("\n[4] donor_min_cells 必须进 DE 缓存键")
    fn = getattr(MODULE, "_de_cache_key", None)
    if not callable(fn):
        check("4a 存在 _de_cache_key", False, "函数缺失")
        return
    ctx = build_ctx()
    base = {"group": "group", "cell_type": "cell_type", "unit": "donor", "donor": "donor"}
    key_low = fn(ctx, dict(base, donor_min_cells=1))
    key_high = fn(ctx, dict(base, donor_min_cells=50))
    check("4a donor_min_cells 不同 → 缓存键不同", key_low != key_high,
          "%s... vs %s..." % (str(key_low)[:24], str(key_high)[:24]))
    if hasattr(MODULE, "_DE_CACHE"):
        MODULE._DE_CACHE.clear()
    ctx2 = build_ctx()
    run_de(ctx2, dict(DONOR_PARAMS), clear_cache=False)
    error_2 = None
    result_2 = None
    try:
        result_2 = run_de(ctx2, dict(DONOR_PARAMS, donor_min_cells=50), clear_cache=False)
    except Exception as error:  # noqa: BLE001
        error_2 = error
    warnings_2 = list((result_2 or {}).get("warnings") or [])
    check("4b 改 donor_min_cells 后不复用旧缓存（重算或明确拒绝）",
          error_2 is not None or not any("复用" in str(item) for item in warnings_2),
          str([w for w in warnings_2 if "复用" in str(w)])[:120])
    check("4c donor_min_cells=50 时没有供体达标 → 明确 prerequisite 而不是给假结果",
          error_2 is not None and str(getattr(error_2, "code", "")) == "prerequisite"
          and ("供体" in str(getattr(error_2, "message", "")) or "细胞" in str(getattr(error_2, "message", ""))),
          "%s: %s" % (type(error_2).__name__, getattr(error_2, "message", error_2)))


# ---------------------------------------------------------------------------
# [5] 缺失 participant 不得变成字符串供体
# ---------------------------------------------------------------------------


def test_missing_participant():
    print("\n[5] 缺失 participant 不得变成字符串供体")
    ctx = build_ctx(missing_donor_for=1)
    result = run_de(ctx, dict(DONOR_PARAMS))
    metrics = result["summary"]["metrics"]
    donors_reported = []
    for key in ("donor_levels_in", "donor_levels_out"):
        donors_reported.extend([str(item) for item in (metrics.get(key) or [])])
    lowered = [item.lower() for item in donors_reported]
    check("5a 报告里没有 'nan' / 'none' / 空字符串供体",
          not any(item in ("nan", "none", "", "na", "<na>") for item in lowered),
          str(donors_reported[:8]))
    total_donors = int(metrics.get("n_donors_in") or 0) + int(metrics.get("n_donors_out") or 0)
    check("5b 被置空的供体已从供体级分析排除（两组有效供体合计 7 个）",
          total_donors == 7,
          "case=%s ctrl=%s" % (metrics.get("n_donors_in"), metrics.get("n_donors_out")))
    warnings = " ".join(str(item) for item in (result.get("warnings") or []))
    check("5c 明确记录缺失 participant 被排除",
          any(token in warnings for token in ("缺失", "排除", "没有供体", "participant")),
          warnings[:160])


# ---------------------------------------------------------------------------
# [6] 单一细胞类型不得被拒
# ---------------------------------------------------------------------------


def test_single_cell_type():
    print("\n[6] 只含单一细胞类型的对象不得被拒")
    ctx = build_ctx(cell_types=("T cell",))
    try:
        result = run_de(ctx, dict(DONOR_PARAMS))
        metrics = result["summary"]["metrics"]
        check("6a 单细胞类型对象可以跑 de_group（不再 prerequisite 拒绝）",
              int(metrics.get("n_cell_types") or 0) == 1, str(metrics.get("n_cell_types")))
        tables = [t.get("id") for t in result["tables"]]
        check("6b 单细胞类型仍产出 de_genes / de_summary 表",
              "de_genes" in tables and "de_summary" in tables, str(tables))
    except Exception as error:  # noqa: BLE001
        check("6a 单细胞类型对象可以跑 de_group（不再 prerequisite 拒绝）", False,
              "%s: %s" % (type(error).__name__, getattr(error, "message", error)))


# ---------------------------------------------------------------------------
# [7] counts 未验证 → 只能叫「均值比较」
# ---------------------------------------------------------------------------


def test_unverified_counts_fallback():
    print("\n[7] counts 未验证时降级为「均值比较」，不得冒充 counts pseudobulk")
    ctx = build_ctx()
    ctx["counts"] = sp.csr_matrix(np.asarray(ctx["X"].todense(), dtype=np.float32))
    ctx["counts_verified"] = False
    result = run_de(ctx, dict(DONOR_PARAMS))
    metrics = result["summary"]["metrics"]
    check("7a method_requested 记录想要的 QL",
          str(metrics.get("method_requested") or "") == "edger_ql",
          str(metrics.get("method_requested")))
    check("7b 实走方法明确标为均值比较（donor_mean_welch）",
          str(metrics.get("method_actual") or "") == "donor_mean_welch",
          str(metrics.get("method_actual")))
    check("7c fallback=True 且有 fallback_reason",
          bool(metrics.get("fallback")) and bool(metrics.get("fallback_reason")),
          str(metrics.get("fallback_reason"))[:80])
    text = " ".join(str(item) for item in (result.get("warnings") or []))
    check("7d 警告里明确说“均值比较 / 未验证 counts 不能做 QL”",
          ("均值比较" in text) and ("未验证" in text or "counts" in text), text[:160])


# ---------------------------------------------------------------------------
# [8] 配对设计
# ---------------------------------------------------------------------------


def test_paired_design():
    print("\n[8] 同一供体跨两组的配对设计")
    ctx = build_ctx(paired=True, donors_per_group=4)
    params = dict(DONOR_PARAMS, paired=True)
    result = run_de(ctx, params)
    metrics = result["summary"]["metrics"]
    check("8a paired=True 被识别并记录", bool(metrics.get("paired")) is True,
          str(metrics.get("paired")))
    check("8b design 里包含供体阻断项（配对模型）",
          "donor" in str(metrics.get("design") or "").lower()
          or "pair" in str(metrics.get("design") or "").lower(),
          str(metrics.get("design"))[:100])
    check("8c 配对下 n 仍按供体数计（4 对）",
          int(metrics.get("n_donors_in") or 0) == 4 and int(metrics.get("n_donors_out") or 0) == 4,
          "%s/%s" % (metrics.get("n_donors_in"), metrics.get("n_donors_out")))
    unpaired_ctx = build_ctx(paired=False)
    result2 = run_de(unpaired_ctx, params)
    metrics2 = result2["summary"]["metrics"]
    check("8d 无法配对时不得声称配对模型，并记录原因",
          bool(metrics2.get("paired")) is False
          and ("配对" in " ".join(str(w) for w in result2.get("warnings") or [])
               or "paired" in str(metrics2.get("fallback_reason") or "").lower()),
          str([w for w in result2.get("warnings") or [] if "配对" in str(w)])[:120])


def main() -> int:
    print("== P0-2 供体级统计（donor-level pseudobulk）回归 ==")
    for fn in (test_pseudobulk_sums_counts, test_full_integer_scan, test_donor_ql_primary,
               test_donor_min_cells_cache_key, test_missing_participant,
               test_single_cell_type, test_unverified_counts_fallback, test_paired_design):
        try:
            fn()
        except Exception:  # noqa: BLE001
            FAILED.append(fn.__name__)
            print("  FAIL %s 抛异常" % fn.__name__)
            traceback.print_exc()
    print("\n== 结果：%d 通过 / %d 失败 ==" % (len(PASSED), len(FAILED)))
    for item in FAILED:
        print("  x %s" % item)
    if os.environ.get("DS_SC_KEEP_FIXTURES") != "1":
        shutil.rmtree(TMP_ROOT, ignore_errors=True)
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
