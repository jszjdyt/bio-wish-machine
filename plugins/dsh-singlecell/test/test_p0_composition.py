#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P0-3 回归：细胞类型组成分析不能再拿「细胞」当重复。

反例（修复前实测）：
  1. ``_run_proportion_test`` 以细胞为独立观测做卡方，每细胞类型 p 也是细胞级；
  2. CI 是细胞计数的 Beta 后验（细胞越多区间越窄，但生物学 n 没变）；
  3. 每个细胞类型单独做 2×G 卡方，分母各自不同、且零计数格子被丢弃；
  4. 表/图/p/q/CI/n 不是同一张表：图用细胞计数占比、表用另一套富集组；
  5. 同一患者的细胞复制一遍，卡方 n 翻倍 → 假显著。

本测试驱动 ``sc_analyses_sets.run("proportion_test", ...)``：
  * [1] 样本层分析：同一分母、保留零、n=生物样本数；
  * [2] 结构性缺失（某类型只在 0/1 个组出现）必须显式区分，不得给 p=1 假装检验；
  * [3] p/q/CI/n 来自同一张 composition_stats 表；FDR 跨细胞类型校正；
  * [4] 复制同一患者细胞不得增加生物学 n；
  * [5] 配对设计（同供体跨两组）走配对模型；
  * [6] 无供体列时降级为细胞级探索性卡方并明确标注（向后兼容）；
  * [7] 图/表/指标一致（n_significant 与表里 q<0.05 行数相等）；
  * [8] 分母定义可回查：当前只有单一分母（total_cells=该样本全部合格细胞），
        分母定义与统计单位必须写在表注/metrics 里。backlog §2 要求的
        「NK/全部合格细胞 vs NK/全部合格免疫细胞」双分母尚未实现，属 P1 缺口；
        将来实现时必须同步改写分母声明，否则 [8] 会红。

运行：cd plugins/dsh-singlecell && python3 test/test_p0_composition.py
退出码 0 = 全部 PASS。
"""

from __future__ import annotations

import importlib.util
import math
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
TMP_ROOT = Path(tempfile.mkdtemp(prefix="sc-p0-composition-"))
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
        print("  PASS %s%s" % (name, ("  [" + str(detail)[:220] + "]") if detail else ""))
        return True
    FAILED.append(name)
    print("  FAIL %s%s" % (name, ("  [" + str(detail)[:220] + "]") if detail else ""))
    return False


def _load(module_name, filename):
    spec = importlib.util.spec_from_file_location(module_name, str(LIB / filename))
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


MODULE = _load("_p0_composition_analyses", "sc_analyses_sets.py")


def _log_norm(counts):
    counts = sp.csr_matrix(np.asarray(counts, dtype=np.float64))
    row_sum = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.zeros_like(row_sum)
    np.divide(1e4, row_sum, out=scale, where=row_sum > 0)
    matrix = (sp.diags(scale) @ counts).tocsr()
    matrix.data = np.log1p(matrix.data)
    matrix.eliminate_zeros()
    return matrix.astype(np.float32)


CELL_TYPES = ["T cell", "B cell", "NK cell", "Monocyte"]
# 组间组成差异明显：case 的 T cell 下降、NK/Monocyte 上升
PROPORTIONS = {
    "ctrl": {"T cell": 0.50, "B cell": 0.25, "NK cell": 0.15, "Monocyte": 0.10},
    "case": {"T cell": 0.20, "B cell": 0.25, "NK cell": 0.30, "Monocyte": 0.25},
}


def build_composition_ctx(n_donors=6, cells_per_sample=120, seed=5, paired=False,
                          duplicate_donor_cells=0, missing_type_only_ctrl=False,
                          with_donor_column=True, donor_nan=0):
    """每供体一个样本；样本内细胞类型按 PROPORTIONS 精确分配。"""
    rng = np.random.default_rng(seed)
    samples = []
    if paired:
        for i in range(n_donors):
            samples.append(("p%d" % (i + 1), i % 2 == 0 and "ctrl" or "case", "p%d" % (i + 1)))
            samples.append(("p%d_case" % (i + 1), "case", "p%d" % (i + 1)))
        # 每个供体同时出现在 ctrl 与 case：配对
        samples = [(("s%d_ctrl" % (i + 1)), "ctrl", "s%d" % (i + 1)) for i in range(n_donors)]
        samples += [(("s%d_case" % (i + 1)), "case", "s%d" % (i + 1)) for i in range(n_donors)]
    else:
        samples = [("ctrl_%d" % (i + 1), "ctrl", "dA%d" % (i + 1)) for i in range(n_donors)]
        samples += [("case_%d" % (i + 1), "case", "dB%d" % (i + 1)) for i in range(n_donors)]
    n_genes = 80
    genes = ["G%03d" % i for i in range(n_genes)]
    obs_rows = []
    count_rows = []
    for sample, group, donor in samples:
        props = dict(PROPORTIONS[group])
        if missing_type_only_ctrl and group != "ctrl":
            props.pop("Monocyte", None)   # case 组完全没有 Monocyte → 结构性缺失
            total = sum(props.values())
            props = {key: value / total for key, value in props.items()}
        n_cells = cells_per_sample
        # 加一点每样本抖动（总细胞数不变）：恒定比例会让两组方差为 0，p 值不可计算
        jitter = rng.normal(0.0, 3.5, size=len(props))
        raw = np.asarray([value * n_cells for value in props.values()], dtype=np.float64) + jitter
        raw = np.maximum(raw, 0.0)
        raw = raw / max(raw.sum(), 1.0) * n_cells
        type_counts = {key: int(math.floor(raw[index]))
                       for index, key in enumerate(props)}
        remainder = n_cells - sum(type_counts.values())
        order = np.argsort(-(raw - np.floor(raw)))
        for index in order[:max(0, remainder)]:
            type_counts[list(props)[index]] += 1
        for cell_type, count in type_counts.items():
            for _ in range(max(0, count)):
                row = rng.poisson(2.0, size=n_genes).astype(np.float64)
                count_rows.append(row)
                obs_rows.append({
                    "sample": sample,
                    "group": group,
                    "cell_type": cell_type,
                    "donor": donor if with_donor_column else None,
                })
        # 结构性缺失：case 组 Monocyte 0 个（无行），ctrl 组有
    if duplicate_donor_cells:
        first_donor = pd.DataFrame(obs_rows)
        mask = first_donor["donor"] == first_donor["donor"].iloc[0]
        index = np.flatnonzero(mask.to_numpy())
        for position in index:
            obs_rows.append(dict(obs_rows[position]))
            count_rows.append(count_rows[position].copy())
    counts = np.vstack(count_rows) if count_rows else np.zeros((1, n_genes))
    obs = pd.DataFrame(obs_rows, index=["cell_%05d" % i for i in range(len(obs_rows))])
    if with_donor_column and donor_nan:
        bad = obs["donor"].dropna().unique()[:donor_nan]
        obs.loc[obs["donor"].isin(bad).to_numpy(), "donor"] = np.nan
    if not with_donor_column:
        # 去掉所有生物学重复线索（sample/donor），否则会被自动识别成样本列
        obs = obs.drop(columns=[column for column in ("donor", "sample") if column in obs.columns])
    matrix = sp.csr_matrix(counts.astype(np.float32))
    return {
        "dataset_id": "ds-composition",
        "dataset_version": 0,
        "dataset_revision": 0,
        "content_revision": "composition-fixture",
        "name": "composition-fixture",
        "n_obs": int(counts.shape[0]),
        "n_vars": n_genes,
        "X": _log_norm(counts),
        "X_kind": "lognorm",
        "counts": matrix,
        "obs": obs,
        "var_names": genes,
        "obsm": {},
        "clusters": None,
        "group_keys": [column for column in obs.columns if obs[column].dtype == object],
        "engine_mode": "builtin",
        "rng": np.random.default_rng(0),
    }


def run_proportion(ctx, params=None):
    params = params or {"group": "group", "cell_type": "cell_type", "min_cells": 1}
    return MODULE.run("proportion_test", ctx, params, str(OUT / ("comp-%d" % len(os.listdir(OUT)))))


def _table(result, table_id):
    return next((table for table in result.get("tables") or [] if table.get("id") == table_id), None)


def _row_keys(table):
    return [str(column.get("key")) for column in table.get("columns") or []]


# ---------------------------------------------------------------------------
# [1] 样本层分析：同分母 / 保留零 / n=样本数
# ---------------------------------------------------------------------------


def test_sample_level_composition():
    print("\n[1] 样本层组成分析（同一分母、保留零、只有 2×6 个样本）")
    ctx = build_composition_ctx()
    result = run_proportion(ctx)
    metrics = result["summary"]["metrics"]
    unit = str(metrics.get("statistical_unit") or "")
    check("1a 统计单位是生物样本（sample/donor），不是 cell",
          unit in ("sample", "donor"), unit)
    check("1b 两组样本数 n=6/6 与设计一致",
          dict(metrics.get("n_samples_per_group") or {}) == {"ctrl": 6, "case": 6},
          str(metrics.get("n_samples_per_group")))
    by_sample = _table(result, "composition_by_sample")
    check("1c 有 composition_by_sample 表（每样本 × 细胞类型计数）", by_sample is not None,
          str([t.get("id") for t in result.get("tables") or []]))
    if by_sample is not None:
        keys = _row_keys(by_sample)
        rows = by_sample["rows"]
        check("1d 每样本列出全部细胞类型（零也保留）",
              len(keys) >= 1 + len(CELL_TYPES) + 2 and all(isinstance(row, dict) for row in rows),
              str(keys[:8]))
        totals = [float(row.get("total_cells") or 0) for row in rows]
        check("1e 每样本 total_cells 相同（同一分母）",
              len(set(int(value) for value in totals)) == 1, str(sorted(set(totals))[:4]))
        type_keys = [key for key in keys if key.startswith("n_")]
        if type_keys:
            sums = [sum(float(row.get(key) or 0) for key in type_keys) for row in rows]
            check("1f 每样本各类型计数之和 = total_cells（零未被丢弃）",
                  all(abs(one - two) < 1e-9 for one, two in zip(sums, totals)),
                  "sums=%s totals=%s" % (sums[:3], totals[:3]))
    stats = _table(result, "composition_stats")
    check("1g 有 composition_stats 表（p/q/CI/n 同一张表）", stats is not None,
          str([t.get("id") for t in result.get("tables") or []]))
    if stats is not None:
        keys = set(_row_keys(stats))
        required = {"cell_type", "n1", "n2", "pct1", "pct2", "diff", "ci_low", "ci_high",
                    "pvalue", "pval_adj"}
        check("1h 表里同时含 n/p/q/CI/组占比", required <= keys, str(sorted(required - keys)))
        rows = stats["rows"]
        check("1i 每个细胞类型一行，n1=n2=6（供体/样本数，不是细胞数）",
              len(rows) == len(CELL_TYPES)
              and all(int(row.get("n1") or 0) == 6 and int(row.get("n2") or 0) == 6 for row in rows),
              str([(row.get("cell_type"), row.get("n1"), row.get("n2")) for row in rows][:4]))


# ---------------------------------------------------------------------------
# [2] 结构性缺失
# ---------------------------------------------------------------------------


def test_structural_missing():
    print("\n[2] 结构性缺失（某细胞类型只在 0/1 个组出现）")
    ctx = build_composition_ctx(missing_type_only_ctrl=True)
    result = run_proportion(ctx)
    stats = _table(result, "composition_stats")
    check("2a 结构缺失时仍产出 composition_stats", stats is not None, "")
    if stats is not None:
        rows = {str(row.get("cell_type")): row for row in stats["rows"]}
        mono = rows.get("Monocyte")
        check("2b 只在 ctrl 出现的类型被显式标注（status/note 非空）",
              mono is not None and str(mono.get("status") or "") not in ("", "tested", "None"),
              str(mono and {k: mono.get(k) for k in ("cell_type", "status", "pvalue")}))
        check("2c 这类类型不得给出伪造 p 值（pvalue 为 None）",
              mono is not None and mono.get("pvalue") is None,
              str(mono and mono.get("pvalue")))
    warnings = " ".join(str(item) for item in (result.get("warnings") or []))
    check("2d 警告里说明结构性缺失/只在单组出现",
          any(token in warnings for token in ("结构性", "只在", "仅一组", "单一组", "structur")),
          warnings[:180])


# ---------------------------------------------------------------------------
# [3] 同一张表 + BH
# ---------------------------------------------------------------------------


def test_same_table_and_fdr():
    print("\n[3] p/q/CI/n 同一张表 + 跨细胞类型 BH")
    ctx = build_composition_ctx()
    result = run_proportion(ctx)
    stats = _table(result, "composition_stats")
    if stats is None:
        check("3a composition_stats 存在", False, "")
        return
    rows = stats["rows"]
    check("3a 每行 pvalue/pval_adj/CI/n 齐全（同一行）",
          all(all(row.get(key) is not None for key in ("pvalue", "pval_adj", "ci_low", "ci_high", "n1", "n2"))
              for row in rows),
          str([{k: row.get(k) for k in ("cell_type", "pvalue", "pval_adj")} for row in rows][:2]))
    tested = [row for row in rows if row.get("pvalue") is not None]
    adjusted = sorted(float(row["pval_adj"]) for row in tested)
    raw = sorted(float(row["pvalue"]) for row in tested)
    check("3b BH 校正后至少不小于原始 p（q>=p）",
          all(q >= p - 1e-12 for p, q in zip(raw, adjusted)), "%s vs %s" % (raw, adjusted))
    check("3c FDR 跨细胞类型做（不是每类各自 1 个检验）",
          len(tested) >= 2 and bool(metrics_fdr_used(result)), str(len(tested)))
    tcell = next((row for row in rows if str(row.get("cell_type")) == "T cell"), None)
    check("3d T cell 组成差异被检出（FDR<0.05）",
          tcell is not None and tcell.get("pval_adj") is not None and float(tcell["pval_adj"]) < 0.05,
          str(tcell and {k: tcell.get(k) for k in ("pct1", "pct2", "pvalue", "pval_adj")}))
    metrics = result["summary"]["metrics"]
    group1 = str(metrics.get("group1") or "")
    expected_high = "ctrl" if group1 == "ctrl" else "case"
    if tcell is not None and group1 == "ctrl":
        direction_ok = float(tcell.get("pct1") or 0) > float(tcell.get("pct2") or 0)
    elif tcell is not None:
        direction_ok = float(tcell.get("pct2") or 0) > float(tcell.get("pct1") or 0)
    else:
        direction_ok = False
    check("3e 差异方向正确（ctrl 的 T cell 占比高于 case）", direction_ok,
          "%s pct1=%s pct2=%s" % (group1, tcell and tcell.get("pct1"), tcell and tcell.get("pct2")))


def metrics_fdr_used(result):
    metrics = result["summary"]["metrics"]
    return bool(metrics.get("fdr_across_cell_types")) or int(metrics.get("n_tested") or 0) >= 2


# ---------------------------------------------------------------------------
# [4] 复制同一患者细胞不得增加生物学 n
# ---------------------------------------------------------------------------


def test_duplicated_cells_do_not_inflate_n():
    print("\n[4] 复制同一患者细胞不得增加生物学 n")
    base = run_proportion(build_composition_ctx())
    dup = run_proportion(build_composition_ctx(duplicate_donor_cells=1))
    base_n = dict(base["summary"]["metrics"].get("n_samples_per_group") or {})
    dup_n = dict(dup["summary"]["metrics"].get("n_samples_per_group") or {})
    base_cells = int(base["summary"]["metrics"].get("n_cells_total") or 0)
    dup_cells = int(dup["summary"]["metrics"].get("n_cells_total") or 0)
    check("4a 细胞数确实翻倍（复制生效）", dup_cells > base_cells, "%d → %d" % (base_cells, dup_cells))
    check("4b 生物学 n（每组件数）不因复制细胞而增加", base_n == dup_n, "%s vs %s" % (base_n, dup_n))
    base_stats = _table(base, "composition_stats")
    dup_stats = _table(dup, "composition_stats")
    if base_stats and dup_stats:
        base_p = {row["cell_type"]: row.get("pvalue") for row in base_stats["rows"]}
        dup_p = {row["cell_type"]: row.get("pvalue") for row in dup_stats["rows"]}
        same = all(
            (base_p[key] is None and dup_p.get(key) is None)
            or (base_p[key] is not None and dup_p.get(key) is not None
                and abs(float(base_p[key]) - float(dup_p[key])) < 1e-9)
            for key in base_p
        )
        check("4c 复制细胞后 p 值不变（n 没被当成细胞数）", same,
              str({key: (base_p[key], dup_p.get(key)) for key in list(base_p)[:3]}))


# ---------------------------------------------------------------------------
# [5] 配对设计
# ---------------------------------------------------------------------------


def test_paired_model():
    print("\n[5] 同供体跨两组的配对模型")
    ctx = build_composition_ctx(paired=True)
    result = run_proportion(ctx)
    metrics = result["summary"]["metrics"]
    check("5a paired=True 被识别并记录", bool(metrics.get("paired")) is True,
          str(metrics.get("paired")))
    check("5b 配对模型下 n 是配对数 6/6",
          dict(metrics.get("n_samples_per_group") or {}) == {"ctrl": 6, "case": 6},
          str(metrics.get("n_samples_per_group")))
    check("5c 模型字段说明走配对检验",
          any(token in str(metrics.get("test") or "").lower() for token in ("paired", "配对")),
          str(metrics.get("test")))


# ---------------------------------------------------------------------------
# [6] 无供体列 → 细胞级探索性并明确标注
# ---------------------------------------------------------------------------


def test_cell_level_fallback():
    print("\n[6] 无生物学重复列时降级为细胞级探索性分析")
    ctx = build_composition_ctx(with_donor_column=False)
    result = run_proportion(ctx)
    metrics = result["summary"]["metrics"]
    check("6a 统计单位明确标为 cell", str(metrics.get("statistical_unit") or "") == "cell",
          str(metrics.get("statistical_unit")))
    check("6b 明确标注细胞级结果不能作患者组间结论",
          any(token in " ".join(str(w) for w in result.get("warnings") or [])
              for token in ("细胞级", "探索", "不能", "重复")),
          str((result.get("warnings") or [""])[:2])[:180])


# ---------------------------------------------------------------------------
# [7] 图/表/指标一致
# ---------------------------------------------------------------------------


def test_figure_table_metrics_agree():
    print("\n[7] 图 / 表 / 指标对应同一张表")
    ctx = build_composition_ctx()
    result = run_proportion(ctx)
    stats = _table(result, "composition_stats")
    metrics = result["summary"]["metrics"]
    if stats is None:
        check("7a composition_stats 存在", False, "")
        return
    significant = [row for row in stats["rows"]
                   if row.get("pval_adj") is not None and float(row["pval_adj"]) < 0.05]
    check("7a metrics.n_significant 与表里 FDR<0.05 行数一致",
          int(metrics.get("n_significant") or 0) == len(significant),
          "%s vs %s" % (metrics.get("n_significant"), len(significant)))
    figure_rows = result["summary"].get("figure_data") or metrics.get("figure_rows")
    check("7b 图数据来自同一张表（figure_rows 与表行数一致）",
          figure_rows is not None and len(figure_rows) == len(stats["rows"]),
          None if figure_rows is None else "%d vs %d" % (len(figure_rows), len(stats["rows"])))
    figures = result.get("figures") or []
    check("7c 仍产出 PNG 图（composition 堆叠柱）", len(figures) >= 1,
          str([figure.get("id") for figure in figures]))


# ---------------------------------------------------------------------------
# [8] 分母定义可回查（P0-3 的 P1 尾巴：双分母尚未实现）
# ---------------------------------------------------------------------------


def test_denominator_contract_is_recorded():
    """把「当前单一分母口径」钉住，避免将来改分母时静默错标。

    backlog §2 的细化要求是「分别定义 NK/全部合格细胞 与 NK/全部合格免疫细胞」。
    当前实现只有**单一分母**（total_cells = 该样本全部合格细胞），属于已记录的
    P1 实现缺口。这里不去写会红的用例，而是要求：分母定义必须写进表注、统计单位
    必须写进 metrics。将来若实现 params.denominator=immune 却忘了同步改写分母
    声明，本组断言会红，从而防止「换了分母但文档仍写全部细胞」。
    """
    print("\n[8] 分母定义可回查（单一分母口径的显式声明）")
    ctx = build_composition_ctx()
    result = run_proportion(ctx)
    metrics = result["summary"]["metrics"]
    by_sample = _table(result, "composition_by_sample")
    stats = _table(result, "composition_stats")
    sample_note = str((by_sample or {}).get("note") or "")
    stats_note = str((stats or {}).get("note") or "")
    check("8a composition_by_sample 表注声明分母口径（全部细胞 / 同一分母）",
          by_sample is not None and "全部细胞" in sample_note and "同一分母" in sample_note,
          sample_note)
    check("8b composition_stats 表注声明统计单位是生物样本",
          stats is not None and "生物样本" in stats_note, stats_note)
    check("8c metrics 记录统计单位与供体列，便于回查分母定义来源",
          str(metrics.get("statistical_unit") or "") == "sample"
          and bool(metrics.get("donor_column")),
          "%s / %s" % (metrics.get("statistical_unit"), metrics.get("donor_column")))
    if by_sample is not None:
        type_keys = [key for key in _row_keys(by_sample) if key.startswith("n_")]
        rows = by_sample["rows"]
        consistent = bool(type_keys) and all(
            abs(sum(float(row.get(key) or 0) for key in type_keys)
                - float(row.get("total_cells") or 0)) < 1e-9
            for row in rows
        )
        check("8d 分母与分子同源（各类型计数之和 = total_cells，未中途换分母）",
              consistent, "n_types=%d" % len(type_keys))


def main() -> int:
    print("== P0-3 样本层组成分析（composition）回归 ==")
    for fn in (test_sample_level_composition, test_structural_missing, test_same_table_and_fdr,
               test_duplicated_cells_do_not_inflate_n, test_paired_model,
               test_cell_level_fallback, test_figure_table_metrics_agree,
               test_denominator_contract_is_recorded):
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
