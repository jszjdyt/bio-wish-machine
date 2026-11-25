#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — P0-4 回归测试：ORA 校正家族 + 免疫分数/比例

反例来源（BIO_P0_BACKLOG §1 GEO Bulk 第 4 条）：
  * ``lib/bg_enrich.py:121-125,166-199`` 的 ``ora()`` 在 BH 之前用 ``if k < 2: continue``
    剔除「命中不足 2 个」的通路，校正家族被缩小：100 条预先合格通路的合成对照里
    q≈0.0155，按完整家族应为 1（本测试用可手算的构造验证家族被缩小 100 倍）；
  * ``lib/bg_immune.py:506-540,898-905`` 的 ``_build_proportions()`` 在没有天然比例
    方法时，把 ssGSEA/MCP-counter 的得分行归一化成「比例」——分数量纲被伪装成
    有分母的比例。

要求：
  * 所有预先合格（映射后大小达标）的通路先做 BH，再按展示阈值筛；universe = 实际
    受检且可映射基因；结果保存库来源/版本/物种/映射率/校正家族；
  * 免疫分数与比例分开输出；比例图只接含明确分母的比例数据。

用法：cd plugins/dsh-geobulk && python3 test/test_p0_4_enrich_immune.py
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="bg-mpl-"))

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bg_enrich  # noqa: E402
import bg_immune  # noqa: E402
from bg_common import EngineError  # noqa: E402

OUTDIR = tempfile.mkdtemp(prefix="bg-p0-4-")
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


def _hypergeom_sf(k: int, K: int, n: int, N: int) -> float:
    from scipy.stats import hypergeom

    return float(hypergeom.sf(k - 1, N, K, n))


def _bh(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(q, 0, 1)
    return out


# ---------------------------------------------------------------------------
# A. ORA：完整校正家族
# ---------------------------------------------------------------------------

SET_SIZE = 5  # 每条通路 5 个基因（>= MIN_SET_SIZE）


def _family_fixture(n_terms=100, n_universe=1000, n_query=5, seed=5):
    """100 条预先合格通路：99 条各命中 1 个 query 基因，1 条命中 2 个。

    universe 1000；每条通路 5 个基因（1~2 个 query 基因 + 互不重叠背景），
    保证映射后 K=5 ≥ MIN_SET_SIZE，全部「预先合格」。
    """
    universe = [f"G{i:04d}" for i in range(n_universe)]
    query = universe[:n_query]
    background = universe[n_query:]
    sets: dict[str, list[str]] = {}
    names: dict[str, str] = {}
    cursor = 0
    for t in range(n_terms - 1):
        members = [query[t % n_query]]
        members += list(background[cursor:cursor + (SET_SIZE - 1)])
        cursor += SET_SIZE - 1
        key = f"TERM_{t:03d}"
        sets[key] = members
        names[key] = f"synthetic term {t}"
    hit = "TERM_HIT"
    sets[hit] = [query[0], query[1]] + list(background[cursor:cursor + (SET_SIZE - 2)])
    names[hit] = "synthetic hit term"
    return query, universe, sets, names


def test_ora_family():
    print("\n[A] ORA 校正家族（先全家族 BH，再筛展示）")
    query, universe, sets, names = _family_fixture()
    frame = bg_enrich.ora(query, universe, sets, names, "SYNTHETIC")
    before_rows = len(frame)
    hit = frame[frame["ID"] == "TERM_HIT"] if not frame.empty else frame
    before_q = None if hit.empty else float(hit.iloc[0]["p.adjust"])
    BEFORE["A1"] = f"展示 {before_rows} 条；TERM_HIT q={before_q}（家族被 k<2 缩小）"

    # 手算：完整家族 = 100 条（99 条 k=1 + 1 条 k=2），N=1000、K=SET_SIZE、n=5
    p_other = _hypergeom_sf(1, SET_SIZE, 5, 1000)
    p_hit = _hypergeom_sf(2, SET_SIZE, 5, 1000)
    full_family = np.concatenate([np.full(99, p_other), [p_hit]])
    expected_q = _bh(full_family)
    expected_hit_q = float(expected_q[-1])   # TERM_HIT 的 p 最小，BH 后排第 1
    check("A1 展示表仍只保留命中 ≥2 的通路（展示阈值在 BH 之后）",
          len(frame) == 1 and str(frame.iloc[0]["ID"]) == "TERM_HIT",
          f"rows={len(frame)} ids={list(frame['ID'])[:3]}")
    family_size = int(getattr(frame, "attrs", {}).get("n_tested") or 0)
    check("A2 记录校正家族大小 = 全部预先合格通路（100）", family_size == 100,
          f"n_tested={family_size}")
    after_q = float(frame.iloc[0]["p.adjust"]) if len(frame) else None
    BEFORE["A2"] = f"TERM_HIT q={after_q}（完整家族手算 {expected_hit_q:.6f}）"
    check("A3 TERM_HIT 的 q 值与完整家族 BH 手算结果一致",
          after_q is not None and abs(after_q - expected_hit_q) < 1e-12,
          f"q={after_q} expected={expected_hit_q:.6f}")
    # 旧实现（家族被缩小成 1 条）会让 q = p_hit；完整家族应约放大 100 倍
    single_family_q = p_hit
    BEFORE["A4"] = f"旧实现（家族=1）q≈{single_family_q:.8f}；完整家族 q={after_q:.6f}"
    check("A4 完整家族使 q 相对「只算命中≥2 的家族」放大约 100 倍",
          after_q is not None and after_q > single_family_q * 50,
          f"single={single_family_q:.8f} after={after_q:.6f}")
    check("A5 结果保存家族/映射率等溯源字段",
          all(k in getattr(frame, "attrs", {}) for k in
              ("n_tested", "universe_size", "query_size", "library_hash", "species")),
          str(sorted(getattr(frame, "attrs", {}) or {}))[:120])

    # A6 universe = 实测基因 ∩ 库可映射基因
    measured = [f"H{i:04d}" for i in range(1200)]   # 实测 1200 个，其中 1000 个能被库映射
    lib_genes = measured[:1000]
    universe_sets = {f"LIB_{i:02d}": lib_genes[i * 100:(i + 1) * 100] for i in range(10)}
    de = pd.DataFrame({"gene": measured, "logFC": np.ones(1200)})
    ctx = {"genes": measured, "X": np.zeros((1, 1200))}
    info = bg_enrich._universe(ctx, de, universe_sets)
    universe_set = info[0] if isinstance(info, tuple) else info
    universe_meta = info[1] if isinstance(info, tuple) else {}
    BEFORE["A6"] = f"universe={len(universe_set)}（库可映射 1000；修复前直接用 1200 个实测基因）"
    check("A6 universe = 实际受检且可映射基因（1000，而不是 1200）",
          len(universe_set) == 1000, f"universe={len(universe_set)}")
    check("A7 记录 ID 映射率", abs(float(universe_meta.get("mapping_rate", 0)) - 1000 / 1200) < 1e-6,
          f"rate={universe_meta.get('mapping_rate')}")

    # A8 步骤级：metrics/store 带家族/库版本/物种
    step_ctx = _enrich_ctx()
    payload = bg_enrich.run_enrich(step_ctx, {"db": ["GO_BP"], "method": "ora"}, os.path.join(OUTDIR, "a8"))
    metrics = payload["summary"]["metrics"]
    BEFORE["A8"] = f"metrics.ora_families={metrics.get('ora_families')!r}"
    families = metrics.get("ora_families")
    ok = isinstance(families, list) and families and all(
        {"family", "n_tested", "universe_size", "mapping_rate", "library_source"} <= set(f)
        for f in families)
    check("A8 步骤 metrics 保存每个校正家族的 n_tested/universe/映射率/库来源", bool(ok),
          str(families)[:160])
    check("A8b 家族非空（确实有预先合格通路进入 BH）",
          bool(families) and sum(int(f.get("n_tested") or 0) for f in families) > 0,
          f"n_tested_sum={sum(int(f.get('n_tested') or 0) for f in families)}")
    check("A9 结果标注库版本（内容 hash）、物种",
          bool(metrics.get("library_version")) and metrics.get("species") == "human",
          f"version={metrics.get('library_version')} species={metrics.get('species')}")
    store = ((payload.get("store") or {}).get("enrich") or {})
    check("A10 store 保存家族与库溯源", bool(store.get("ora_families")) and bool(store.get("library_version")),
          f"families={len(store.get('ora_families') or [])} version={store.get('library_version')}")
    notes = " ".join(payload.get("notes") or [])
    check("A11 方法说明写明「先全家族 BH 再筛展示」", "家族" in notes and "BH" in notes, notes[:120])
    columns = payload["tables"][0]["columns"] if payload.get("tables") else []
    col_keys = [c.get("key") for c in columns]
    check("A12 导出表带 Family/FamilySize 列（可核对校正家族）",
          "Family" in col_keys and "FamilySize" in col_keys, str(col_keys))


def _enrich_ctx(n_extra=800, seed=9):
    """GO_BP 库基因全部纳入表达矩阵 + 额外真实基因，保证 ORA 家族非空。"""
    rng = np.random.default_rng(seed)
    pool = [line.strip() for line in (ROOT / "lib" / "demo_gene_pool.txt").read_text(encoding="utf-8").splitlines()
            if line.strip()]
    prepared = bg_enrich._prepare_sets("GO_BP", workdir=OUTDIR)
    library: set[str] = set()
    for members in prepared["sets"].values():
        library.update(str(g) for g in members)
    lib_genes = [g for g in pool if g in library]
    extra = list(rng.choice([g for g in pool if g not in library], size=min(n_extra, len(pool)), replace=False))
    genes = lib_genes + extra
    n_genes = len(genes)
    samples = [f"S{i}" for i in range(8)]
    X = rng.normal(7, 1, (8, n_genes))
    logfc = rng.normal(0, 0.6, n_genes)
    logfc[:40] += 1.5
    p = np.clip(rng.uniform(0.2, 1.0, n_genes), 1e-6, 1)
    p[:40] = 1e-4
    de = pd.DataFrame({
        "gene": genes, "logFC": logfc, "AveExpr": X.mean(0), "t": logfc,
        "P.Value": p, "adj.P.Val": np.clip(p * 2, 0, 1), "B": np.zeros(n_genes),
    })
    de["regulation"] = np.where(de["adj.P.Val"] <= 0.05, "up", "ns")
    return {
        "dataset_id": "ds_p0_4",
        "name": "synthetic",
        "genes": genes, "samples": samples, "X": X, "X_kind": "lognorm",
        "sample_info": pd.DataFrame({"group": ["case"] * 4 + ["control"] * 4}, index=samples),
        "group_key": "group", "groups": ["case", "control"],
        "results": {}, "de": {"table": de, "case": "case", "control": "control"},
        "workdir": OUTDIR,
    }


# ---------------------------------------------------------------------------
# B. 免疫：分数 ≠ 比例
# ---------------------------------------------------------------------------

def _score_frames(n=12, seed=3):
    rng = np.random.default_rng(seed)
    samples = [f"S{i:02d}" for i in range(n)]
    ssgsea = pd.DataFrame(
        rng.normal(3000, 800, (n, 6)),
        index=samples, columns=[f"Immune_{i}" for i in range(6)],
    )
    mcp = pd.DataFrame(
        rng.normal(500, 120, (n, 5)),
        index=samples, columns=[f"MCP_{i}" for i in range(5)],
    )
    return samples, {"ssgsea": ssgsea, "mcpcounter": mcp}


def test_immune_scores():
    print("\n[B] 免疫分数与比例分开")
    samples, per_method = _score_frames()
    warnings: list[str] = []
    notes: list[str] = []
    frame, meta = bg_immune._build_proportions(per_method, "relative", warnings, notes)
    row_sums = frame.sum(axis=1).values
    before_kind = (meta or {}).get("kind", "<未记录>")
    BEFORE["B1"] = (f"kind={before_kind}; 行和范围 {row_sums.min():.4f}–{row_sums.max():.4f}；"
                    f"derived={(meta or {}).get('derived')}")
    check("B1 只有 ssGSEA/MCP 得分时输出 kind=score（不是比例）",
          str((meta or {}).get("kind")) == "score", f"kind={(meta or {}).get('kind')}")
    check("B2 得分表不再被行归一化成「比例」（行和不应≈1）",
          not np.allclose(row_sums, 1.0, atol=1e-3), f"row_sums={row_sums[:3]}")
    check("B3 得分表与原始 ssGSEA 得分一致（未做 min-max+行归一）",
          np.allclose(frame.values, per_method["ssgsea"].values),
          f"max|diff|={np.nanmax(np.abs(frame.values - per_method['ssgsea'].values)):.3g}")
    check("B4 记录「无明确分母」", (meta or {}).get("proportion_available") is False,
          f"proportion_available={(meta or {}).get('proportion_available')}")

    # 天然比例方法：保持比例语义
    prop_frame = pd.DataFrame(
        np.random.default_rng(7).dirichlet(np.ones(6), size=10),
        index=[f"S{i}" for i in range(10)], columns=[f"C{j}" for j in range(6)],
    )
    pframe, pmeta = bg_immune._build_proportions({"cibersort": prop_frame}, "relative", [], [])
    BEFORE["B2"] = f"cibersort kind={(pmeta or {}).get('kind')} 行和={pframe.sum(axis=1).values[:3]}"
    check("B5 有天然比例方法时 kind=proportion 且行和≈1",
          str((pmeta or {}).get("kind")) == "proportion"
          and np.allclose(pframe.sum(axis=1).values, 1.0, atol=1e-6),
          f"kind={(pmeta or {}).get('kind')}")

    # B6 步骤级：score-only 时表/图/指标都要分开（强制 Python 回退，保证只有得分方法）
    ctx = _immune_ctx()
    original_r = bg_immune._run_r_immune
    bg_immune._run_r_immune = lambda *a, **k: (_ for _ in ()).throw(EngineError("no_r", "测试强制回退"))
    try:
        payload = bg_immune.run_immune(ctx, {"methods": ["estimate", "ssgsea", "mcpcounter"]},
                                       os.path.join(OUTDIR, "b6"))
    finally:
        bg_immune._run_r_immune = original_r
    metrics = payload["summary"]["metrics"]
    table_ids = [t.get("id") for t in payload.get("tables") or []]
    store = ((payload.get("store") or {}).get("immune") or {})
    BEFORE["B6"] = (f"metrics.table_kind={metrics.get('table_kind')}; tables={table_ids}; "
                    f"store.table_kind={store.get('table_kind')}")
    check("B6 score-only 时 metrics 标注 table_kind=score", metrics.get("table_kind") == "score",
          f"table_kind={metrics.get('table_kind')}")
    check("B7 score-only 时不再输出「免疫细胞比例」表，只输出特征/得分表",
          "immune_proportions" not in table_ids
          and any(tid in table_ids for tid in ("immune_feature_scores", "immune_scores")),
          str(table_ids))
    check("B8 store 标记 proportion_available=False + table_kind=score",
          store.get("proportion_available") is False and store.get("table_kind") == "score",
          f"available={store.get('proportion_available')} kind={store.get('table_kind')}")
    fig_ids = [str(f.get("id") or f.get("name") or "") for f in payload.get("figures") or []]
    check("B9 score-only 时不生成「组成比例」堆叠图", not any("immune_stack" in f for f in fig_ids),
          str(fig_ids))
    bullets = " ".join(payload["summary"].get("bullets") or [])
    check("B10 摘要明确是得分、不写成「比例/行和≈1」",
          "得分" in bullets and "行和" not in bullets and "不是细胞比例" in bullets,
          bullets[:140])
    check("B11 记录免疫表类型与分母来源",
          "proportion_denominator" in metrics,
          f"denominator={metrics.get('proportion_denominator')}")

    # B12 比例图闸门：score-only 的 store 不能画组成堆叠图
    import bg_plots

    score_ctx = dict(ctx)
    score_ctx["results"] = {"immune": dict(store)}
    raised = None
    try:
        bg_plots.plot_immune_stack(score_ctx, {}, os.path.join(OUTDIR, "b12"))
    except EngineError as error:
        raised = error
    BEFORE["B12"] = f"plot_immune_stack(score store) raised={None if raised is None else raised.code}"
    message = str(getattr(raised, "message", ""))
    check("B12 score-only 时组成堆叠图被闸门拦下（不把分数画成比例）",
          raised is not None and raised.code == "prerequisite" and "分母" in message,
          message[:100])

    # B13 真比例数据仍可画
    prop_store = {
        "table_kind": "proportion",
        "proportion_available": True,
        "table": bg_immune._frame_records(pframe),
        "scores": [],
    }
    prop_ctx = dict(ctx)
    prop_ctx["results"] = {"immune": prop_store}
    figure = bg_plots.plot_immune_stack(prop_ctx, {}, os.path.join(OUTDIR, "b13"))
    check("B13 含明确分母的比例数据仍能画组成图",
          bool(figure) and os.path.exists(str(figure.get("file"))), str(figure.get("file")))


def _immune_ctx(n=24, seed=11):
    """用仓库内置 demo 基因池造免疫数据：必须含真实 marker 基因，Python 回退才可算。"""
    rng = np.random.default_rng(seed)
    pool = [line.strip() for line in (ROOT / "lib" / "demo_gene_pool.txt").read_text(encoding="utf-8").splitlines()
            if line.strip()]
    pool_upper = {g.upper() for g in pool}
    markers: set[str] = set()
    for members in (bg_immune.IMMUNE_28 or {}).values():
        markers.update(str(g).upper() for g in members)
    for panel in bg_immune.MCP_PANELS.values():
        markers.update(str(g).upper() for g in panel)
    must = [g for g in pool if g.upper() in markers][:200]
    rest = list(rng.choice([g for g in pool if g.upper() not in markers], size=300, replace=False))
    genes = must + rest
    samples = [f"P{i:02d}" for i in range(n)]
    X = rng.normal(7, 1, (n, len(genes)))
    info = pd.DataFrame({"group": ["case"] * (n // 2) + ["control"] * (n - n // 2)}, index=samples)
    return {
        "dataset_id": "ds_immune",
        "name": "synthetic",
        "genes": genes, "samples": samples, "X": X, "X_kind": "lognorm",
        "sample_info": info, "group_key": "group", "groups": ["case", "control"],
        "results": {}, "workdir": OUTDIR,
    }


def main() -> int:
    print("== dsh-geobulk P0-4 ORA 家族 + 免疫分数/比例 ==")
    test_ora_family()
    test_immune_scores()
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
