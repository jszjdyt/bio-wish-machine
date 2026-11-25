#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — P0-5 回归测试：结果溯源 / 旧结果不得冒充当前结果

反例来源（BIO_P0_BACKLOG §1 GEO Bulk 第 5 条）：
  * ``lib/bg_registry.py`` 的 ``_de_consistency`` 只写进 Methods 描述，
    ``export_dataset()`` 直接读 ``dataset.results`` 旧缓存 → 改组/改矩阵后旧 DE 仍能导出；
  * 修订版 ID 不含元数据 + 原始层 + 协变量，无法判断「同一版本」；
  * 样本集合比较把「三组取二」误判成 stale（store 只含 8 个样本，数据集有 12 个）。

本测试要求：
  * 修订版 ID = 元数据 + 原始层（X/counts）+ 协变量；
  * 每个结果缓存记录 revision_id / result_fingerprint / run_id / 实际纳入样本；
  * 导出与出图前校验同一版本与实际纳入样本，旧缓存不得冒充当前结果；
  * 三组取二不被误判 stale；旧 run_id 可检索并导出历史版本。

用法：cd plugins/dsh-geobulk && python3 test/test_p0_5_result_provenance.py
"""

from __future__ import annotations

import json
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

import bg_registry  # noqa: E402
import bg_state  # noqa: E402
from bg_common import EngineError  # noqa: E402

OUTDIR = tempfile.mkdtemp(prefix="bg-p0-5-")
PASSED: list[str] = []
FAILED: list[str] = []
BEFORE: dict[str, str] = {}
_COUNTER = {"n": 0}


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


def _engine() -> dict:
    return {"python": "3.14.2", "r": {"available": True, "version": "4.4.2"}, "packages": {}}


def _pool_genes(n: int) -> list[str]:
    pool = [line.strip() for line in (ROOT / "lib" / "demo_gene_pool.txt").read_text(encoding="utf-8").splitlines()
            if line.strip()]
    return pool[:n]


def make_dataset(*, n_case=6, n_ctrl=6, n_genes=300, groups=None, clinical=False, seed=11,
                 extra_meta=None) -> "bg_state.Dataset":
    """构造计数矩阵数据集（真实基因 symbol，可跑 DE/富集/生存）。"""
    rng = np.random.default_rng(seed)
    genes = _pool_genes(n_genes)
    samples = [f"T{i}" for i in range(n_case)] + [f"N{i}" for i in range(n_ctrl)]
    base = rng.uniform(80, 400, n_genes)
    X = np.zeros((n_case + n_ctrl, n_genes))
    for i in range(n_case + n_ctrl):
        mu = base * rng.uniform(0.9, 1.1)
        if i < n_case:
            mu[:10] *= 2.5
        X[i] = rng.negative_binomial(4.0, 4.0 / (4.0 + mu))
    label = groups or (["tumor"] * n_case + ["normal"] * n_ctrl)
    info = pd.DataFrame({"group": label}, index=samples)
    if clinical:
        info["os_time"] = np.linspace(10, 90, len(samples))
        info["os_event"] = [1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1][: len(samples)]
    if extra_meta:
        for column, values in extra_meta.items():
            info[column] = values
    _COUNTER["n"] += 1
    ds = bg_state.Dataset(
        f"p5ds{_COUNTER['n']}", "p0-5", X=X, genes=genes, samples=samples,
        counts=X.copy(), sample_info=info, source_path=None, source_format="matrix_text",
        species="human", value_kind="counts",
    )
    # 直接标记前置步骤，专注测试 registry 的版本/导出逻辑
    ds.steps = ["qc", "filter", "normalize"]
    return ds


def run_de(ds, workdir, params=None):
    return bg_registry.run_step(ds, "de", dict(params or {"method": "welch"}), workdir, _engine())


def export_json(ds, workdir, what, out_path) -> dict:
    res = bg_registry.export_dataset(ds, what, str(out_path), workdir, _engine())
    return res


def expect_raise(fn):
    try:
        fn()
    except EngineError as error:
        return f"{error.code}: {error.message}"
    except Exception as error:  # noqa: BLE001
        return f"{type(error).__name__}: {error}"
    return ""


# ---------------------------------------------------------------------------
# A. 修订版 ID = 元数据 + 原始层 + 协变量
# ---------------------------------------------------------------------------

def test_revision_id():
    print("\n[A] 修订版 ID 组成")
    workdir = os.path.join(OUTDIR, "rev")
    os.makedirs(workdir, exist_ok=True)
    ds = make_dataset()
    fn = getattr(bg_registry, "_revision_id", None)
    BEFORE["A0"] = "缺少 _revision_id" if fn is None else "有 _revision_id"
    check("A1 存在修订版 ID（16 位 hex）",
          callable(fn) and isinstance(fn(ds), str) and len(fn(ds)) == 16,
          f"{fn(ds) if callable(fn) else None}")
    if not callable(fn):
        for name in ("A2 修订版 ID 可复现（无改动时相同）",
                     "A3 改分组标签 → 修订版 ID 变化",
                     "A4 改协变量列 → 修订版 ID 变化",
                     "A5 改原始 counts 层 → 修订版 ID 变化",
                     "A6 改矩阵 log 层 → 修订版 ID 变化"):
            check(name, False, "缺少 _revision_id")
        return
    first = fn(ds)
    check("A2 修订版 ID 可复现（无改动时相同）", fn(ds) == first, f"{first}")
    ds.sample_info["group"] = ["normal" if v == "tumor" else "tumor" for v in ds.sample_info["group"]]
    after_group = fn(ds)
    check("A3 改分组标签 → 修订版 ID 变化", after_group != first, f"{first} → {after_group}")
    ds.sample_info["group"] = ["tumor"] * 6 + ["normal"] * 6
    ds.sample_info["age"] = np.linspace(40, 70, len(ds.sample_info))   # 新增协变量
    after_covar = fn(ds)
    check("A4 改协变量列 → 修订版 ID 变化", after_covar != first, f"{first} → {after_covar}")
    raw_before = fn(ds)
    counts = ds.counts.copy()
    counts[0, 0] = counts[0, 0] + 7
    ds.counts = counts
    after_raw = fn(ds)
    check("A5 改原始 counts 层 → 修订版 ID 变化", after_raw != raw_before, f"{raw_before} → {after_raw}")
    ds.counts = counts.copy()
    ds.counts[0, 0] -= 7
    matrix_before = fn(ds)
    X = ds.X.copy()
    X[1, 1] = X[1, 1] + 0.5
    ds.X = X
    after_x = fn(ds)
    check("A6 改矩阵 log 层 → 修订版 ID 变化", after_x != matrix_before, f"{matrix_before} → {after_x}")
    parts = getattr(bg_registry, "_revision_parts", None)
    detail = parts(ds) if callable(parts) else {}
    metadata = (detail or {}).get("metadata") or {}
    BEFORE["A7"] = f"covariates={metadata.get('covariates')}"
    check("A7 修订版 ID 记录协变量列表（含新加的 age）",
          "age" in (metadata.get("covariates") or []), str(metadata.get("covariates")))


# ---------------------------------------------------------------------------
# B. 结果缓存打上版本戳
# ---------------------------------------------------------------------------

def test_stamping():
    print("\n[B] 结果缓存版本戳")
    workdir = os.path.join(OUTDIR, "stamp")
    os.makedirs(workdir, exist_ok=True)
    ds = make_dataset()
    run_de(ds, workdir, {"method": "welch", "covars": []})
    store = ds.results.get("de") or {}
    BEFORE["B1"] = f"store keys={sorted(k for k in store if k.endswith('_id') or 'fingerprint' in k)}"
    check("B1 DE 缓存带 revision_id / result_fingerprint / run_id",
          bool(store.get("revision_id")) and bool(store.get("result_fingerprint")) and bool(store.get("run_id")),
          f"{store.get('run_id')}")
    inputs = store.get("result_inputs") or {}
    check("B2 记录实际纳入样本与分组/矩阵版本",
          isinstance(inputs.get("samples"), list) and len(inputs["samples"]) == 12
          and inputs.get("group_col") == "group" and bool(inputs.get("matrix_sha256")),
          f"n={len(inputs.get('samples') or [])} group={inputs.get('group_col')}")
    history = [h for h in ds.history if h.get("step") == "de" and h.get("ok")]
    entry = history[-1] if history else {}
    check("B3 history 记录 revision_id / run_ids", bool(entry.get("revision_id")) and bool(entry.get("run_ids")),
          f"revision={entry.get('revision_id')} run_ids={entry.get('run_ids')}")
    protocol_path = os.path.join(OUTDIR, "stamp", "protocol.json")
    export_json(ds, workdir, "protocol", protocol_path)
    payload = json.loads(Path(protocol_path).read_text(encoding="utf-8"))
    execution = payload.get("execution_record") or {}
    BEFORE["B4"] = f"protocol.revision_id={payload.get('revision_id')}"
    check("B4 protocol 导出带当前修订版 ID",
          bool(payload.get("revision_id")) and payload.get("revision_id") == bg_registry._revision_id(ds),
          f"{payload.get('revision_id')}")
    csv_path = os.path.join(OUTDIR, "stamp", "de.csv")
    res = export_json(ds, workdir, "de_csv", csv_path)
    check("B5 de_csv 导出返回 revision_id / run_id 且与缓存一致",
          res.get("revision_id") == store.get("revision_id") and res.get("run_id") == store.get("run_id"),
          f"{res.get('revision_id')} {res.get('run_id')}")


# ---------------------------------------------------------------------------
# C. 一致性判定
# ---------------------------------------------------------------------------

def test_consistency():
    print("\n[C] 结果 vs 当前数据一致性")
    workdir = os.path.join(OUTDIR, "cons")
    os.makedirs(workdir, exist_ok=True)

    # C1 三组取二：不应误判 stale
    ds = make_dataset(n_case=6, n_ctrl=6)
    ds.sample_info = pd.concat([
        ds.sample_info,
        pd.DataFrame({"group": ["met"] * 4}, index=[f"M{i}" for i in range(4)]),
    ])
    ds.samples = list(ds.samples) + [f"M{i}" for i in range(4)]
    rng = np.random.default_rng(3)
    ds.X = np.vstack([ds.X, rng.negative_binomial(4.0, 0.5, (4, ds.X.shape[1]))])
    ds.counts = np.vstack([ds.counts, rng.negative_binomial(4.0, 0.5, (4, ds.counts.shape[1]))])
    ds.sample_info = ds.sample_info.reindex(ds.samples)
    ds.steps = ["qc", "filter", "normalize"]
    run_de(ds, workdir, {"method": "welch", "case": "tumor", "control": "normal"})
    status = bg_registry._de_consistency(ds)
    BEFORE["C1"] = f"三组取二 _de_consistency={status}"
    check("C1 三组取二不再被误判 stale",
          status.get("status") == "ok", str(status))

    # C2 改分组标签 → stale
    ds2 = make_dataset(seed=12)
    run_de(ds2, workdir)
    ds2.sample_info["group"] = ["normal"] * 6 + ["tumor"] * 6
    status2 = bg_registry._de_consistency(ds2)
    phrase = "does not match the current sample groups"
    methods_text = bg_registry._methods_markdown(ds2, _engine())
    BEFORE["C2"] = f"分组改动后 _de_consistency={status2}；Methods 短语={'有' if phrase in methods_text else '无'}"
    check("C2 改分组标签 → stale 且 Methods 保留 EXPORT-01p 失效文案",
          status2.get("status") == "stale" and phrase in methods_text,
          f"{status2.get('detail')} | phrase={'Y' if phrase in methods_text else 'N'}")

    # C3 改矩阵（重新标准化/批次校正）→ stale
    ds3 = make_dataset(seed=13)
    run_de(ds3, workdir)
    ds3.X = np.asarray(ds3.X, dtype=np.float64) * 1.01
    status3 = bg_registry._de_consistency(ds3)
    BEFORE["C3"] = f"矩阵改动后 _de_consistency={status3}"
    check("C3 矩阵版本变化 → stale", status3.get("status") == "stale", str(status3.get("detail")))

    # C4 富集结果依赖上游 DE 版本
    ds4 = make_dataset(seed=14)
    run_de(ds4, workdir)
    bg_registry.run_step(ds4, "enrich", {"method": "ora", "db": ["GO_BP"]}, workdir, _engine())
    enrich_store = ds4.results.get("enrich") or {}
    inputs = enrich_store.get("result_inputs") or {}
    BEFORE["C4"] = f"enrich depends_on={inputs.get('depends_on')}"
    check("C4 富集缓存记录上游 DE 的 result_fingerprint",
          isinstance(inputs.get("depends_on"), dict) and bool(inputs["depends_on"].get("de")),
          str(inputs.get("depends_on")))
    enrich_store["result_inputs"]["depends_on"]["de"] = "deadbeefdeadbeef"
    status4 = bg_registry._result_consistency(ds4, "enrich")
    check("C5 上游 DE 版本变了 → 富集 stale", status4.get("status") == "stale", str(status4.get("detail")))

    # C6 生存：改事件列 → stale
    ds6 = make_dataset(seed=15, clinical=True)
    bg_registry.run_step(ds6, "survival", {"group_col": "group", "time_col": "os_time",
                                          "event_col": "os_event"}, workdir, _engine())
    ds6.sample_info.loc["T0", "os_event"] = 0 if int(ds6.sample_info.loc["T0", "os_event"]) == 1 else 1
    status6 = bg_registry._result_consistency(ds6, "survival")
    BEFORE["C6"] = f"改事件列后 survival consistency={status6}"
    check("C6 事件列改动 → 生存结果 stale",
          status6.get("status") == "stale" and ("事件" in str(status6.get("detail") or "")
                                               or "event" in str(status6.get("detail") or "").lower()),
          str(status6.get("detail")))


# ---------------------------------------------------------------------------
# D. 导出闸门 + 旧 run_id 检索
# ---------------------------------------------------------------------------

def test_export_gate():
    print("\n[D] 导出闸门与旧 run_id")
    workdir = os.path.join(OUTDIR, "gate")
    os.makedirs(workdir, exist_ok=True)
    ds = make_dataset(seed=21)
    run_de(ds, workdir, {"method": "welch", "case": "tumor", "control": "normal"})
    run1 = (ds.results.get("de") or {}).get("run_id")
    csv1 = os.path.join(OUTDIR, "gate", "run1.csv")
    res1 = export_json(ds, workdir, "de_csv", csv1)
    check("D1 当前版本 de_csv 正常导出", os.path.exists(csv1) and bool(res1.get("run_id")),
          f"run_id={res1.get('run_id')}")
    frame1 = pd.read_csv(csv1)
    lfc1 = float(frame1.set_index("gene").loc["A2M", "logFC"]) if "A2M" in set(frame1["gene"]) else None
    if lfc1 is None:
        lfc1 = float(frame1["logFC"].iloc[0])

    # 旧结果在改分组后必须被拒
    ds.sample_info["group"] = ["normal"] * 6 + ["tumor"] * 6
    os.makedirs(os.path.join(OUTDIR, "gate", "stale"), exist_ok=True)
    msg = expect_raise(lambda: export_json(ds, workdir, "de_csv",
                                           os.path.join(OUTDIR, "gate", "stale", "de.csv")))
    BEFORE["D2"] = f"改分组后导出 de_csv: {msg[:110]}"
    check("D2 改分组后旧 DE 缓存不得导出（prerequisite + 提示重跑）",
          msg.startswith("prerequisite") and ("重跑" in msg or "不一致" in msg), msg[:140])

    # 还原分组，重跑 DE（方向相反）→ 形成第二个 run
    ds.sample_info["group"] = ["tumor"] * 6 + ["normal"] * 6
    run_de(ds, workdir, {"method": "welch", "case": "normal", "control": "tumor"})
    run2 = (ds.results.get("de") or {}).get("run_id")
    BEFORE["D3"] = f"run1={run1} run2={run2}"
    check("D3 重跑产生新 run_id（旧结果进历史）",
          bool(run2) and run2 != run1, f"{run1} / {run2}")
    history = bg_registry._result_history(ds)
    run_ids = [item.get("run_id") for item in history]
    check("D4 _result_history 可按旧 run_id 检索",
          run1 in run_ids and run2 in run_ids, str(run_ids))
    csv_old = os.path.join(OUTDIR, "gate", "run1_old.csv")
    res_old = bg_registry.export_dataset(ds, "de_csv", csv_old, workdir, _engine(), run_id=run1)
    frame_old = pd.read_csv(csv_old) if os.path.exists(csv_old) else pd.DataFrame()
    lfc_old = None
    if not frame_old.empty:
        gene = "A2M" if "A2M" in set(frame_old["gene"]) else str(frame_old["gene"].iloc[0])
        lfc_old = float(frame_old.set_index("gene").loc[gene, "logFC"])
    BEFORE["D5"] = f"历史导出 historical={res_old.get('historical')} run_id={res_old.get('run_id')}"
    check("D5 用旧 run_id 导出历史版本（内容与原 run 一致，方向未反转）",
          res_old.get("historical") is True and res_old.get("run_id") == run1
          and abs(float(lfc_old) - float(lfc1)) < 1e-9,
          f"historical={res_old.get('historical')} lfc={lfc_old} vs {lfc1}")
    history_path = os.path.join(OUTDIR, "gate", "result_history.json")
    export_json(ds, workdir, "result_history", history_path)
    payload = json.loads(Path(history_path).read_text(encoding="utf-8"))
    entries = payload.get("history") if isinstance(payload, dict) else payload
    check("D6 result_history 导出包含旧 run_id / revision_id",
          isinstance(entries, list) and any(item.get("run_id") == run1 for item in entries)
          and all(item.get("revision_id") for item in entries),
          f"n={len(entries) if isinstance(entries, list) else '?'}")

    # 出图闸门：分组改动后旧 DE 不能画火山图
    ds.sample_info["group"] = ["normal"] * 6 + ["tumor"] * 6
    os.makedirs(os.path.join(OUTDIR, "gate", "plot"), exist_ok=True)
    msg_plot = expect_raise(lambda: bg_registry.run_plot(ds, "volcano", {}, workdir, _engine()))
    BEFORE["D7"] = f"改分组后画火山图: {msg_plot[:110]}"
    check("D7 改分组后旧 DE 不能出图（出图前校验版本）",
          msg_plot.startswith("prerequisite") and ("重跑" in msg_plot or "不一致" in msg_plot),
          msg_plot[:140])


def main() -> int:
    print("== dsh-geobulk P0-5 结果溯源 ==")
    test_revision_id()
    test_stamping()
    test_consistency()
    test_export_gate()
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
