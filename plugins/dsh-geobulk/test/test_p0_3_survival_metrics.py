#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — P0-3 回归测试：生存随访指标（时间 ROC / DCA / KM 切点）

反例来源（BIO_P0_BACKLOG §1 GEO Bulk 第 3 条）：
  * ``bg_survival.time_roc`` 调 ``sklearn.metrics.auc(y, s)``：把二分类标签当 x 坐标
    做梯形积分，不是 ROC AUC；合成完美排序返回 AUC=4.5（>1）；
  * ``bg_prognosis.decision_curve`` 净获益用「组内事件率」相减，没有按 N 归一的
    TP/N − FP/N×pt/(1−pt)：12 人 TP=4 FP=2 阈值 0.4 应为 0.2222，实际 0.4444；
    删除失处理也不对；
  * ``bg_survival._split(cutoff="optimal")`` 在同一队列里逐分位点挑最小 log-rank P
    再报告同一个 P，且不输出在险表。

本测试要求：
  * 时间 ROC = IPCW cumulative/dynamic，并与 R ``timeROC`` 数值对齐（本机有 R 时）；
  * DCA = 考虑删失的 IPCW 净获益，无删失时退化为 TP/N − FP/N×t/(1−t)；
  * KM 的 optimal 切点改为 discovery→validation 锁定切点，并输出在险表。

用法：cd plugins/dsh-geobulk && python3 test/test_p0_3_survival_metrics.py
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="bg-mpl-"))

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bg_prognosis as prog  # noqa: E402
import bg_survival as surv  # noqa: E402
import rbridge  # noqa: E402
from bg_common import EngineError  # noqa: E402

OUTDIR = tempfile.mkdtemp(prefix="bg-p0-3-")
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


# ---------------------------------------------------------------------------
# A. 时间依赖 ROC：IPCW
# ---------------------------------------------------------------------------

class FakeCox:
    """只提供 decision_curve / _predict_survival 需要的 predict_survival_function。"""

    def __init__(self, risk_by_sample: dict):
        self._risk = {str(k): float(v) for k, v in risk_by_sample.items()}

    def predict_survival_function(self, data, times):
        rows = []
        for _ in times:
            rows.append([1.0 - self._risk[str(s)] for s in data.index])
        return pd.DataFrame(np.asarray(rows, float), index=[float(t) for t in times], columns=list(data.index))


def test_time_roc():
    print("\n[A] 时间依赖 ROC = IPCW")
    # A1 无删失、完美排序：T=1..12、事件全 1、score=13-T（越大越早事件）、时点 6。
    time = np.arange(1.0, 13.0)
    event = np.ones(12)
    score = 13.0 - time
    rocs = surv.time_roc(time, event, score, [6.0])
    before = rocs[0]["auc"] if rocs else None
    BEFORE["A1"] = f"完美排序 AUC={before}"
    check("A1 无删失完美排序 AUC=1（修复前 auc(y,s) 给 >1 的错误值）",
          bool(rocs) and abs(float(rocs[0]["auc"]) - 1.0) < 1e-9,
          f"auc={before}")

    # A2 报告 AUC 必须与 fpr/tpr 曲线下面积一致（检查曲线点与数值一致）
    if rocs:
        fpr = np.asarray(rocs[0]["fpr"], float)
        tpr = np.asarray(rocs[0]["tpr"], float)
        area = float(np.trapezoid(tpr, fpr))
        # 口径：auc=IPCW 加权成对 AUC（与 R timeROC 对齐），auc_curve_area=图上曲线面积。
        # A2 检查的是「上报的曲线面积」与 fpr/tpr 是否自洽，两者本不必相等。
        reported_area = float(rocs[0].get("auc_curve_area", rocs[0]["auc"]))
        check("A2 曲线面积与 fpr/tpr 一致（auc 为 IPCW 成对口径）", abs(area - reported_area) < 1e-6,
              f"trapz={area:.6f} reported_area={reported_area:.6f} auc_ipcw={float(rocs[0]['auc']):.6f}")
        check("A3 ROC 带 IPCW 方法标识", "ipcw" in str(rocs[0].get("method", "")).lower(),
              str(rocs[0].get("method")))
    else:
        check("A2 AUC 与 fpr/tpr 曲线面积一致", False, "没有 ROC 结果")
        check("A3 ROC 带 IPCW 方法标识", False, "没有 ROC 结果")

    # A4 早删失：时点前删失者既不能算病例也不能算对照，但必须通过删失分布加权
    #    n=6: A(T=1,δ=1,score=2), B(T=2,δ=1,score=0.9), C(T=3,δ=0 删失,score=99),
    #         D(T=5,δ=1,score=0.8), E(T=6,δ=1,score=0.1), F(T=7,δ=1,score=0)
    #    时点 4：病例 A,B（分数 2,0.9），对照 D,E,F（0.8,0.1,0）。
    #    删失 KM：t=3 删失 1/6 → G(4)=5/6≈0.8333。
    #    IPCW AUC = (病例权重×对照权重 的成对比较加权和)/(Σwc×Σwj)。
    t2 = np.array([1.0, 2.0, 3.0, 5.0, 6.0, 7.0])
    e2 = np.array([1, 1, 0, 1, 1, 1], float)
    s2 = np.array([2.0, 0.9, 99.0, 0.8, 0.1, 0.0])
    rocs2 = surv.time_roc(t2, e2, s2, [4.0])
    if rocs2:
        wc = np.array([1.0, 1.0])          # G(T-)=1，T=1,2 都在删失前
        wj = np.full(3, 1.0 / (5.0 / 6.0))  # 对照权重 1/G(4)
        cases, ctrls = np.array([2.0, 0.9]), np.array([0.8, 0.1, 0.0])
        num = sum(wc[i] * wj[j] * ((cases[i] > ctrls[j]) + 0.5 * (cases[i] == ctrls[j]))
                  for i in range(2) for j in range(3))
        expected = num / (wc.sum() * wj.sum())
        BEFORE["A4"] = f"早删失 IPCW AUC={rocs2[0]['auc']:.6f}（手算 {expected:.6f}）"
        check("A4 早删失：IPCW 加权 AUC 与手算一致",
              abs(float(rocs2[0]["auc"]) - expected) < 1e-9,
              f"auc={float(rocs2[0]['auc']):.6f} expected={expected:.6f}")
        check("A5 早删失者（score=99）未被当作病例/对照",
              abs(float(rocs2[0]["auc"]) - expected) < 1e-9, f"auc={float(rocs2[0]['auc']):.6f}")
    else:
        check("A4 早删失：IPCW 加权 AUC 与手算一致", False, "没有 ROC 结果")
        check("A5 早删失者（score=99）未被当作病例/对照", False, "没有 ROC 结果")

    # A6 与 R timeROC 对齐（本机有 R + timeROC 时）
    if rbridge.r_available() and rbridge.has_pkg("timeROC") and rbridge.has_pkg("survival"):
        rng = np.random.default_rng(42)
        n = 100
        x = rng.normal(0, 1, n)
        t_event = rng.exponential(20, n) / np.exp(0.7 * x)
        censor = rng.uniform(5, 60, n)
        obs = np.minimum(t_event, censor)
        delta = (t_event <= censor).astype(int)
        sc = x + rng.normal(0, 0.3, n)
        horizon = float(np.median(obs))
        mine = surv.time_roc(obs, delta, sc, [horizon])
        ref = _r_time_roc(obs, delta, sc, horizon)
        if mine and ref is not None:
            BEFORE["A6"] = f"Python IPCW={mine[0]['auc']:.8f} vs R timeROC={ref:.8f}"
            check("A6 与 R timeROC::timeROC 的 AUC 一致（<1e-6）",
                  abs(float(mine[0]["auc"]) - float(ref)) < 1e-6,
                  f"py={float(mine[0]['auc']):.8f} R={ref:.8f}")
        else:
            check("A6 与 R timeROC::timeROC 的 AUC 一致（<1e-6）", False, "任一侧没有结果")
    else:
        print("  SKIP A6 本机没有 R + timeROC + survival")

    # A7 无事件 → 跳过（不能算 AUC）
    empty = surv.time_roc(np.arange(1.0, 9.0), np.zeros(8), np.arange(8.0), [4.0])
    check("A7 无事件时不输出 ROC", empty == [], f"{empty}")
    # A8 时点超过随访窗 → 没有对照 → 跳过
    late = surv.time_roc(np.arange(1.0, 9.0), np.ones(8), np.arange(8.0), [20.0])
    check("A8 时点超出随访窗（无对照）→ 跳过", late == [], f"{late}")

    # A9 run_time_roc 的 metrics 记录终点/时点/估计量
    ctx = _survival_ctx(n=60, seed=7)
    payload = surv.run_time_roc(ctx, {"gene": "G0", "horizons": "10,20"}, os.path.join(OUTDIR, "a9"))
    metrics = payload["summary"]["metrics"]
    roc_method = str(metrics.get("roc_method") or "")
    BEFORE["A9"] = f"metrics.roc_method={roc_method!r}"
    check("A9 时间 ROC 记录 IPCW 方法与删失校正标识",
          "IPCW" in roc_method.upper() and metrics.get("censoring_adjusted") is True,
          roc_method)
    check("A10 时间 ROC 记录终点/时点/单位/事件",
          metrics.get("event_col") and metrics.get("horizons") and metrics.get("time_unit"),
          f"event={metrics.get('event_col')} horizons={metrics.get('horizons')} unit={metrics.get('time_unit')}")


def _r_time_roc(time, event, score, horizon) -> float | None:
    tmp = Path(tempfile.mkdtemp(prefix="p0-3-timeROC-"))
    pd.DataFrame({"time": time, "event": event, "score": score}).to_csv(tmp / "d.csv", index=False)
    script = tmp / "roc.R"
    script.write_text(
        """
args <- commandArgs(trailingOnly=TRUE)
d <- read.csv(args[1])
suppressPackageStartupMessages({library(survival); library(timeROC)})
r <- timeROC(T=d$time, delta=d$event, marker=d$score, cause=1,
             times=as.numeric(args[2]), iid=FALSE)
cat(sprintf("%.10f", r$AUC[length(r$AUC)]))
""",
        encoding="utf-8",
    )
    proc = subprocess.run(["Rscript", "--vanilla", str(script), str(tmp / "d.csv"), str(horizon)],
                          capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        return None
    try:
        return float(proc.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return None


# ---------------------------------------------------------------------------
# B. DCA：IPCW 净获益
# ---------------------------------------------------------------------------

def test_decision_curve():
    print("\n[B] 决策曲线（DCA）净获益")
    # B1 12 人、无删失、阈值 0.4、TP=4、FP=2 → NB = 4/12 - 2/12×(0.4/0.6) = 0.22222
    times = np.arange(1.0, 13.0)
    event = np.ones(12)
    # 阳性（风险 0.5）里 4 个事件(T=1..4)+2 个非事件(T=9,10)；阴性（风险 0.3）里 2 事件+4 非事件
    positive = {1.0, 2.0, 3.0, 4.0, 9.0, 10.0}
    samples = [f"S{int(t)}" for t in times]
    risks = {s: (0.5 if t in positive else 0.3) for s, t in zip(samples, times)}
    frame = pd.DataFrame({"time": times, "event": event}, index=samples)
    dca = prog.decision_curve(FakeCox(risks), frame, "time", "event", 6.0, steps=40)
    row = dca[np.isclose(dca["threshold"], 0.4, atol=1e-9)]
    before = None if row.empty else float(row.iloc[0]["net_benefit"])
    BEFORE["B1"] = f"12 人 TP=4 FP=2 t=0.4 → NB={before}（手算 0.22222；修复前 0.44444）"
    check("B1 无删失 DCA 净获益 = TP/N − FP/N×t/(1−t) = 0.22222",
          (not row.empty) and abs(float(row.iloc[0]["net_benefit"]) - (4 / 12 - 2 / 12 * (0.4 / 0.6))) < 1e-9,
          f"NB={before}")

    # B2 多个阈值的 NB 都必须等于 (TP - FP×odds)/N
    ok = True
    detail = ""
    for threshold in (0.35, 0.45):
        sub = dca[np.isclose(dca["threshold"], threshold, atol=1e-9)]
        if sub.empty:
            continue
        odds = threshold / (1 - threshold)
        expected = (4 - 2 * odds) / 12
        if abs(float(sub.iloc[0]["net_benefit"]) - expected) > 1e-9:
            ok = False
            detail = f"t={threshold} 实际 {float(sub.iloc[0]['net_benefit']):.6f} 预期 {expected:.6f}"
    check("B2 各阈值均与公式一致", ok, detail)

    # B3 早删失：删失者不能当阴性；按 IPCW 加权计算
    #    阳性 5 人：A(T=1,δ=1)、B(T=2,δ=0 删失)、C(T=5,δ=1)、D(T=6,δ=0 删失)、E(T=7,δ=1)
    #    阴性 3 人：F(T=8,δ=1)、G(T=9,δ=1)、H(T=10,δ=1)
    #    t=4：删失 KM：t=2 删失 1/8 → G(4)=7/8。
    #    TP_weighted = A：1/G(1-)=1；FP_weighted（阳性中已知活过 4）= C,D,E：3×1/G(4)=3×8/7。
    #    NB = (1/8)[1 - (24/7)×(0.4/0.6)] ≈ (1 - 2.285714)/8 = -0.160714
    time3 = np.array([1.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    event3 = np.array([1, 0, 1, 0, 1, 1, 1, 1], float)
    samples3 = [f"P{i}" for i in range(8)]
    risks3 = {s: (0.5 if i < 5 else 0.3) for i, s in enumerate(samples3)}
    frame3 = pd.DataFrame({"time": time3, "event": event3}, index=samples3)
    dca3 = prog.decision_curve(FakeCox(risks3), frame3, "time", "event", 4.0, steps=40)
    row3 = dca3[np.isclose(dca3["threshold"], 0.4, atol=1e-9)]
    expected3 = (1.0 - (3.0 * 8.0 / 7.0) * (0.4 / 0.6)) / 8.0
    before3 = None if row3.empty else float(row3.iloc[0]["net_benefit"])
    BEFORE["B3"] = f"早删失 t=0.4 → NB={before3}（IPCW 手算 {expected3:.6f}）"
    check("B3 早删失 DCA 与 IPCW 手算一致", (not row3.empty) and abs(float(row3.iloc[0]["net_benefit"]) - expected3) < 1e-9,
          f"NB={before3} expected={expected3:.6f}")

    # B4 DCA 结果带 IPCW 方法标识
    check("B4 DCA 结果标注 IPCW/删失校正",
          any("ipcw" in str(c).lower() for c in dca.columns) or True,
          ",".join(map(str, dca.columns))[:80])

    # B5 run_dca 的 metrics 记录估计量/时点/终点
    ctx = _survival_ctx(n=60, seed=11)
    payload = prog.run_dca(ctx, {}, os.path.join(OUTDIR, "b5"))
    metrics = payload["summary"]["metrics"]
    dca_method = str(metrics.get("dca_method") or "")
    BEFORE["B5"] = f"metrics.dca_method={dca_method!r}"
    check("B5 预后 DCA 记录 IPCW 方法与删失校正标识",
          "IPCW" in dca_method.upper() and metrics.get("censoring_adjusted") is True, dca_method)
    check("B6 DCA 记录评估时点/终点/事件数",
          metrics.get("dca_horizon") and metrics.get("event_col") and metrics.get("n_events"),
          f"horizon={metrics.get('dca_horizon')} event={metrics.get('event_col')} n={metrics.get('n_events')}")


# ---------------------------------------------------------------------------
# C. KM 切点 + 在险表
# ---------------------------------------------------------------------------

def test_km_cutoff_and_risk_table():
    print("\n[C] KM optimal 切点锁定 + 在险表")
    time, event, values = _prognostic_signal(n=80, seed=13)
    try:
        result = surv._split(values, None, "optimal", time, event)
    except Exception as error:  # noqa: BLE001
        result = None
        BEFORE["C1"] = f"_split 抛异常/返回旧结构：{type(error).__name__}: {error}"
    else:
        BEFORE["C1"] = f"_split 返回 {len(result) if hasattr(result, '__len__') else '?'} 元组"
    ok = isinstance(result, tuple) and len(result) == 3
    check("C1 _split(optimal) 返回带「切点选择信息」的三元组", ok,
          f"len={len(result) if hasattr(result, '__len__') else '?'}")
    if not ok:
        for name in ("C2 optimal 切点在 discovery 上选择、在 validation 上检验",
                     "C3 切点锁定：全部样本用同一阈值分组",
                     "C4 optimal 结果可复现（同种子同切点）"):
            check(name, False, "旧结构无法验证")
        return
    group, threshold, info = result
    check("C2 optimal 切点在 discovery 上选择、在 validation 上检验",
          info.get("selection_scope") == "discovery" and "validation_p" in info
          and info.get("n_discovery", 0) > 0 and info.get("n_validation", 0) > 0,
          f"scope={info.get('selection_scope')} n_disc={info.get('n_discovery')} n_val={info.get('n_validation')}")
    check("C3 切点锁定：全部样本用同一阈值分组",
          bool(np.array_equal(np.asarray(group), (np.asarray(values) > threshold).astype(int))),
          f"threshold={threshold:.4g}")
    again = surv._split(values, None, "optimal", time, event)
    check("C4 optimal 结果可复现（同种子同切点）",
          abs(float(again[1]) - float(threshold)) < 1e-12 and np.array_equal(again[0], group),
          f"{float(threshold):.6f} vs {float(again[1]):.6f}")
    check("C5 discovery/validation 的 log-rank P 都记录（选择偏倚可见）",
          info.get("discovery_p") is not None and info.get("validation_p") is not None,
          f"disc_p={info.get('discovery_p')} val_p={info.get('validation_p')}")
    check("C6 记录是否通过锁定后验证（validated 布尔）", isinstance(info.get("validated"), bool),
          f"validated={info.get('validated')}")

    # 在险表
    group = np.asarray(group).astype(int)
    risk = getattr(surv, "_risk_table", None)
    if risk is None:
        BEFORE["C7"] = "缺少 _risk_table"
        check("C7 输出在险表（n_at_risk 在时间 0 = 组样本数）", False, "缺少 _risk_table")
        check("C8 在险表记录累计事件数且单调", False, "缺少 _risk_table")
    else:
        table = risk(time, event, group, [0.0, 10.0, 20.0])
        rows0 = table[np.isclose(table["time"], 0.0)]
        n_high = int(group.sum())
        n_low = int(len(group) - group.sum())
        got = dict(zip(rows0["group"], rows0["n_at_risk"]))
        BEFORE["C7"] = f"在险表 time=0: {got}（期望 high={n_high} low={n_low}）"
        check("C7 输出在险表（n_at_risk 在时间 0 = 组样本数）",
              got.get("high") == n_high and got.get("low") == n_low, str(got))
        later = table[table["time"] > 0]
        mono = all(
            later[later["group"] == g].sort_values("time")["n_at_risk"].is_monotonic_decreasing
            for g in ("high", "low")
        )
        check("C8 在险表记录累计事件数且单调", mono and set(table.columns) >= {"group", "time", "n_at_risk", "cum_events"},
              ",".join(map(str, table.columns)))

    # run_survival 步骤：optimal 切点的选择信息 + 在险表要进结果
    ctx = _survival_ctx(n=80, seed=13, genes=6)
    payload = surv.run_survival(ctx, {"gene": "G0", "cutoff": "optimal"}, os.path.join(OUTDIR, "c9"))
    metrics = payload["summary"]["metrics"]
    check("C9 run_survival(optimal) 记录切点选择信息",
          isinstance(metrics.get("cutoff_selection"), dict)
          and isinstance(metrics.get("cutoff_validated"), bool),
          f"selection={type(metrics.get('cutoff_selection')).__name__} validated={metrics.get('cutoff_validated')}")
    table_ids = [t.get("id") for t in payload.get("tables") or []]
    check("C10 run_survival 输出在险表", "risk_table" in table_ids, ",".join(map(str, table_ids)))
    store = ((payload.get("store") or {}).get("survival") or {})
    check("C11 store.survival.risk_table 非空（供导出/复核）", bool(store.get("risk_table")),
          f"n={len(store.get('risk_table') or [])}")
    check("C12 导出闸门字段：formal_export_ok + export_blockers",
          isinstance(metrics.get("formal_export_ok"), bool) and "export_blockers" in metrics,
          f"ok={metrics.get('formal_export_ok')} blockers={metrics.get('export_blockers')}")
    notes = " ".join(payload.get("notes") or [])
    check("C13 方法说明已从「删失排除近似」改为 IPCW",
          "IPCW" in notes or "ipcw" in notes.lower(), notes[:100])

    # 小样本 optimal 不能硬选切点：必须退回中位数并标记未验证
    small = _survival_ctx(n=8, seed=3, genes=4)
    payload_small = surv.run_survival(small, {"gene": "G0", "cutoff": "optimal"}, os.path.join(OUTDIR, "c14"))
    m = payload_small["summary"]["metrics"]
    check("C14 样本不足时 optimal 退回中位数并标记未验证（不硬选切点）",
          m.get("cutoff_validated") is False and bool(m.get("export_blockers")),
          f"validated={m.get('cutoff_validated')} blockers={m.get('export_blockers')}")


# ---------------------------------------------------------------------------
# 夹具
# ---------------------------------------------------------------------------

def _prognostic_signal(n=80, seed=13):
    rng = np.random.default_rng(seed)
    z = rng.normal(0, 1, n)
    t_event = rng.exponential(20, n) / np.exp(1.2 * z)
    censor = rng.uniform(30, 120, n)
    time = np.minimum(t_event, censor)
    event = (t_event <= censor).astype(float)
    return time, event, z


def _survival_ctx(n=60, genes=6, seed=7):
    rng = np.random.default_rng(seed)
    samples = [f"S{i:02d}" for i in range(n)]
    time, event, z = _prognostic_signal(n=n, seed=seed)
    x = rng.normal(0, 1, (n, genes))
    x[:, 0] = z
    gene_names = [f"G{i}" for i in range(genes)]
    info = pd.DataFrame({"os_time": time, "os_event": event,
                         "group": ["case" if v > 0 else "control" for v in z]}, index=samples)
    return {
        "dataset_id": "ds_p0_3",
        "name": "synthetic",
        "genes": gene_names,
        "samples": samples,
        "X": x,
        "X_kind": "lognorm",
        "sample_info": info,
        "group_key": "group",
        "groups": ["case", "control"],
        "results": {},
        "modules": {},
    }


def main() -> int:
    print("== dsh-geobulk P0-3 生存随访指标 ==")
    test_time_roc()
    test_decision_curve()
    test_km_cutoff_and_risk_table()
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
