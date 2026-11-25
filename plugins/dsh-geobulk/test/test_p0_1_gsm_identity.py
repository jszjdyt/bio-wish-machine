#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — P0-1 回归测试：GSM 元数据错配 / 数据集身份

反例来源（BIO_P0_BACKLOG §1 GEO Bulk 第 1 条）：
  * lib/bg_io.py 解析 GEO series matrix 时按「!Sample_* 行的列位置」复制元数据，
    矩阵列与 !Sample_ 字段同一组 GSM 但顺序相反时，两组标签对调；
  * lib/bg_engine.py / lib/bg_state.py 只按路径复用旧数据对象，同一文件换新样本表
    仍返回旧分组；同一路径文件内容变了也照样复用。

本测试全部离线（不联网），直接构造合成 series matrix。

用法：cd plugins/dsh-geobulk && python3 test/test_p0_1_gsm_identity.py
失败打印 FAIL 并以 1 退出，成功打印 PASS 并以 0 退出。
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

import bg_io  # noqa: E402
from bg_common import EngineError  # noqa: E402

PASSED: list[str] = []
FAILED: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


def write_series(path: Path, *, header: list[str], accessions: list[str], titles: list[str],
                 groups: list[str]) -> None:
    """写一个最小 series matrix：!Sample_ 行按 accessions 顺序，表头按 header 顺序。"""
    lines = [
        '!Series_title\t"synthetic"',
        '!Sample_geo_accession\t' + "\t".join(f'"{a}"' for a in accessions),
        '!Sample_title\t' + "\t".join(f'"{t}"' for t in titles),
        '!Sample_characteristics_ch1\t' + "\t".join(f'"group: {g}"' for g in groups),
        '!Sample_organism_ch1\t' + "\t".join('"Homo sapiens"' for _ in accessions),
        "!series_matrix_table_begin",
        '"ID_REF"\t' + "\t".join(f'"{h}"' for h in header),
    ]
    # 每个 GSM 一个可区分的固定均值：GSM1=1, GSM2=11（基因 TP53/CD8A/EPCAM）
    base = {"GSM1": 1.0, "GSM2": 11.0, "GSM3": 21.0}
    for gene in ("TP53", "CD8A", "EPCAM"):
        cells = []
        for column in header:
            value = base.get(column, 5.0)
            cells.append(f"{value:.1f}")
        lines.append(f'"{gene}"\t' + "\t".join(cells))
    lines.append("!series_matrix_table_end")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_matrix(path: Path, **options) -> dict:
    return bg_io.load_source(str(path), workdir=str(path.parent), options=dict(options))


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="bg-p0-1-"))
    print("== P0-1 GSM 按 accession 一对一连接 ==")

    # ---- 反例 1：表头顺序与 !Sample_ 顺序相反 -------------------------------
    flipped = tmp / "GSE001_series_matrix.txt"
    write_series(flipped, header=["GSM2", "GSM1"], accessions=["GSM1", "GSM2"],
                 titles=["ctrl_A", "tumor_B"], groups=["ctrl", "tumor"])
    loaded = load_matrix(flipped)
    info = loaded["sample_info"]
    title_gsm2 = str(info.loc["GSM2", "title"]) if "title" in info.columns else "<no title col>"
    group_gsm2 = str(info.loc["GSM2", "group"]) if "group" in info.columns else "<no group col>"
    check("1a 乱序表头：GSM2 的 title 来自 GSM2（不是按位置抄 GSM1）",
          title_gsm2 == "tumor_B", f"title={title_gsm2}")
    check("1b 乱序表头：GSM2 的分组 = tumor", group_gsm2 == "tumor", f"group={group_gsm2}")
    check("1c 乱序表头：GSM1 的分组 = ctrl",
          str(info.loc["GSM1", "group"]) == "ctrl", str(info.loc["GSM1", "group"]))
    # 矩阵列本身不能被重排（列值仍属该 GSM）
    X = np.asarray(loaded["X"], dtype=float)
    gsm2_col = list(loaded["samples"]).index("GSM2")
    gene_row = list(loaded["genes"]).index("TP53")
    check("1d 矩阵列值仍属 GSM2（均值 11）", abs(X[gsm2_col, gene_row] - 11.0) < 1e-6,
          f"value={X[gsm2_col, gene_row]}")
    joined_notes = " ".join(str(n) for n in loaded.get("notes") or [])
    check("1e 记录按 accession 一对一连接", "accession" in joined_notes.lower() or "GSM" in joined_notes,
          joined_notes[:80])

    # ---- 反例 2：重复 GSM（!Sample_ 里同一 accession 出现两次）---------------
    dup = tmp / "GSE002_series_matrix.txt"
    write_series(dup, header=["GSM1", "GSM1"], accessions=["GSM1", "GSM1"],
                 titles=["a", "b"], groups=["ctrl", "tumor"])
    dup_warn: list[str] = []
    try:
        dup_loaded = load_matrix(dup)
        dup_warn = [str(w) for w in dup_loaded.get("warnings") or []]
        dup_ok = True
    except EngineError as error:
        dup_warn = [str(error)]
        dup_ok = True
    dup_text = " ".join(dup_warn)
    check("2a 重复 GSM 明确列出（警告或报错，不静默）",
          dup_ok and ("重复" in dup_text or "duplicate" in dup_text.lower()),
          dup_text[:120])
    check("2b 重复 GSM 的警告点名 GSM1", "GSM1" in dup_text, dup_text[:120])

    # ---- 反例 3：缺失 accession（表头有 GSM2，!Sample_ 只有 GSM1）-----------
    missing = tmp / "GSE003_series_matrix.txt"
    write_series(missing, header=["GSM1", "GSM2"], accessions=["GSM1"],
                 titles=["only_one"], groups=["ctrl"])
    miss_warn: list[str] = []
    try:
        miss_loaded = load_matrix(missing)
        miss_warn = [str(w) for w in miss_loaded.get("warnings") or []]
    except EngineError as error:
        miss_warn = [str(error)]
    miss_text = " ".join(miss_warn)
    check("3a 缺失 GSM 明确列出", "GSM2" in miss_text, miss_text[:140])

    # ---- 反例 4：显式映射（表头不是 GSM 时可一对一指定）---------------------
    named = tmp / "GSE004_series_matrix.txt"
    write_series(named, header=["S1", "S2"], accessions=["GSM1", "GSM2"],
                 titles=["ctrl_A", "tumor_B"], groups=["ctrl", "tumor"])
    mapped = load_matrix(named, sample_map={"S1": "GSM2", "S2": "GSM1"})
    minfo = mapped["sample_info"]
    check("4a 显式映射：S1 取 GSM2 的 title",
          str(minfo.loc["S1", "title"]) == "tumor_B", str(minfo.loc["S1", "title"]))
    check("4b 显式映射：S2 取 GSM1 的 title",
          str(minfo.loc["S2", "title"]) == "ctrl_A", str(minfo.loc["S2", "title"]))

    # ---- 反例 5：同一路径换样本表 / 换字节，不得复用旧数据集 ----------------
    print("== P0-1 数据集身份 = 文件字节 + 导入选项 + 样本表版本 ==")
    import bg_engine

    work = tmp / "work"
    work.mkdir()
    matrix = tmp / "GSE005_series_matrix.txt"
    write_series(matrix, header=["GSM1", "GSM2"], accessions=["GSM1", "GSM2"],
                 titles=["ctrl_A", "tumor_B"], groups=["ctrl", "tumor"])
    sheet_a = tmp / "sheet_a.csv"
    sheet_a.write_text("sample,group\nGSM1,ctrl\nGSM2,tumor\n", encoding="utf-8")
    sheet_b = tmp / "sheet_b.csv"
    sheet_b.write_text("sample,group\nGSM1,ctrl\nGSM2,ctrl\n", encoding="utf-8")

    worker = bg_engine.Worker(str(work))
    worker.load_state()
    first = worker.cmd_open({"path": str(matrix), "sample_sheet": str(sheet_a)})
    ds1 = first["dataset"]["id"]
    group_1 = first["dataset"].get("groups")
    second = worker.cmd_open({"path": str(matrix), "sample_sheet": str(sheet_a)})
    check("5a 同文件同脚本重复 open 复用同一数据集",
          bool(second.get("reused")) and second["dataset"]["id"] == ds1,
          f"reused={second.get('reused')} id={second['dataset']['id']}")

    third = worker.cmd_open({"path": str(matrix), "sample_sheet": str(sheet_b)})
    ds3 = third["dataset"]["id"]
    check("5b 同一路径换样本表 → 新数据集（不复用旧分组）",
          (not third.get("reused")) and ds3 != ds1,
          f"reused={third.get('reused')} id={ds3}")
    ds3_obj = worker.state.get(ds3)
    labels_b = [str(v) for v in ds3_obj.sample_info["group"].tolist()]
    check("5c 新数据集的分组来自新样本表", labels_b == ["ctrl", "ctrl"], str(labels_b))
    ds1_obj = worker.state.get(ds1)
    check("5d 旧数据集仍保留旧分组", [str(v) for v in ds1_obj.sample_info["group"].tolist()] == group_1,
          str(group_1))

    # 同路径、同选项，但文件内容变化 → 必须重新导入
    write_series(matrix, header=["GSM1", "GSM2"], accessions=["GSM1", "GSM2"],
                 titles=["ctrl_A", "tumor_B"], groups=["ctrl", "tumor"])
    text = matrix.read_text(encoding="utf-8").replace('"TP53"\t1.0\t11.0', '"TP53"\t3.0\t13.0')
    matrix.write_text(text, encoding="utf-8")
    fourth = worker.cmd_open({"path": str(matrix), "sample_sheet": str(sheet_b)})
    check("5e 同路径但字节改变 → 新数据集", (not fourth.get("reused")) and fourth["dataset"]["id"] != ds3,
          f"reused={fourth.get('reused')} id={fourth['dataset']['id']}")
    ds4_obj = worker.state.get(fourth["dataset"]["id"])
    tp53 = list(ds4_obj.genes).index("TP53")
    value = float(ds4_obj.X[0, tp53])
    check("5f 新数据集的矩阵是改后的字节", abs(value - 3.0) < 1e-6, f"TP53={value}")

    # 同一目录、不同 member → 身份不同
    archive_dir = tmp / "multi"
    archive_dir.mkdir()
    write_series(archive_dir / "a_series_matrix.txt", header=["GSM1", "GSM2"],
                 accessions=["GSM1", "GSM2"], titles=["a", "b"], groups=["ctrl", "tumor"])
    write_series(archive_dir / "b_series_matrix.txt", header=["GSM3", "GSM1"],
                 accessions=["GSM1", "GSM3"], titles=["c", "d"], groups=["x", "y"])
    fa = bg_io.source_fingerprint(str(archive_dir / "a_series_matrix.txt"))
    fb = bg_io.source_fingerprint(str(archive_dir / "b_series_matrix.txt"))
    check("5g 不同文件 → 不同身份", fa["fingerprint"] != fb["fingerprint"],
          f"{fa['fingerprint'][:12]} vs {fb['fingerprint'][:12]}")
    f_opt = bg_io.source_fingerprint(str(matrix), options={"transpose": True})
    f_opt2 = bg_io.source_fingerprint(str(matrix), options={"transpose": False})
    check("5h 导入选项进入身份", f_opt["fingerprint"] != f_opt2["fingerprint"], "")
    f_sheet = bg_io.source_fingerprint(str(matrix), options={"sample_sheet": str(sheet_a)})
    f_sheet2 = bg_io.source_fingerprint(str(matrix), options={"sample_sheet": str(sheet_b)})
    check("5i 样本表版本进入身份", f_sheet["fingerprint"] != f_sheet2["fingerprint"], "")

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
