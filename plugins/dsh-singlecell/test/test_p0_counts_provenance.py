#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P0-1 回归：原始层来源被洗白（counts provenance）。

把任务单里的最小反例逐条固化成断言：

  1. 纯 counts H5AD（只有 .X、整数）的 .X 不能被当成「未验证」挡掉；
  2. X 只有 log 标准化、没有 counts 层 → 必须是 processed 模式（counts 不可用），
     依赖原始计数的 normalize 步骤要拒绝执行，而不是拿 log 值当 counts；
  3. counts 层的整数校验必须分块全量扫描，不能只看前 2000 个非零值
     （否则「整数前缀 + 小数尾部」的假 counts 会被判为已验证）；
  4. X 是整数 counts、raw 是 log 值（文件写反）时，要放弃 raw、用 X；
  5. 子集（subset）必须继承父数据的真实尺度/验证状态，不能默认已验证；
  6. 导出 h5ad：未验证的占位矩阵不得写成 layers['counts']；已验证 counts
     必须 round-trip 回来仍是已验证；未验证数据 round-trip 后仍是 processed；
  7. basename 拼 barcode 不能产生重复细胞键：两个同名不同目录的文件合并后
     细胞索引必须唯一，键是「稳定样本 ID×barcode」；部分样本没有 counts 时
     合并不许把表达矩阵洗成已验证 counts。

运行：
    cd plugins/dsh-singlecell && python3 test/test_p0_counts_provenance.py
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
PLUGIN_DIR = HERE.parent
LIB = PLUGIN_DIR / "lib"
TMP_ROOT = Path(tempfile.mkdtemp(prefix="sc-p0-provenance-"))
WORK = TMP_ROOT / "work"
WORK.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(TMP_ROOT / "mpl"))
os.environ["DS_SC_WORKDIR"] = str(WORK)
sys.dont_write_bytecode = True

PASSED: list[str] = []
FAILED: list[str] = []


def check(name, cond, detail=""):
    if cond:
        PASSED.append(name)
        print("  PASS %s%s" % (name, ("  [" + str(detail) + "]") if detail else ""))
        return True
    FAILED.append(name)
    print("  FAIL %s%s" % (name, ("  [" + str(detail) + "]") if detail else ""))
    return False


def _load(module_name, filename):
    spec = importlib.util.spec_from_file_location(module_name, str(LIB / filename))
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


SC_IO = _load("_p0_sc_io", "sc_io.py")
SC_META = _load("_p0_sc_metadata", "sc_metadata.py")
ENGINE = _load("_p0_sc_engine", "sc_engine.py")

try:
    import anndata as ad
except Exception:  # pragma: no cover
    ad = None


def _fresh():
    ENGINE.DATASETS.clear()
    try:
        ENGINE._EXPM1_CACHE.clear()
    except Exception:
        pass


def cmd(command, args=None):
    return ENGINE.COMMANDS[command](dict(args or {}))


# 所有依赖「原始 counts」的步骤都要先过 qc→filter；测试里用最宽松的过滤参数，
# 以免小矩阵被默认阈值（min_genes=200）整体过滤掉——那会让断言误报成「步骤顺序」错。
_PERMISSIVE = {"min_genes": 0, "max_pct_mt": 100.0}


def run_to_normalize(dataset_id):
    """按真实依赖序跑到 normalize，返回 (summary, dataset)。"""
    cmd("run", {"dataset": dataset_id, "step": "qc", "params": dict(_PERMISSIVE)})
    cmd("run", {"dataset": dataset_id, "step": "filter",
                "params": dict(_PERMISSIVE, min_cells=0)})
    return cmd("run", {"dataset": dataset_id, "step": "normalize", "params": {"target_sum": 10000}})


def _counts_semantics_error(error):
    """错误是否在说「counts 语义不可用」，而不是「步骤顺序不对」。"""
    message = str(getattr(error, "message", error) or "")
    return any(token in message for token in ("counts", "原始计数", "计数矩阵"))


def _write_h5ad(path, X, genes, cells=None, layers=None, raw_X=None, raw_genes=None,
                obs=None):
    if ad is None:
        raise RuntimeError("需要 anndata")
    n_obs = np.asarray(X).shape[0] if not sp.issparse(X) else X.shape[0]
    obs_frame = obs if obs is not None else pd.DataFrame(index=["c%d" % i for i in range(n_obs)])
    data = ad.AnnData(X=sp.csr_matrix(np.asarray(X, dtype=np.float32)) if not sp.issparse(X) else sp.csr_matrix(X),
                      obs=obs_frame,
                      var=pd.DataFrame(index=list(genes)))
    if cells is not None:
        data.obs_names = list(cells)
    for name, value in (layers or {}).items():
        data.layers[name] = sp.csr_matrix(value) if not sp.issparse(value) else sp.csr_matrix(value)
    if raw_X is not None:
        raw = ad.AnnData(X=sp.csr_matrix(np.asarray(raw_X, dtype=np.float32)),
                         obs=data.obs.copy(), var=pd.DataFrame(index=list(raw_genes or genes)))
        data.raw = raw
    data.write_h5ad(str(path))
    return path


def _dense(matrix):
    if matrix is None:
        return None
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def _counts_meta(payload):
    meta = payload.get("meta") or {}
    return {
        "verified": meta.get("counts_verified"),
        "source": meta.get("counts_source"),
        "kind": meta.get("matrix_kind"),
    }


# ---------------------------------------------------------------------------
# 1. 纯 counts H5AD 的 .X 不能被挡
# ---------------------------------------------------------------------------


def test_x_only_counts():
    print("\n[1] X-only 整数 counts H5AD")
    genes = ["G1", "G2", "G3"]
    X = np.array([[5, 0, 3], [0, 7, 1], [2, 2, 0], [4, 0, 0]], dtype=np.float32)
    path = _write_h5ad(TMP_ROOT / "x-only-counts.h5ad", X, genes)
    payload = SC_IO.load_source(str(path), workdir=str(TMP_ROOT / "w1"))
    counts = _dense(payload.get("counts"))
    info = _counts_meta(payload)
    check("1a 纯 counts 的 .X 被识别为可用 counts",
          counts is not None and np.array_equal(counts, X), str(counts))
    check("1b counts_verified=True 且来源写明来自 X",
          info["verified"] is True and "X" in str(info["source"]), str(info))
    check("1c matrix_kind=counts", info["kind"] == "counts", str(info["kind"]))

    _fresh()
    opened = cmd("open", {"path": str(path)})
    dataset = ENGINE.DATASETS[opened["dataset"]["id"]]
    check("1d 引擎 open 后 counts_verified=True", dataset.counts_verified is True,
          str(dataset.counts_verified))
    try:
        run_to_normalize(dataset.id)
        normalized = ENGINE.DATASETS[opened["dataset"]["id"]].X
        ok = (normalized.max() < 20) if sp.issparse(normalized) else (np.asarray(normalized).max() < 20)
        check("1e 纯 counts 数据可以正常 normalize", bool(ok), "max=%s" % float(normalized.max()))
    except Exception as error:  # noqa: BLE001
        check("1e 纯 counts 数据可以正常 normalize", False,
              "%s: %s" % (type(error).__name__, getattr(error, "message", error)))


# ---------------------------------------------------------------------------
# 2. X 只有 log 标准化 → processed 模式
# ---------------------------------------------------------------------------


def test_lognorm_only_is_processed():
    print("\n[2] X 只有 log 标准化（无 counts 层）")
    genes = ["G1", "G2", "G3"]
    counts = np.array([[5, 1, 3], [0, 7, 1], [2, 2, 0], [4, 1, 2]], dtype=np.float64)
    lognorm = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4)
    path = _write_h5ad(TMP_ROOT / "log-only.h5ad", lognorm, genes)
    payload = SC_IO.load_source(str(path), workdir=str(TMP_ROOT / "w2"))
    info = _counts_meta(payload)
    check("2a 没有 counts 层时 counts=None", payload.get("counts") is None, str(payload.get("counts")))
    check("2b counts_verified=False", info["verified"] is False, str(info))
    check("2c matrix_kind 不是 counts", info["kind"] != "counts", str(info["kind"]))

    _fresh()
    opened = cmd("open", {"path": str(path)})
    dataset = ENGINE.DATASETS[opened["dataset"]["id"]]
    check("2d 引擎进入 processed 模式（counts_verified=False）", dataset.counts_verified is False,
          str(dataset.counts_verified))
    raised = None
    try:
        run_to_normalize(dataset.id)
    except Exception as error:  # noqa: BLE001
        raised = error
    check("2e 拒绝把 log 值当 counts 做 normalize",
          raised is not None
          and str(getattr(raised, "code", "")) == "prerequisite"
          and _counts_semantics_error(raised),
          "%s: %s" % (type(raised).__name__, getattr(raised, "message", raised)))
    check("2f 警告里说明 counts 语义未验证",
          any("counts" in str(w) and ("未验证" in str(w) or "没有" in str(w))
              for w in dataset.warnings),
          str((dataset.warnings or [""])[-1])[:90])


# ---------------------------------------------------------------------------
# 3. 假 counts：整数前缀 + 小数尾部 → 必须全量扫描识别
# ---------------------------------------------------------------------------


def test_fake_counts_full_scan():
    print("\n[3] counts 层必须分块全量校验（不能只看前 2000 个值）")
    n_genes = 2600
    genes = ["G%04d" % i for i in range(n_genes)]
    row = np.ones((2, n_genes), dtype=np.float64)
    row[:, 2000:] = 0.375          # 前 2000 个非零是整数，尾部是小数
    path = _write_h5ad(TMP_ROOT / "fake-counts.h5ad", np.zeros((2, n_genes)),
                       genes, layers={"counts": row})
    payload = SC_IO.load_source(str(path), workdir=str(TMP_ROOT / "w3"))
    info = _counts_meta(payload)
    check("3a 假 counts 被判为未验证（counts=None）", payload.get("counts") is None,
          str(info))
    check("3b counts_verified=False", info["verified"] is False, str(info))
    check("3c 给出未通过校验的说明",
          any(("整数" in str(w) or "非负" in str(w) or "有限" in str(w)) for w in payload.get("warnings") or []),
          str((payload.get("warnings") or [""])[:2]))


# ---------------------------------------------------------------------------
# 4. X 是 counts、raw 是 log（写反）→ 放弃 raw
# ---------------------------------------------------------------------------


def test_raw_log_with_counts_X():
    print("\n[4] X=整数 counts，raw=log 值（文件写反）")
    genes = ["G1", "G2", "G3"]
    X = np.array([[5, 0, 3], [0, 7, 1], [2, 2, 0]], dtype=np.float32)
    raw_log = np.log1p(X * 0.5 + 0.25)
    path = _write_h5ad(TMP_ROOT / "raw-log.h5ad", X, genes, raw_X=raw_log, raw_genes=genes)
    payload = SC_IO.load_source(str(path), workdir=str(TMP_ROOT / "w4"))
    counts = _dense(payload.get("counts"))
    info = _counts_meta(payload)
    check("4a 放弃 log 值的 raw，改用可验证的 X",
          counts is not None and np.array_equal(counts, X), str(counts))
    check("4b 来源标明来自 X 而不是 raw",
          info["verified"] is True and "X" in str(info["source"]) and "raw" not in str(info["source"]).lower(),
          str(info))
    check("4c 记录 raw 被放弃的原因",
          any("raw" in str(w) for w in payload.get("warnings") or []),
          str((payload.get("warnings") or [""])[:2]))


# ---------------------------------------------------------------------------
# 5. 子集继承真实尺度
# ---------------------------------------------------------------------------


def test_subset_inherits_scale():
    print("\n[5] 子集必须继承父数据尺度")
    genes = ["G1", "G2", "G3"]
    counts = np.array([[5, 1, 3], [0, 7, 1], [2, 2, 0], [4, 1, 2], [1, 0, 6], [3, 3, 3]],
                      dtype=np.float64)
    lognorm = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4)
    obs = pd.DataFrame({"clusters": ["0", "0", "0", "1", "1", "1"]},
                       index=["c%d" % i for i in range(6)])
    path = _write_h5ad(TMP_ROOT / "subset-src.h5ad", lognorm, genes, obs=obs)
    _fresh()
    opened = cmd("open", {"path": str(path)})
    parent_id = opened["dataset"]["id"]
    sub = cmd("subset", {"dataset": parent_id, "clusters": ["1"], "column": "clusters"})
    child = ENGINE.DATASETS[sub["dataset"]["id"]]
    check("5a 未验证父数据 → 子集仍是未验证",
          child.counts_verified is False and child.counts_source == ENGINE.DATASETS[parent_id].counts_source,
          "%s / %s" % (child.counts_verified, child.counts_source))
    check("5b 子集保留 matrix_kind", child.matrix_kind == ENGINE.DATASETS[parent_id].matrix_kind,
          "%s vs %s" % (child.matrix_kind, ENGINE.DATASETS[parent_id].matrix_kind))

    X = np.array([[5, 1, 3], [0, 7, 1], [2, 2, 0], [4, 1, 2], [1, 0, 6], [3, 3, 3]],
                 dtype=np.float32)
    path2 = _write_h5ad(TMP_ROOT / "subset-src-counts.h5ad", X, genes, obs=obs)
    opened2 = cmd("open", {"path": str(path2)})
    parent2 = ENGINE.DATASETS[opened2["dataset"]["id"]]
    sub2 = cmd("subset", {"dataset": parent2.id, "clusters": ["0"], "column": "clusters"})
    child2 = ENGINE.DATASETS[sub2["dataset"]["id"]]
    check("5c 已验证父数据 → 子集保持已验证",
          child2.counts_verified is True and child2.counts_source == parent2.counts_source,
          "%s / %s" % (child2.counts_verified, child2.counts_source))


# ---------------------------------------------------------------------------
# 6. 导出 h5ad 的 layers/counts
# ---------------------------------------------------------------------------


def _read_back(path):
    data = ad.read_h5ad(str(path))
    layers = list(getattr(data, "layers", {}).keys())
    uns = {}
    raw_uns = data.uns.get("dsh_singlecell")
    if raw_uns is not None:
        try:
            import json
            uns = json.loads(str(raw_uns))
        except Exception:
            uns = {}
    return data, layers, uns


def test_export_h5ad_layers():
    print("\n[6] 导出 h5ad 不得把归一化占位写成 layers/counts")
    if ad is None:
        check("6 需要 anndata", False, "缺少 anndata")
        return
    genes = ["G1", "G2", "G3"]
    counts = np.array([[5, 1, 3], [0, 7, 1], [2, 2, 0], [4, 1, 2]], dtype=np.float64)
    lognorm = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4)
    unverified = _write_h5ad(TMP_ROOT / "exp-unverified.h5ad", lognorm, genes)
    _fresh()
    opened = cmd("open", {"path": str(unverified)})
    dataset_id = opened["dataset"]["id"]
    out = TMP_ROOT / "export-unverified.h5ad"
    cmd("export", {"dataset": dataset_id, "what": "h5ad", "path": str(out)})
    data, layers, uns = _read_back(out)
    check("6a 未验证数据导出后没有 layers['counts']", "counts" not in layers, str(layers))
    check("6b 导出记录 counts_verified=False",
          uns.get("counts_verified") is False, str(uns.get("counts_verified")))
    _fresh()
    reopened = cmd("open", {"path": str(out)})
    roundtrip = ENGINE.DATASETS[reopened["dataset"]["id"]]
    check("6c round-trip 后仍是 processed 模式", roundtrip.counts_verified is False,
          str(roundtrip.counts_verified))

    verified = _write_h5ad(TMP_ROOT / "exp-verified.h5ad", counts.astype(np.float32), genes,
                           layers={"counts": counts.astype(np.float32)})
    _fresh()
    opened2 = cmd("open", {"path": str(verified)})
    out2 = TMP_ROOT / "export-verified.h5ad"
    cmd("export", {"dataset": opened2["dataset"]["id"], "what": "h5ad", "path": str(out2)})
    data2, layers2, uns2 = _read_back(out2)
    check("6d 已验证数据导出保留 layers['counts']", "counts" in layers2, str(layers2))
    ok_values = "counts" in layers2 and np.array_equal(
        np.asarray(data2.layers["counts"].todense() if hasattr(data2.layers["counts"], "todense")
                   else data2.layers["counts"]), counts)
    check("6e layers['counts'] 与原始计数逐值一致", bool(ok_values),
          str(uns2.get("counts_verified")))
    _fresh()
    reopened2 = cmd("open", {"path": str(out2)})
    rt2 = ENGINE.DATASETS[reopened2["dataset"]["id"]]
    check("6f round-trip 后仍是已验证 counts", rt2.counts_verified is True, str(rt2.counts_verified))


# ---------------------------------------------------------------------------
# 7. 合并：稳定样本 ID×barcode + counts 语义继承
# ---------------------------------------------------------------------------


def _write_csv(path, matrix, genes, cells):
    """matrix 是 cells×genes；文本矩阵文件是 genes×cells，写盘前转置。"""
    frame = pd.DataFrame(np.asarray(matrix).T, index=genes, columns=cells)
    frame.to_csv(path, index=True)


def test_merge_unique_keys():
    print("\n[7] 合并导入的细胞唯一键与 counts 语义")
    genes = ["G1", "G2", "G3"]
    cells = ["c1", "c2", "c3", "c4"]
    core = np.array([[5, 1, 3], [0, 7, 1], [2, 2, 0], [4, 1, 2]], dtype=np.float64)
    lognorm = np.log1p(core / core.sum(axis=1, keepdims=True) * 1e4)
    dir_a = TMP_ROOT / "merge-a"
    dir_b = TMP_ROOT / "merge-b"
    dir_c = TMP_ROOT / "merge-c"
    for folder in (dir_a, dir_b, dir_c):
        folder.mkdir(parents=True, exist_ok=True)
    # 两个不同目录、同名 basename → basename 拼 barcode 会产生重复键
    _write_csv(dir_a / "sample.csv", core, genes, cells)
    _write_csv(dir_b / "sample.csv", core, genes, cells)
    _write_csv(dir_c / "expr.csv", lognorm, genes, cells)
    sources = [
        {"source": str(dir_a / "sample.csv"), "name": "S1"},
        {"source": str(dir_b / "sample.csv"), "name": "S2"},
    ]
    merged = SC_META.merge_sources(SC_IO, sources, str(TMP_ROOT / "merge-work"))
    index = [str(item) for item in merged["obs"].index]
    check("7a 同名文件合并后细胞键唯一", len(index) == len(set(index)), "%d/%d" % (len(set(index)), len(index)))
    check("7b 细胞键使用稳定样本 ID×barcode",
          {"S1|c1", "S2|c1"} <= set(index), str(index[:4]))
    check("7c sample 列保留自定义样本名",
          sorted(set(merged["obs"]["sample"].tolist())) == ["S1", "S2"],
          str(sorted(set(merged["obs"]["sample"].tolist()))))
    check("7d 样例都验证为 counts 时合并 counts 可用且已验证",
          merged.get("counts") is not None and (merged.get("meta") or {}).get("counts_verified") is True,
          str((merged.get("meta") or {}).get("counts_verified")))

    # 部分样本没有可验证 counts
    mixed = SC_META.merge_sources(
        SC_IO,
        [{"source": str(dir_a / "sample.csv"), "name": "S1"},
         {"source": str(dir_c / "expr.csv"), "name": "S3"}],
        str(TMP_ROOT / "merge-work"),
    )
    check("7e 混合样本时合并不再输出假 counts",
          mixed.get("counts") is None, str(type(mixed.get("counts"))))
    check("7f 混合样本时 meta 标记 counts_verified=False",
          (mixed.get("meta") or {}).get("counts_verified") is False,
          str((mixed.get("meta") or {}).get("counts_verified")))

    _fresh()
    opened = cmd("open", {"paths": [str(dir_a / "sample.csv"), str(dir_c / "expr.csv")]})
    dataset = ENGINE.DATASETS[opened["dataset"]["id"]]
    check("7g 引擎合并后进入 processed 模式", dataset.counts_verified is False,
          str(dataset.counts_verified))
    raised = None
    try:
        run_to_normalize(dataset.id)
    except Exception as error:  # noqa: BLE001
        raised = error
    check("7h 合并的 processed 数据拒绝 normalize",
          raised is not None
          and str(getattr(raised, "code", "")) == "prerequisite"
          and _counts_semantics_error(raised),
          getattr(raised, "message", raised))
    check("7i 合并后的细胞键仍然唯一",
          len(set(dataset.obs.index.tolist())) == dataset.obs.shape[0],
          "%d/%d" % (len(set(dataset.obs.index.tolist())), dataset.obs.shape[0]))


def main() -> int:
    print("== P0-1 原始层来源（counts provenance）回归 ==")
    if ad is None:
        print("  缺少 anndata：本测试需要 h5ad 支持")
        return 1
    for fn in (test_x_only_counts, test_lognorm_only_is_processed, test_fake_counts_full_scan,
               test_raw_log_with_counts_X, test_subset_inherits_scale, test_export_h5ad_layers,
               test_merge_unique_keys):
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
    else:
        print("fixtures: %s" % TMP_ROOT)
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
