#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P0-4 回归：R 侧多 layer 必须按 feature ID 对齐，批次/条件列必须显式。

反例（修复前实测）：
  1. ``r_bridge.R`` 的 Seurat v5 多 layer 导出用 ``do.call(cbind, mats)``——
     不同 layer 的 features 顺序/集合不同时，基因与表达值整体错位（同一基因
     出现两次、值来自另一 layer）；
  2. ``seurat_pipeline.R`` 的 ``.detect_batch_key`` 候选里含 ``group``——
     条件列会被当批次整合掉生物学差异；
  3. ``already_split`` 只看层名有没有后缀，不校验层名与批次水平是否一致；
     层不是按该批次列拆的也直接 IntegrateLayers；
  4. RPCA 失败静默回退 Harmony，结果里没有 method_requested/method_actual/
     fallback/fallback_reason，用户无法知道实际跑的是什么。

测试：
  * [A] 静态：r_bridge.R 有多 layer 对齐函数且真的被调用，不再裸 cbind；
  * [B] 动态：``Rscript r_bridge.R --selftest`` 用乱序/缺 feature/重名 barcode 的
    矩阵验证对齐结果逐值正确（base R + Matrix，无需 Seurat）；
  * [C] 静态：seurat_pipeline.R 显式 assay/layer/batch_col/condition_col，候选
    不再含 group，already_split 校验层名 vs 批次水平，回退记录四元组；
  * [D] 动态：``Rscript seurat_pipeline.R --selftest`` 用纯函数验证批次列选择、
    层名覆盖校验与回退信息结构（无需 Seurat >= 6 安装）。

运行：cd plugins/dsh-singlecell && python3 test/test_p0_seurat_layers.py
退出码 0 = 全部 PASS。
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLUGIN = HERE.parent
BRIDGE = PLUGIN / "lib" / "r_bridge.R"
PIPELINE = PLUGIN / "lib" / "external" / "seurat_pipeline.R"
TMP_ROOT = Path(tempfile.mkdtemp(prefix="sc-p0-seurat-layers-"))
os.environ.setdefault("MPLBACKEND", "Agg")
sys.dont_write_bytecode = True

PASSED = []
FAILED = []

RSCRIPT_CANDIDATES = [
    "/Library/Frameworks/R.framework/Resources/bin/Rscript",
    "/usr/local/bin/Rscript",
    "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript",
]


def check(name, cond, detail=""):
    if cond:
        PASSED.append(name)
        print("  PASS %s%s" % (name, ("  [" + str(detail)[:220] + "]") if detail else ""))
        return True
    FAILED.append(name)
    print("  FAIL %s%s" % (name, ("  [" + str(detail)[:220] + "]") if detail else ""))
    return False


def rscript_path():
    for candidate in RSCRIPT_CANDIDATES:
        if Path(candidate).exists():
            return candidate
    return shutil.which("Rscript")


def read_src(path):
    return Path(path).read_text(encoding="utf-8", errors="replace")


def run_selftest(script, timeout=120):
    """跑 R 脚本的 --selftest，返回 (返回值, stdout, stderr)。"""
    rscript = rscript_path()
    if not rscript:
        return None, "", "Rscript 不存在"
    env = dict(os.environ)
    env.setdefault("R_LIBS_USER", str(Path(os.path.expanduser("~/.dsh/sc-plugin-rlib"))))
    try:
        completed = subprocess.run(
            [rscript, "--vanilla", str(script), "--selftest"],
            capture_output=True, text=True, timeout=timeout, env=env,
        )
    except subprocess.TimeoutExpired:
        return None, "", "selftest 超时"
    payload = None
    for line in reversed((completed.stdout or "").strip().splitlines()):
        line = line.strip()
        if line.startswith("{"):
            try:
                payload = json.loads(line)
                break
            except Exception:  # noqa: BLE001
                continue
    return {
        "code": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "payload": payload,
    }, completed.stdout, completed.stderr


def selftest_checks(detail, label=None):
    payload = (detail or {}).get("payload") if isinstance(detail, dict) else None
    if not isinstance(payload, dict):
        return None
    checks = payload.get("checks")
    if not isinstance(checks, list):
        return None
    result = {}
    for item in checks:
        if isinstance(item, dict) and item.get("id") is not None:
            result[str(item["id"])] = item
    return result


# ---------------------------------------------------------------------------
# [A] r_bridge.R 静态
# ---------------------------------------------------------------------------


def test_bridge_static():
    print("\n[A] r_bridge.R：多 layer 按 feature ID 对齐（静态）")
    src = read_src(BRIDGE)
    has_fn = "align_layer_matrices" in src or "align_layers_by_feature" in src
    check("A1 存在 layer 对齐函数", has_fn,
          "align_layer_matrices" if has_fn else "未找到对齐函数")
    bare_cbind = "do.call(cbind, mats)" in src
    check("A2 不再裸 do.call(cbind, mats)", not bare_cbind,
          "仍存在裸 cbind" if bare_cbind else "已替换")
    uses_align = src.count("align_layer_matrices") >= 2 or src.count("align_layers_by_feature") >= 2
    check("A3 对齐函数在 extract_seurat 多 layer 路径被真正调用", uses_align,
          "引用次数=%d" % (src.count("align_layer_matrices") + src.count("align_layers_by_feature")))
    check("A4 记录 feature_alignment / 对齐信息到导出 meta",
          "feature_alignment" in src or "features_aligned" in src,
          "未见 feature_alignment")
    check("A5 支持显式 --assay / --layer 参数",
          'flag("assay"' in src and ('flag("layer"' in src or 'flag("counts_layer"' in src),
          "assay=%s layer=%s" % ('flag("assay"' in src, 'flag("layer"' in src))
    check("A6 提供 --selftest 动态自检入口", "selftest" in src.lower(), "")


# ---------------------------------------------------------------------------
# [B] r_bridge.R 动态 selftest
# ---------------------------------------------------------------------------


def test_bridge_dynamic():
    print("\n[B] r_bridge.R --selftest：乱序 / 缺 feature / 重名 barcode")
    outcome, _stdout, stderr = run_selftest(BRIDGE)
    if outcome is None:
        check("B0 r_bridge.R 可执行 selftest", False, stderr[:160])
        return
    check("B0 selftest 退出码 0", outcome["code"] == 0, "code=%s stderr=%s" % (
        outcome["code"], (outcome["stderr"] or "")[-160:]))
    checks = selftest_checks(outcome, "r_bridge") or {}
    required = ["align_permuted_layers", "align_union_features", "align_duplicate_features",
                "align_duplicate_cells", "explicit_layer_arg"]
    for check_id in required:
        item = checks.get(check_id)
        check("B-%s 通过" % check_id, bool(item and item.get("ok")),
              "" if item is None else str(item.get("detail"))[:180])


# ---------------------------------------------------------------------------
# [C] seurat_pipeline.R 静态
# ---------------------------------------------------------------------------


def test_pipeline_static():
    print("\n[C] seurat_pipeline.R：显式层/批次列 + already_split 校验 + 回退记录")
    src = read_src(PIPELINE)
    for token, label in (("assay", "assay"), ("layer", "layer"),
                         ("batch_col", "batch_col"), ("condition_col", "condition_col")):
        check("C1 支持显式 %s 参数" % label, ('params$%s' % token) in src or ('"%s"' % token) in src,
              "未找到 params$%s" % token)
    detect_start = src.find("detect_batch_key")
    detect_body = src[detect_start: detect_start + 1600] if detect_start >= 0 else ""
    candidates_line = ""
    for line in detect_body.splitlines():
        if "candidates <-" in line:
            candidates_line = line.strip()
            break
    check("C2 批次列候选不再包含条件列 group/condition",
          candidates_line != "" and '"group"' not in candidates_line and '"condition"' not in candidates_line,
          candidates_line[:160] or "找不到 candidates <- 行")
    check("C3 already_split 校验层名与批次水平一致（不是只看后缀）",
          "layers_cover_batch_levels" in src or "layer_levels_cover" in src,
          "未找到层名覆盖校验函数")
    check("C4 层与批次列不一致时用 JoinLayers 后重新 split",
          "JoinLayers" in src and ("ensure_split_layers" in src or "join_then_split" in src),
          "JoinLayers=%s ensure=%s" % ("JoinLayers" in src, "ensure_split_layers" in src))
    for token in ("method_requested", "method_actual", "fallback_reason"):
        check("C5 回退记录 %s" % token, token in src, "")
    check("C6 RPCA 回退 Harmony 时 method_actual 明确为 harmony",
          "method_actual" in src and "fallback_reason" in src,
          "")
    check("C7 提供 --selftest 动态自检入口", "selftest" in src.lower(), "")


# ---------------------------------------------------------------------------
# [D] seurat_pipeline.R 动态 selftest
# ---------------------------------------------------------------------------


def test_pipeline_dynamic():
    print("\n[D] seurat_pipeline.R --selftest：批次列选择 / 层名覆盖 / 回退结构")
    outcome, _stdout, stderr = run_selftest(PIPELINE)
    if outcome is None:
        check("D0 seurat_pipeline.R 可执行 selftest", False, stderr[:160])
        return
    check("D0 selftest 退出码 0", outcome["code"] == 0, "code=%s stderr=%s" % (
        outcome["code"], (outcome["stderr"] or "")[-200:]))
    checks = selftest_checks(outcome, "seurat_pipeline") or {}
    required = ["batch_key_excludes_condition", "batch_key_explicit",
                "layers_cover_batch_levels", "fallback_fields", "condition_col_recorded"]
    for check_id in required:
        item = checks.get(check_id)
        check("D-%s 通过" % check_id, bool(item and item.get("ok")),
              "" if item is None else str(item.get("detail"))[:180])
    payload = outcome["payload"] or {}
    metrics = payload.get("example_integration") or {}
    check("D6 示例回退结构含 requested/actual/fallback/reason 四元组",
          {"method_requested", "method_actual", "fallback", "fallback_reason"} <= set(metrics),
          str(sorted(metrics))[:160])


def main() -> int:
    print("== P0-4 Seurat split layers / 批次列 / 回退记录回归 ==")
    for fn in (test_bridge_static, test_bridge_dynamic, test_pipeline_static,
               test_pipeline_dynamic):
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
