#!/usr/bin/env Rscript
# dsh-geobulk — 差异表达 R 脚本（limma / DESeq2 官方实现）
#
# 调用：Rscript --vanilla lib/bg_de.R --in <bundle> --out <out> --params '<json>'
# bundle: expr.tsv（样本×基因）、info.csv（样本注释）、params.json
# 输出：<out>/result.json、<out>/tables/de.csv
# 退出码：0 成功；3 缺前置（stderr 给中文原因）；4 运行期程序错误（不可自动换方法）；其它失败
#
# 约定见 CONTRACT.md §3.2。图由 Python 侧统一绘制，本脚本只产出差异结果表。
#
# 修复记录（见《DSH_GEO_Singlecell_审查与优化任务单》）：
#   GEO-02  DESeq2 结果列名先统一成 logFC/P.Value/adj.P.Val 再过滤，且 NA 行按语义保留；
#           程序错误用退出码 4，与「缺包 / 数据不满足」区分，避免上层静默换检验。
#   GEO-03  设计矩阵用显式 data.frame（数值协变量保持 numeric），model.matrix 传 data；
#           配对以 pair 固定效应进入 limma/DESeq2，并校验缺失配对与设计满秩。

suppressPackageStartupMessages({
  ok_json <- requireNamespace("jsonlite", quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(key) {
  idx <- match(key, args)
  if (is.na(idx) || idx >= length(args)) return(NULL)
  args[idx + 1]
}

fail_prereq <- function(message) {
  cat(message, "\n", file = stderr())
  quit(save = "no", status = 3)
}
# 运行期程序错误：上层不得据此改换统计方法
fail_runtime <- function(message) {
  cat(message, "\n", file = stderr())
  quit(save = "no", status = 4)
}
fail_error <- function(message) {
  cat(message, "\n", file = stderr())
  quit(save = "no", status = 1)
}

in_dir <- getarg("--in")
out_dir <- getarg("--out")
params_raw <- getarg("--params")
if (is.null(in_dir) || is.null(out_dir)) fail_prereq("R 脚本缺少 --in / --out 参数")
if (!dir.exists(in_dir)) fail_prereq(paste0("bundle 目录不存在：", in_dir))
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

params <- list()
if (!is.null(params_raw)) {
  params <- tryCatch({
    if (requireNamespace("jsonlite", quietly = TRUE)) jsonlite::fromJSON(params_raw, simplifyVector = TRUE)
    else list()
  }, error = function(e) list())
}

read_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.delim(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
}

expr <- read_tsv(file.path(in_dir, "expr.tsv"))
if (is.null(expr) || nrow(expr) == 0 || ncol(expr) == 0) fail_prereq("bundle 里没有可用的 expr.tsv")
info <- read_tsv(file.path(in_dir, "info.csv"))
if (is.null(info)) {
  info <- data.frame(row.names = rownames(expr))
}

method <- if (!is.null(params$method)) tolower(as.character(params$method)) else "limma"
group_col <- if (!is.null(params$group_col)) as.character(params$group_col) else ""
case_label <- if (!is.null(params$case)) as.character(params$case) else ""
control_label <- if (!is.null(params$control)) as.character(params$control) else ""
covars <- if (!is.null(params$covars)) as.character(unlist(params$covars)) else character(0)
paired <- isTRUE(params$paired)
pair_col <- if (!is.null(params$pair_col)) as.character(params$pair_col) else ""
# 矩阵尺度由 Python 侧显式传入：counts=原始计数（DESeq2 / limma-voom），log=已标准化尺度
matrix_kind <- if (!is.null(params$matrix_kind)) tolower(as.character(params$matrix_kind)) else ""
if (!nzchar(matrix_kind)) matrix_kind <- if (method == "deseq2") "counts" else "log"

if (group_col == "" || !(group_col %in% colnames(info))) {
  fail_prereq(paste0("临床/样本表里没有分组列：", group_col))
}
sample_order <- rownames(expr)
if (!all(sample_order %in% rownames(info))) {
  missing <- setdiff(sample_order, rownames(info))
  fail_prereq(paste0("样本注释缺少这些样本的分组：", paste(head(missing, 5), collapse = ", ")))
}
groups <- as.character(info[[group_col]])[match(sample_order, rownames(info))]

if (case_label == "" || !(case_label %in% groups)) {
  fail_prereq(paste0("实验组不存在或为空：", case_label))
}
if (control_label == "" || !(control_label %in% groups)) {
  fail_prereq(paste0("对照组不存在或为空：", control_label))
}

# 只保留参与比较的两组（Python 侧已按需筛过配对，这里再兜一层）
keep_samples <- groups %in% c(case_label, control_label)
expr <- expr[keep_samples, , drop = FALSE]
info <- info[match(rownames(expr), rownames(info)), , drop = FALSE]
groups <- groups[keep_samples]

if (sum(groups == case_label) < 2 || sum(groups == control_label) < 2) {
  fail_prereq("实验组或对照组样本数不足 2 个")
}
if (nrow(expr) < 4) fail_prereq("可用样本不足 4 个")

# ---------------------------------------------------------------------------
# 设计数据框（GEO-03）
#   * 显式 data.frame：数值协变量保持 numeric，分类协变量 factor（不再一律 factor 化）
#   * 缺失值直接拒绝，不让 NA 悄悄进模型
# ---------------------------------------------------------------------------
design_data <- data.frame(
  group = factor(groups, levels = c(control_label, case_label)),
  row.names = rownames(expr),
  stringsAsFactors = FALSE
)
covar_terms <- character(0)
for (covariate in covars) {
  if (!(covariate %in% colnames(info))) {
    fail_prereq(paste0("协变量不在样本表里：", covariate))
  }
  raw <- info[[covariate]]
  if (any(is.na(raw))) {
    fail_prereq(paste0("协变量「", covariate, "」含缺失值，请先补齐或剔除这些样本"))
  }
  col_name <- paste0("cov_", make.names(covariate))
  if (is.numeric(raw) && !is.logical(raw)) {
    values <- as.numeric(raw)
    if (length(unique(values)) < 2) {
      fail_prereq(paste0("协变量「", covariate, "」在参与比较的样本里只有一个取值，设计矩阵不可识别；",
                         "请去掉该协变量或补充样本（不会静默丢弃协变量）"))
    }
    design_data[[col_name]] <- values
  } else {
    values <- factor(as.character(raw))
    if (nlevels(values) < 2) {
      fail_prereq(paste0("协变量「", covariate, "」在参与比较的样本里只有一个取值，设计矩阵不可识别；",
                         "请去掉该协变量或补充样本（不会静默丢弃协变量）"))
    }
    design_data[[col_name]] <- values
  }
  covar_terms <- c(covar_terms, col_name)
}

if (paired) {
  if (pair_col == "" || !(pair_col %in% colnames(info))) {
    fail_prereq("配对分析需要提供 pair_col 且该列存在于样本表")
  }
  pair_values <- as.character(info[[pair_col]])
  if (any(is.na(pair_values)) || any(trimws(pair_values) == "") || any(tolower(pair_values) %in% c("na", "nan", "none"))) {
    fail_prereq(paste0("配对列「", pair_col, "」存在缺失的配对 ID，无法配对；请补齐或显式剔除这些样本"))
  }
  design_data$pair <- factor(pair_values)
  if (nlevels(design_data$pair) < 2) {
    fail_prereq(paste0("配对列「", pair_col, "」只有 1 个取值，无法做配对分析"))
  }
  # 每一对必须同时含两组样本（Python 侧一般已剔除，这里兜底检查）
  pair_table <- table(design_data$pair, design_data$group)
  incomplete <- rownames(pair_table)[pair_table[, as.character(control_label)] == 0 |
                                       pair_table[, as.character(case_label)] == 0]
  if (length(incomplete) > 0) {
    fail_prereq(paste0("这些配对缺少一组样本，无法配对分析：", paste(head(incomplete, 5), collapse = ", ")))
  }
}

design_formula <- paste(
  c("group", covar_terms, if (paired) "pair"),
  collapse = " + "
)

n_dropped_p <- 0L
frame <- NULL
method_used <- ""
design_text <- ""

check_full_rank <- function(design, label) {
  q <- qr(design)
  rank <- q$rank
  if (rank < ncol(design)) {
    aliased <- character(0)
    if (!is.null(q$pivot) && length(q$pivot) == ncol(design)) {
      aliased <- colnames(design)[q$pivot[(rank + 1):ncol(design)]]
    }
    detail <- if (length(aliased) > 0) paste0("（与前面项的共线项：", paste(aliased, collapse = ", "), "）") else ""
    fail_prereq(sprintf(
      "%s 设计矩阵不满秩（rank=%d < 列数=%d）：这些变量在当前数据上不可同时识别%s。常见原因：协变量在配对/分组内部恒定，或样本数不足；请去掉该协变量或补充样本",
      label, rank, ncol(design), detail
    ))
  }
  invisible(TRUE)
}

# P0-2：计数模型（DESeq2 / limma-voom）的 count 校验 —— 负数不 clip、分数不 round
check_count_matrix <- function(mat, label) {
  if (any(!is.finite(mat))) {
    fail_prereq(paste0(label, " 需要有限的整数 counts；检测到 NaN/Inf，请先处理缺失值"))
  }
  if (any(mat < 0)) {
    fail_prereq(sprintf("%s 需要非负整数 counts；检测到负值（最小 %s），不会把负值截断成 0",
                        label, format(min(mat))))
  }
  frac <- abs(mat - round(mat)) > 1e-6
  if (any(frac)) {
    first_bad <- which(frac)[1]
    fail_prereq(sprintf(paste0("%s 需要整数 counts；检测到 %d 个非整数值（如 %s），",
                               "不会静默四舍五入；若这是标准化表达值请改用 limma"),
                        label, sum(frac), format(mat[first_bad])))
  }
  invisible(TRUE)
}

if (method == "deseq2") {
  if (!requireNamespace("DESeq2", quietly = TRUE)) fail_prereq("本机没有安装 DESeq2")
  suppressPackageStartupMessages(library(DESeq2))
  mat_raw <- t(as.matrix(expr))  # genes × samples
  # P0-2：负值/分数/非有限值一律拒绝，绝不 round() 静默圆整后照常出结果
  check_count_matrix(mat_raw, "DESeq2")
  mat <- round(mat_raw)
  storage.mode(mat) <- "integer"
  formula_text <- paste(c("group", covar_terms, if (paired) "pair"), collapse = " + ")
  check_full_rank(model.matrix(as.formula(paste("~", formula_text)), data = design_data), "DESeq2")
  design <- tryCatch(
    as.formula(paste("~", formula_text)),
    error = function(e) fail_runtime(paste0("DESeq2 设计公式无效：", conditionMessage(e)))
  )
  dds <- tryCatch(
    DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = design_data, design = design),
    error = function(e) fail_prereq(paste0("DESeq2 构建对象失败：", conditionMessage(e)))
  )
  dds <- tryCatch(
    DESeq2::DESeq(dds, quiet = TRUE),
    error = function(e) fail_prereq(paste0("DESeq2 拟合失败：", conditionMessage(e)))
  )
  res <- tryCatch(
    DESeq2::results(dds, contrast = c("group", case_label, control_label)),
    error = function(e) fail_runtime(paste0("DESeq2::results 失败：", conditionMessage(e)))
  )
  # GEO-02：先把 DESeq2 的结果列统一成内部列名，再做任何过滤/统计
  frame <- data.frame(
    gene = rownames(res),
    logFC = as.numeric(res$log2FoldChange),
    AveExpr = as.numeric(res$baseMean),
    t = as.numeric(res$stat),
    P.Value = as.numeric(res$pvalue),
    adj.P.Val = as.numeric(res$padj),
    B = rep(0, nrow(res)),
    stringsAsFactors = FALSE
  )
  method_used <- "deseq2"
  design_text <- paste("~", formula_text)
} else {
  if (!requireNamespace("limma", quietly = TRUE)) fail_prereq("本机没有安装 limma")
  suppressPackageStartupMessages(library(limma))
  # P0-2：真 counts 走 limma-voom（不能再 log2(x+1) 后当普通 limma）
  use_voom <- identical(matrix_kind, "counts")
  mat_raw <- t(as.matrix(expr))  # genes × samples
  voom_norm <- "none"
  voom_expr <- NULL
  if (use_voom) {
    check_count_matrix(mat_raw, "limma-voom")
    mat <- round(mat_raw)
    storage.mode(mat) <- "integer"
    # 标准 limma-voom 流程：TMM 归一化后再 voom（只按总计数归一在组成偏移时会有偏）
    if (requireNamespace("edgeR", quietly = TRUE)) {
      dge <- edgeR::DGEList(counts = mat)
      dge <- edgeR::calcNormFactors(dge, method = "TMM")
      voom_expr <- dge
      voom_norm <- "TMM"
    } else {
      voom_expr <- mat
    }
  } else {
    mat <- mat_raw
  }
  design <- tryCatch(
    model.matrix(as.formula(paste("~ 0 +", design_formula)), data = design_data),
    error = function(e) fail_runtime(paste0("limma 设计矩阵构建失败：", conditionMessage(e)))
  )
  check_full_rank(design, "limma")
  fit <- tryCatch(
    if (use_voom) limma::lmFit(limma::voom(voom_expr, design), design) else limma::lmFit(mat, design),
    error = function(e) fail_prereq(paste0("limma::lmFit 失败：", conditionMessage(e)))
  )
  # 分组列名形如 groupTumor（model.matrix 会把水平改成语法合法名）
  group_cols <- grep("^group", colnames(design), value = TRUE)
  pick_group_col <- function(label) {
    wanted <- make.names(label)
    hit <- group_cols[make.names(sub("^group", "", group_cols)) == wanted]
    if (length(hit) != 1) return(NA_character_)
    hit
  }
  case_col <- pick_group_col(case_label)
  ctrl_col <- pick_group_col(control_label)
  if (is.na(case_col) || is.na(ctrl_col)) {
    fail_prereq(sprintf("无法在 limma 设计矩阵里定位分组列（%s / %s）：%s",
                        case_label, control_label, paste(colnames(design), collapse = ", ")))
  }
  contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
  rownames(contrast_matrix) <- colnames(design)
  contrast_matrix[case_col, 1] <- 1
  contrast_matrix[ctrl_col, 1] <- -1
  colnames(contrast_matrix) <- paste0(case_label, "-", control_label)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  table <- limma::topTable(fit2, number = Inf, sort.by = "none")
  frame <- data.frame(
    gene = rownames(table),
    logFC = as.numeric(table$logFC),
    AveExpr = as.numeric(table$AveExpr),
    t = as.numeric(table$t),
    P.Value = as.numeric(table$P.Value),
    adj.P.Val = as.numeric(table$adj.P.Val),
    B = as.numeric(table$B),
    stringsAsFactors = FALSE
  )
  method_used <- if (use_voom) "limma-voom" else "limma"
  design_text <- paste("~ 0 +", design_formula)
}

if (is.null(frame) || nrow(frame) == 0) fail_runtime("R 差异分析没有产生结果行")

# GEO-02：DESeq2 独立过滤产生的 NA 行按语义保留并标记，不当成 0 行清空
n_dropped_p <- sum(is.na(frame$P.Value))
frame$regulation <- ifelse(
  !is.na(frame$adj.P.Val) & frame$adj.P.Val <= 0.05 & !is.na(frame$logFC) & frame$logFC > 0, "up",
  ifelse(!is.na(frame$adj.P.Val) & frame$adj.P.Val <= 0.05 & !is.na(frame$logFC) & frame$logFC < 0, "down", "ns")
)
if (nrow(frame) == 0 || ncol(frame) == 0) fail_runtime("R 差异结果在过滤后为空（不应发生，请报告此问题）")

write.csv(frame, file.path(out_dir, "tables", "de.csv"), row.names = FALSE, quote = FALSE)

n_up <- sum(frame$regulation == "up")
n_down <- sum(frame$regulation == "down")
rdir <- if (method_used == "deseq2") "DESeq2" else if (method_used == "limma-voom") "limma-voom" else "limma"
bullets <- c(
  sprintf("%s：%s vs %s（log2FC 方向 = %s - %s；%d vs %d 样本）",
          rdir, case_label, control_label, case_label, control_label,
          sum(groups == case_label), sum(groups == control_label)),
  sprintf("差异基因：上调 %d、下调 %d（FDR<0.05）", n_up, n_down)
)
if (paired) bullets <- c(bullets, sprintf("配对设计：%d 对（%s）已进入模型", nlevels(design_data$pair), pair_col))
if (n_dropped_p > 0) bullets <- c(bullets, sprintf("有 %d 个基因无可用 P 值（DESeq2 独立过滤/全零计数），保留在表中但判为 ns", n_dropped_p))

result <- list(
  summary = list(
    title = if (method_used == "deseq2") "差异表达（DESeq2）" else if (method_used == "limma-voom") "差异表达（limma-voom）" else "差异表达（limma）",
    bullets = bullets,
    metrics = list(
      method = method_used,
      requested_method = method,
      actual_method = method_used,
      design = design_text,
      contrast = paste0(case_label, "-", control_label),
      matrix_kind = matrix_kind,
      voom_normalization = if (method_used == "limma-voom") voom_norm else NA_character_,
      value_kind = if (method_used %in% c("deseq2", "limma-voom")) "counts" else "lognorm",
      n_effective = nrow(expr),
      paired = paired,
      pair_col = pair_col,
      covars = I(covar_terms),
      n_genes = nrow(frame),
      n_na_p = n_dropped_p,
      n_up = n_up,
      n_down = n_down,
      n_case = sum(groups == case_label),
      n_ctrl = sum(groups == control_label),
      n_pairs = if (paired) nlevels(design_data$pair) else 0
    )
  ),
  tables = list(list(
    id = "de",
    title = "差异基因",
    columns = lapply(c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "regulation"),
                     function(name) list(key = name, label = name,
                                         type = if (name %in% c("gene", "regulation")) "string" else "number")),
    rows = lapply(seq_len(min(nrow(frame), 500)), function(i) as.list(frame[i, ])),
    total = nrow(frame)
  )),
  warnings = list(),
  notes = list(sprintf("%s 官方 R 实现；design = %s；对比方向 = %s - %s",
                       rdir, design_text, case_label, control_label))
)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(result, file.path(out_dir, "result.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
} else {
  # 没有 jsonlite 时手写最小 JSON
  writeLines(
    sprintf('{"summary":{"title":"%s","bullets":["R 分析完成"],"metrics":{"method":"%s","n_genes":%d}}}',
            result$summary$title, method_used, nrow(frame)),
    file.path(out_dir, "result.json")
  )
}

quit(save = "no", status = 0)
