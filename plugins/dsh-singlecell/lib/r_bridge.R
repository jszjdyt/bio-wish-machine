#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell — R bridge (CONTRACT §8.3)
#
# 把 R 的单细胞对象（.rds / .RData / .rda）导出成 Python 侧能直接读的中立文件：
#   <out>/matrix.mtx     MatrixMarket，cells × genes（orientation 见 meta.json）
#   <out>/features.tsv   每行一个基因名
#   <out>/barcodes.tsv   每行一个细胞名
#   <out>/obs.csv        row.names = 细胞名，列 = 细胞元数据
#   <out>/meta.json      {status,class,assay,assays,layers,orientation,dim,n_obs,n_vars,
#                         has_counts,warnings,notes[,candidates]}
#
# 支持对象：Seurat v3/v4/v5、SingleCellExperiment / SummarizedExperiment、
#           matrix / Matrix / dgCMatrix、data.frame（数值列）。
# 多对象歧义：meta.status="ambiguous" + candidates（object/class/dim/label），
#            Python 侧转成 EngineError("ambiguous")；--object 可指定对象名。
#
# 调用（由 lib/sc_io.py 负责）：
#   Rscript --vanilla r_bridge.R --in <file> --out <dir> [--object <name>] [--inspect]
#
# 约束：只用 base R + 已安装的 Seurat/SeuratObject/SingleCellExperiment/
#       SummarizedExperiment/Matrix；不依赖 jsonlite（手写 JSON，控制转义）。
#       所有文本输出 UTF-8、quote=FALSE（obs 里会破坏 CSV 的逗号/换行先替换成空格）。
# =============================================================================

options(warn = 1)
options(encoding = "UTF-8")
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || (is.character(a) && length(a) == 1 && is.na(a))) b else a

# ---------------------------------------------------------------------------
# 0. CLI
# ---------------------------------------------------------------------------

argv <- commandArgs(trailingOnly = TRUE)

flag <- function(name, default = NULL) {
  eq <- paste0("--", name, "=")
  hit <- which(startsWith(argv, eq))
  if (length(hit) > 0) return(substring(argv[hit[1]], nchar(eq) + 1L))
  hit <- which(argv == paste0("--", name))
  if (length(hit) > 0) {
    idx <- hit[1]
    if (idx < length(argv) && !startsWith(argv[idx + 1], "--")) return(argv[idx + 1])
    return("1")
  }
  default
}

IN <- flag("in")
OUT <- flag("out")
OBJ <- flag("object")
ASSAY <- flag("assay", "")
LAYER <- flag("layer", "")
SELFTEST <- !is.null(flag("selftest"))
INSPECT <- !is.null(flag("inspect"))

if (!SELFTEST && (is.null(IN) || is.null(OUT))) {
  cat("[r-bridge] 缺少参数：--in <file> --out <dir>\n", file = stderr())
  quit(status = 2L)
}

if (!SELFTEST && !dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. JSON 输出（手写；不依赖 jsonlite）
# ---------------------------------------------------------------------------

json_escape <- function(s) {
  s <- enc2utf8(as.character(s))
  if (is.na(s)) return("null")
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  paste0("\"", s, "\"")
}

json_num <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1 || is.na(x) || !is.finite(x)) return("null")
  if (x == round(x) && abs(x) < 1e15) return(format(x, scientific = FALSE, trim = TRUE))
  format(x, scientific = FALSE, digits = 15, trim = TRUE)
}

json_value <- function(x) {
  if (is.null(x)) return("null")
  if (is.data.frame(x)) x <- as.list(x)
  if (is.list(x)) {
    nms <- names(x)
    if (!is.null(nms) && length(nms) == length(x) && all(nzchar(nms))) {
      keep <- !vapply(x, is.null, logical(1))
      x <- x[keep]
      nms <- nms[keep]
      parts <- vapply(seq_along(x),
                      function(i) paste0(json_escape(nms[i]), ":", json_value(x[[i]])),
                      character(1))
      return(paste0("{", paste(parts, collapse = ","), "}"))
    }
    parts <- vapply(x, json_value, character(1))
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    if (length(x) == 0) return("[]")
    parts <- vapply(x, json_escape, character(1))
    if (length(parts) == 1) return(parts)
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (is.logical(x)) {
    if (length(x) == 0) return("[]")
    parts <- ifelse(is.na(x), "null", ifelse(x, "true", "false"))
    if (length(parts) == 1) return(parts)
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (is.numeric(x)) {
    if (length(x) == 0) return("[]")
    parts <- vapply(x, json_num, character(1))
    if (length(parts) == 1) return(parts)
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  json_escape(paste(as.character(x), collapse = " "))
}

write_meta <- function(meta, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(json_value(meta), con, useBytes = TRUE)
}

bridge_error <- function(msg) {
  cat(sprintf("[r-bridge] ERROR: %s\n", msg), file = stderr())
  try(write_meta(list(status = "error", error = msg), file.path(OUT, "meta.json")), silent = TRUE)
  quit(save = "no", status = 1L, runLast = FALSE)
}

# ---------------------------------------------------------------------------
# 2. 对象识别工具
# ---------------------------------------------------------------------------

class_label <- function(x) paste(class(x), collapse = ",")

is_supported_object <- function(x) {
  inherits(x, "Seurat") ||
    inherits(x, "SingleCellExperiment") ||
    inherits(x, "SummarizedExperiment") ||
    is.matrix(x) ||
    inherits(x, "Matrix") ||
    inherits(x, "DelayedMatrix") ||
    is.data.frame(x)
}

is_high_level <- function(x) {
  inherits(x, "Seurat") || inherits(x, "SingleCellExperiment") || inherits(x, "SummarizedExperiment")
}

obj_dim <- function(x) {
  d <- tryCatch(dim(x), error = function(e) NULL)
  if (is.null(d) || length(d) != 2L) return(NULL)
  as.integer(d)
}

# 统一表达成 c(cells, genes)；data.frame 按「行=细胞」估计（只是候选提示用）。
cell_gene_dim <- function(x) {
  d <- obj_dim(x)
  if (is.null(d)) return(NULL)
  if (is.data.frame(x)) return(c(cells = d[1], genes = d[2]))
  c(cells = d[2], genes = d[1])
}

candidate_entry <- function(name, x) {
  d <- cell_gene_dim(x)
  edim <- if (is.null(d)) list(NA_integer_, NA_integer_) else list(unname(d["cells"]), unname(d["genes"]))
  lab <- if (is.null(d)) {
    sprintf("%s · %s · 维度未知", name, class_label(x))
  } else {
    sprintf("%s · %s · %d 细胞 × %d 基因", name, class_label(x), d["cells"], d["genes"])
  }
  list(object = name, class = class_label(x), dim = edim, label = lab)
}

is_barcode_vec <- function(x) {
  x <- as.character(x)
  grepl("^[ACGTN]{10,}", x) | grepl("-[0-9]+$", x)
}

looks_barcode <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) return(FALSE)
  mean(is_barcode_vec(x), na.rm = TRUE) >= 0.6
}

looks_gene_symbol <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) return(FALSE)
  mask <- grepl("^[A-Za-z][A-Za-z0-9._-]{0,40}$", x) &
    !is_barcode_vec(x) &
    !grepl("^(GSM|GSE|SRR|ERR|DRR|SAM|SAMPLE)[-_]?[0-9]+", x, ignore.case = TRUE)
  mean(mask, na.rm = TRUE) >= 0.6
}

looks_sample_id <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) return(FALSE)
  mask <- grepl("^(GSM|GSE|SRR|ERR|DRR|SAM|SAMPLE)[-_]?[0-9]+", x, ignore.case = TRUE) |
    grepl("^[0-9]+$", x)
  mean(mask, na.rm = TRUE) >= 0.6
}

# 返回当前矩阵方向："genes_cells" 或 "cells_genes"。
guess_orientation <- function(rows, cols, default = "genes_cells") {
  if (is.null(rows) || is.null(cols) || length(rows) == 0 || length(cols) == 0) return(default)
  rb <- looks_barcode(rows)
  cb <- looks_barcode(cols)
  if (rb && !cb) return("cells_genes")
  if (cb && !rb) return("genes_cells")
  rg <- looks_gene_symbol(rows)
  cg <- looks_gene_symbol(cols)
  if (rg && !cg) return("genes_cells")
  if (cg && !rg) return("cells_genes")
  rs <- looks_sample_id(rows)
  cs <- looks_sample_id(cols)
  if (rs && !cs) return("cells_genes")
  if (cs && !rs) return("genes_cells")
  default
}

is_default_names <- function(x, n) {
  if (is.null(x) || length(x) != n) return(TRUE)
  x <- as.character(x)
  identical(x, as.character(seq_len(n)))
}

axis_names <- function(x, n, prefix, warnings) {
  bad <- is.null(x) || length(x) != n || any(!nzchar(as.character(x))) || is_default_names(x, n)
  if (bad) {
    warnings <- c(warnings, sprintf("%s 名缺失或无效，已自动生成 %s1..%sN", prefix, prefix, prefix))
    return(list(names = paste0(prefix, seq_len(n)), warnings = warnings))
  }
  list(names = enc2utf8(as.character(x)), warnings = warnings)
}

# ---------------------------------------------------------------------------
# 3. 读入 .rds / .RData
# ---------------------------------------------------------------------------

read_any_r_object <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    return(list(kind = "rds", objects = list(object = obj)))
  }
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (length(loaded) == 0) return(list(kind = "rdata", objects = list()))
  objs <- lapply(loaded, function(nm) get(nm, envir = env))
  names(objs) <- loaded
  list(kind = "rdata", objects = objs)
}

# rds 里常见 list(seurat_a=..., seurat_b=...) 的打包方式，摊平成候选对象。
# 只展开「纯 list」（class(x) 恰好是 list），避免把 lm 之类的有类 list 拆坏。
expand_objects <- function(objs) {
  out <- list()
  for (nm in names(objs)) {
    x <- objs[[nm]]
    if (!is_supported_object(x) && is.list(x) && identical(class(x), "list")) {
      inner <- names(x)
      if (is.null(inner) || any(!nzchar(inner))) inner <- paste0(nm, "_", seq_along(x))
      supported_inner <- vapply(x, is_supported_object, logical(1))
      if (any(supported_inner)) {
        for (i in which(supported_inner)) out[[inner[i]]] <- x[[i]]
        next
      }
    }
    out[[nm]] <- x
  }
  out
}

# 返回 list(status="ok", object=x) 或 list(status="ambiguous", candidates=[...])
select_object <- function(objs) {
  names_all <- names(objs)
  if (length(names_all) == 0) bridge_error("R 文件里没有找到任何对象")
  supported <- names_all[vapply(objs, is_supported_object, logical(1))]

  if (!is.null(OBJ) && nzchar(OBJ) && OBJ != "1") {
    if (!(OBJ %in% names_all)) {
      bridge_error(sprintf("对象 %s 不存在；可用对象：%s", OBJ, paste(names_all, collapse = ", ")))
    }
    if (!is_supported_object(objs[[OBJ]])) {
      bridge_error(sprintf("对象 %s 的类型 %s 不支持（支持 Seurat/SCE/matrix/Matrix/data.frame）",
                           OBJ, class_label(objs[[OBJ]])))
    }
    return(list(status = "ok", object = objs[[OBJ]]))
  }

  if (length(supported) == 0) {
    kinds <- paste(sprintf("%s(%s)", names_all, vapply(objs, class_label, character(1))), collapse = ", ")
    bridge_error(sprintf("R 文件里没有可识别的单细胞对象（找到：%s）", kinds))
  }

  high <- supported[vapply(objs[supported], is_high_level, logical(1))]
  pool <- if (length(high) > 0) high else supported
  if (length(pool) == 1) return(list(status = "ok", object = objs[[pool]]))

  list(status = "ambiguous",
       candidates = lapply(pool, function(nm) candidate_entry(nm, objs[[nm]])))
}

# ---------------------------------------------------------------------------
# 4. 提取：统一返回 genes × cells 的稀疏矩阵 + obs
# ---------------------------------------------------------------------------

seurat_assay_layers <- function(obj, assay) {
  layers <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character(0))
  if (length(layers) == 0) {
    # v3/v4 的 Assay 没有 Layers()，退化成固定槽名
    layers <- c("counts", "data", "scale.data")
  }
  layers
}

seurat_get_layer <- function(obj, assay, layer) {
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      "LayerData" %in% getNamespaceExports("SeuratObject")) {
    m <- tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = layer),
                  error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) return(m)
  }
  m <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, slot = layer),
                error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) return(m)
  m <- tryCatch({
    a <- obj@assays[[assay]]
    switch(layer,
           counts = a@counts,
           data = a@data,
           scale.data = a@scale.data,
           NULL)
  }, error = function(e) NULL)
  if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) return(m)
  NULL
}

# 显式 --layer 请求：精确层名优先；没有精确名时按前缀匹配拆分后的多 layer
pick_requested_layer <- function(layer_names, requested) {
  requested <- trimws(as.character(requested %||% ""))
  if (!nzchar(requested)) return(list(matched = character(0), strategy = "auto"))
  if (requested %in% layer_names) return(list(matched = requested, strategy = "exact"))
  prefixed <- layer_names[startsWith(layer_names, paste0(requested, "."))]
  if (length(prefixed)) return(list(matched = prefixed, strategy = "prefix"))
  list(matched = character(0), strategy = "missing")
}

# Seurat v5 多 layer 合并：必须按 feature ID 对齐（不能直接 cbind）
#   * feature 顺序/集合不同 → 按第一层的 feature 顺序取并集（缺失补 0）；
#   * 同一层内 feature 重名 → 无法可靠对齐，返回 NULL 并说明；
#   * 跨层 barcode 重名 → 用「层名|barcode」去重并记录。
align_layer_matrices <- function(mats, layer_names = NULL, min_cells = 1L) {
  empty <- list(matrix = NULL, features = character(0), alignment = "none",
                union_added = 0L, cells_renamed = 0L, warnings = character(0))
  mats <- mats[!vapply(mats, is.null, logical(1))]
  if (!length(mats)) return(empty)
  feature_lists <- lapply(mats, rownames)
  if (any(vapply(feature_lists, is.null, logical(1)))) {
    empty$warnings <- c(empty$warnings, "多 layer 中有矩阵缺少基因名（rownames），无法按 feature ID 对齐")
    return(empty)
  }
  for (index in seq_along(feature_lists)) {
    if (anyDuplicated(feature_lists[[index]])) {
      label <- if (!is.null(layer_names) && index <= length(layer_names)) layer_names[[index]] else as.character(index)
      empty$warnings <- c(empty$warnings,
                          sprintf("layer %s 的基因名有重复，无法可靠按 feature ID 对齐；已放弃该多 layer 合并", label))
      return(empty)
    }
  }
  features <- feature_lists[[1]]
  union_added <- 0L
  for (other in feature_lists[-1]) {
    extra <- setdiff(other, features)
    if (length(extra)) {
      features <- c(features, extra)
      union_added <- union_added + length(extra)
    }
  }
  rows <- list()
  for (index in seq_along(mats)) {
    mat <- mats[[index]]
    names <- feature_lists[[index]]
    if (!identical(names, features)) {
      aligned <- Matrix::Matrix(0, nrow = length(features), ncol = ncol(mat), sparse = TRUE)
      rownames(aligned) <- features
      colnames(aligned) <- colnames(mat)
      aligned[names, ] <- mat
      mat <- aligned
    }
    rows[[index]] <- mat
  }
  combined <- do.call(cbind, rows)
  cells_renamed <- 0L
  if (anyDuplicated(colnames(combined))) {
    labels <- rep(seq_along(rows), vapply(rows, ncol, integer(1)))
    combined_names <- colnames(combined)
    for (index in seq_along(combined_names)) {
      if (duplicated(combined_names)[index]) {
        layer_label <- if (!is.null(layer_names) && labels[index] <= length(layer_names)) {
          layer_names[[labels[index]]]
        } else {
          as.character(labels[index])
        }
        combined_names[index] <- sprintf("%s|%s", layer_label, combined_names[index])
        cells_renamed <- cells_renamed + 1L
      }
    }
    colnames(combined) <- combined_names
  }
  if (!is.null(min_cells) && ncol(combined) < min_cells) return(empty)
  list(matrix = combined, features = features,
       alignment = "feature_id_union", union_added = as.integer(union_added),
       cells_renamed = as.integer(cells_renamed), warnings = character(0))
}

extract_seurat <- function(obj, assay_requested = "", layer_requested = "") {
  warnings <- character(0)
  notes <- character(0)
  assays <- tryCatch(names(obj@assays), error = function(e) character(0))
  assay <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) NULL)
  if (is.null(assay) || !nzchar(assay) || !(assay %in% assays)) {
    assay <- if (length(assays) > 0) assays[1] else NULL
  }
  if (is.null(assay)) bridge_error("Seurat 对象里没有 assay")
  assay_requested <- trimws(as.character(assay_requested %||% ""))
  if (nzchar(assay_requested)) {
    if (assay_requested %in% assays) {
      assay <- assay_requested
    } else {
      warnings <- c(warnings, sprintf("请求的 assay %s 不存在（可用：%s），已回退 %s",
                                      assay_requested, paste(assays, collapse = ","), assay))
    }
  }
  if (length(assays) > 1 && !nzchar(assay_requested)) {
    notes <- c(notes, sprintf("Seurat 对象含 %d 个 assay（%s），已导出默认 assay %s",
                              length(assays), paste(assays, collapse = ","), assay))
  }
  layers <- seurat_assay_layers(obj, assay)

  m <- NULL
  has_counts <- FALSE
  layer_used <- NULL
  requested_pick <- pick_requested_layer(layers, layer_requested)
  if (nzchar(trimws(as.character(layer_requested %||% "")))) {
    if (requested_pick$strategy == "missing") {
      warnings <- c(warnings, sprintf("请求的 layer %s 不存在（可用：%s），已回退自动选择",
                                      layer_requested, paste(layers, collapse = ",")))
    } else if (requested_pick$strategy == "exact") {
      m <- seurat_get_layer(obj, assay, requested_pick$matched)
      if (!is.null(m)) {
        layer_used <- requested_pick$matched
        has_counts <- grepl("^counts", requested_pick$matched)
      }
    } else if (requested_pick$strategy == "prefix") {
      picked <- lapply(requested_pick$matched, function(name) seurat_get_layer(obj, assay, name))
      picked <- picked[!vapply(picked, is.null, logical(1))]
      if (length(picked)) {
        aligned <- align_layer_matrices(picked, requested_pick$matched)
        if (!is.null(aligned$matrix)) {
          m <- aligned$matrix
          layer_used <- paste(requested_pick$matched, collapse = "+")
          has_counts <- grepl("^counts", layer_requested)
          warnings <- c(warnings, aligned$warnings)
        }
      }
    }
  }
  if (is.null(m)) {
    m <- seurat_get_layer(obj, assay, "counts")
    has_counts <- !is.null(m)
    layer_used <- if (has_counts) "counts" else "data"
    if (is.null(m)) {
      m <- seurat_get_layer(obj, assay, "data")
      warnings <- c(warnings, "Seurat 对象没有 counts 层，已使用 data 层（X 可能已标准化）")
    }
  }
  feature_alignment <- "single_layer"

  # Seurat v5 多 layer：counts.1 / counts.2 ... 按列拼接
  if (is.null(m)) {
    count_layers <- grep("^counts(\\.|$)", layers, value = TRUE)
    if (length(count_layers) == 0) count_layers <- grep("^((data)|(scale\\.data))(\\.|$)", layers, value = TRUE)
    if (length(count_layers) > 0) {
      mats <- lapply(count_layers, function(l) seurat_get_layer(obj, assay, l))
      mats <- mats[!vapply(mats, is.null, logical(1))]
      used_layers <- count_layers[vapply(count_layers, function(l) !is.null(seurat_get_layer(obj, assay, l)), logical(1))]
      if (length(mats) > 0) {
        aligned <- align_layer_matrices(mats, used_layers)
        warnings <- c(warnings, aligned$warnings)
        if (!is.null(aligned$matrix)) {
          m <- aligned$matrix
          layer_used <- paste(used_layers, collapse = "+")
          has_counts <- grepl("^counts", layer_used)
          feature_alignment <- aligned$alignment
          notes <- c(notes, sprintf("Seurat v5 多 layer（%s）按 feature ID 对齐后合并：feature 并集共 %d 个（新增 %d），barcode 重名重命名 %d 个",
                                    paste(used_layers, collapse = ","), length(aligned$features),
                                    aligned$union_added, aligned$cells_renamed))
          if (aligned$union_added > 0 || aligned$cells_renamed > 0) {
            warnings <- c(warnings, sprintf("多 layer 的 feature/barcode 不一致（新增 feature %d 个、重命名 barcode %d 个），已按 feature ID 对齐",
                                            aligned$union_added, aligned$cells_renamed))
          }
        }
      }
    }
  }

  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) {
    bridge_error("Seurat 对象里没有可导出的表达矩阵（counts/data 都为空）")
  }

  genes <- rownames(m)
  cells <- colnames(m)
  obs <- tryCatch(obj@meta.data, error = function(e) NULL)

  list(matrix = m, genes = genes, cells = cells, obs = obs,
       class = class_label(obj), assay = assay, assays = assays, layers = layers,
       has_counts = has_counts, layer_used = layer_used,
       assay_requested = assay_requested %||% "",
       layer_requested = trimws(as.character(layer_requested %||% "")),
       feature_alignment = feature_alignment,
       warnings = warnings, notes = notes)
}

extract_sce <- function(obj) {
  warnings <- character(0)
  notes <- character(0)
  assays <- tryCatch(SummarizedExperiment::assayNames(obj), error = function(e) character(0))
  m <- NULL
  layer_used <- NULL
  has_counts <- FALSE
  if ("counts" %in% assays) {
    m <- tryCatch(SummarizedExperiment::assay(obj, "counts"), error = function(e) NULL,
                  warning = function(w) NULL)
    if (!is.null(m)) {
      has_counts <- TRUE
      layer_used <- "counts"
    }
  }
  if (is.null(m)) {
    first <- if (length(assays) > 0) assays[1] else 1
    m <- tryCatch(SummarizedExperiment::assay(obj, first), error = function(e) NULL,
                  warning = function(w) NULL)
    layer_used <- if (length(assays) > 0) first else "assay1"
    warnings <- c(warnings, "SingleCellExperiment 没有 counts assay，已使用第一个 assay（X 可能已标准化）")
  }
  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) {
    bridge_error("SingleCellExperiment 对象里没有可导出的 assay")
  }
  if (length(assays) > 1) {
    notes <- c(notes, sprintf("SingleCellExperiment 含 %d 个 assay（%s），已导出 %s",
                              length(assays), paste(assays, collapse = ","), layer_used))
  }
  obs <- tryCatch(as.data.frame(SummarizedExperiment::colData(obj)), error = function(e) NULL)
  list(matrix = m, genes = rownames(m), cells = colnames(m), obs = obs,
       class = class_label(obj), assay = layer_used, assays = assays, layers = assays,
       has_counts = has_counts, layer_used = layer_used,
       warnings = warnings, notes = notes)
}

extract_dataframe <- function(obj) {
  warnings <- character(0)
  notes <- character(0)
  numeric_cols <- vapply(obj, is.numeric, logical(1))
  dropped <- names(obj)[!numeric_cols]
  if (!any(numeric_cols)) bridge_error("data.frame 里没有数值列，无法作为表达矩阵")
  if (length(dropped) > 0) {
    warnings <- c(warnings, sprintf("data.frame 有 %d 个非数值列已忽略（%s）",
                                    length(dropped), paste(utils::head(dropped, 8), collapse = ",")))
  }
  raw <- as.matrix(obj[, numeric_cols, drop = FALSE])
  rownames(raw) <- rownames(obj)
  if (is.null(colnames(raw))) colnames(raw) <- names(obj)[numeric_cols]
  row_names <- rownames(raw)
  col_names <- colnames(raw)
  orient <- guess_orientation(row_names, col_names, default = "cells_genes")
  warnings <- c(warnings, sprintf("data.frame 方向启发式判定：%s",
                                  if (orient == "cells_genes") "行=细胞（默认）" else "行=基因，已转置"))
  if (orient == "cells_genes") {
    cells <- row_names
    genes <- col_names
    m <- Matrix::t(Matrix::Matrix(raw, sparse = TRUE))   # 统一成 genes × cells
  } else {
    genes <- row_names
    cells <- col_names
    m <- Matrix::Matrix(raw, sparse = TRUE)
  }
  list(matrix = m, genes = genes, cells = cells,
       obs = data.frame(row.names = cells),
       class = class_label(obj), assay = "", assays = character(0), layers = character(0),
       has_counts = TRUE, layer_used = "matrix",
       warnings = c(warnings, "矩阵没有独立的 counts 层，已把唯一矩阵同时作为 counts 与 X"),
       notes = notes)
}

extract_matrix_obj <- function(obj) {
  warnings <- character(0)
  notes <- character(0)
  m <- obj
  if (!is.matrix(m) && is.data.frame(m)) {
    return(extract_dataframe(m))
  }
  if (!is.matrix(m)) {
    m <- tryCatch(methods::as(m, "CsparseMatrix"),
                  error = function(e) tryCatch(as.matrix(m), error = function(e2) NULL))
  }
  if (is.null(m)) bridge_error(sprintf("无法把 %s 转成矩阵", class_label(obj)))
  row_names <- rownames(m)
  col_names <- colnames(m)
  orient <- guess_orientation(row_names, col_names, default = "genes_cells")
  warnings <- c(warnings, sprintf("矩阵方向启发式判定：%s",
                                  if (orient == "genes_cells") "行=基因（默认）" else "行=细胞，已转置"))
  if (is.null(row_names) || is.null(col_names)) {
    warnings <- c(warnings, "矩阵 dimnames 缺失，已自动生成基因名/细胞名")
  }
  if (orient == "cells_genes") {
    cells <- row_names
    genes <- col_names
    m <- Matrix::t(if (is.matrix(m)) Matrix::Matrix(m, sparse = TRUE) else m)  # genes × cells
  } else {
    genes <- row_names
    cells <- col_names
  }
  list(matrix = m, genes = genes, cells = cells,
       obs = data.frame(row.names = cells),
       class = class_label(obj), assay = "", assays = character(0), layers = character(0),
       has_counts = TRUE, layer_used = "matrix",
       warnings = c(warnings, "矩阵没有独立的 counts 层，已把唯一矩阵同时作为 counts 与 X"),
       notes = notes)
}

extract_any <- function(obj) {
  if (inherits(obj, "Seurat")) return(extract_seurat(obj, ASSAY %||% "", LAYER %||% ""))
  if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) return(extract_sce(obj))
  if (is.data.frame(obj)) return(extract_dataframe(obj))
  if (is.matrix(obj) || inherits(obj, "Matrix") || inherits(obj, "DelayedMatrix")) return(extract_matrix_obj(obj))
  bridge_error(sprintf("不支持的 R 对象类型：%s", class_label(obj)))
}

# 统一补基因名/细胞名、obs 对齐
finalize_extract <- function(result) {
  warnings <- result$warnings %||% character(0)
  n_genes <- nrow(result$matrix)
  n_cells <- ncol(result$matrix)

  genes <- axis_names(result$genes, n_genes, "gene_", warnings)
  warnings <- genes$warnings
  genes <- genes$names

  cells <- axis_names(result$cells, n_cells, "cell_", warnings)
  warnings <- cells$warnings
  cells <- cells$names

  obs <- result$obs
  if (is.null(obs) || nrow(obs) != n_cells) {
    if (!is.null(obs)) warnings <- c(warnings, "obs 行数与矩阵细胞数不一致，已只保留细胞名")
    obs <- data.frame(row.names = cells)
  } else {
    obs <- as.data.frame(obs, stringsAsFactors = FALSE)
    rn <- rownames(obs)
    if (is.null(rn) || any(!nzchar(rn))) {
      rownames(obs) <- cells
    } else if (setequal(rn, cells)) {
      obs <- obs[cells, , drop = FALSE]
    } else {
      warnings <- c(warnings, "obs 行名与细胞名不完全一致，已按顺序对齐")
      rownames(obs) <- cells
    }
  }

  result$matrix <- result$matrix
  result$genes <- genes
  result$cells <- cells
  result$obs <- obs
  result$warnings <- warnings
  result
}

# ---------------------------------------------------------------------------
# 5. 写输出文件
# ---------------------------------------------------------------------------

write_lines_utf8 <- function(x, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(as.character(x)), con, useBytes = TRUE)
}

to_csparse <- function(m) {
  if (is.matrix(m)) m <- Matrix::Matrix(m, sparse = TRUE)
  if (!inherits(m, "sparseMatrix")) m <- methods::as(m, "CsparseMatrix")
  if (inherits(m, "nMatrix") || inherits(m, "lMatrix") || !is.numeric(m@x)) m <- m * 1
  if (!inherits(m, "dgCMatrix")) {
    m <- methods::as(methods::as(m, "CsparseMatrix"), "dgCMatrix")
  }
  m
}

sanitize_matrix <- function(m) {
  warnings <- character(0)
  if (inherits(m, "sparseMatrix")) {
    if (!is.numeric(m@x)) m <- m * 1
    bad <- !is.finite(m@x)
    if (any(bad)) {
      m@x[bad] <- 0
      warnings <- c(warnings, sprintf("矩阵中 %d 个非有限值已替换为 0", sum(bad)))
    }
  } else {
    bad <- !is.finite(m)
    if (any(bad)) {
      m[bad] <- 0
      warnings <- c(warnings, sprintf("矩阵中 %d 个非有限值已替换为 0", sum(bad)))
    }
  }
  list(matrix = m, warnings = warnings)
}

write_obs_csv <- function(obs, cells, path) {
  warnings <- character(0)
  obs <- as.data.frame(obs, stringsAsFactors = FALSE)
  if (nrow(obs) != length(cells)) obs <- data.frame(row.names = cells)
  rownames(obs) <- cells
  if (ncol(obs) == 0) {
    file.create(path)
    return(warnings)
  }
  # 只保留可写成一列标量的列（list 列、matrix 列丢掉）
  keep <- vapply(obs, function(col) is.atomic(col) || is.factor(col), logical(1))
  dropped <- names(obs)[!keep]
  if (length(dropped) > 0) {
    warnings <- c(warnings, sprintf("obs 有 %d 个非标量列已忽略（%s）",
                                    length(dropped), paste(utils::head(dropped, 8), collapse = ",")))
    obs <- obs[, keep, drop = FALSE]
  }
  for (col in names(obs)) if (is.factor(obs[[col]])) obs[[col]] <- as.character(obs[[col]])
  bad_values <- 0
  for (col in names(obs)) {
    if (is.character(obs[[col]])) {
      hit <- grepl("[,\r\n]", obs[[col]])
      bad_values <- bad_values + sum(hit, na.rm = TRUE)
      obs[[col]] <- gsub("[,\r\n]", " ", obs[[col]])
    }
  }
  if (bad_values > 0) {
    warnings <- c(warnings, sprintf("obs 中 %d 个值含逗号/换行，已替换为空格以保证 CSV 可解析", bad_values))
  }
  names(obs) <- gsub("[,\r\n]", "_", names(obs))
  names(obs) <- make.unique(ifelse(nzchar(names(obs)), names(obs), "obs"))
  utils::write.table(obs, path, sep = ",", row.names = TRUE, col.names = NA,
                     quote = FALSE, na = "", fileEncoding = "UTF-8")
  warnings
}

export_extra <- function(obj, out_dir) {
  warnings <- character(0)
  state <- list(reductions = list(), graphs = list(), variable_features = list(),
                rank_genes = FALSE, embed_kind = "", class = class(obj)[1])
  reductions <- character(0)
  embeddings <- list()
  write_embedding_csv <- function(emb, file, max_dim, prefix) {
    if (is.null(emb) || is.null(rownames(emb)) || nrow(emb) < 1 || ncol(emb) < 1) return(FALSE)
    dims <- seq_len(min(max_dim, ncol(emb)))
    frame <- data.frame(cell = rownames(emb), emb[, dims, drop = FALSE], check.names = FALSE)
    names(frame)[-1] <- paste0(prefix, dims)
    write.csv(frame, file, row.names = FALSE, quote = FALSE)
    TRUE
  }
  if (inherits(obj, "Seurat")) {
    reductions <- tryCatch(as.character(SeuratObject::Reductions(obj)), error = function(e) character(0))
    state$reductions <- as.list(reductions)
    vf <- tryCatch(SeuratObject::VariableFeatures(obj), error = function(e) character(0))
    state$variable_features <- as.list(as.character(vf))
    assay <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) "")
    state$graphs <- as.list(tryCatch(as.character(names(SeuratObject::Graphs(obj, assay))), error = function(e) character(0)))
    misc_names <- tryCatch(as.character(names(obj@misc)), error = function(e) character(0))
    state$rank_genes <- any(grepl("marker|rank|deg|diff", misc_names, ignore.case = TRUE))
    for (name in reductions) {
      emb <- tryCatch(SeuratObject::Embeddings(obj, name), error = function(e) NULL)
      if (!is.null(emb)) embeddings[[name]] <- emb
    }
  } else if (inherits(obj, "SingleCellExperiment")) {
    reductions <- tryCatch(as.character(SingleCellExperiment::reducedDimNames(obj)), error = function(e) character(0))
    state$reductions <- as.list(reductions)
    for (name in reductions) {
      emb <- tryCatch(SingleCellExperiment::reducedDim(obj, name), error = function(e) NULL)
      if (!is.null(emb)) embeddings[[name]] <- emb
    }
    rd <- tryCatch(SummarizedExperiment::rowData(obj), error = function(e) NULL)
    if (!is.null(rd) && "highly_variable" %in% colnames(rd)) {
      state$variable_features <- as.list(rownames(rd)[which(as.logical(rd[["highly_variable"]]))])
    }
  }
  pca_names <- reductions[tolower(reductions) == "pca"]
  if (length(pca_names) > 0 && !is.null(embeddings[[pca_names[1]]])) {
    write_embedding_csv(embeddings[[pca_names[1]]], file.path(out_dir, "pca.csv"), 30, "PC")
  }
  kind <- ""
  for (candidate in c("umap", "UMAP", "tsne", "tSNE", "TSNE")) {
    if (candidate %in% reductions && !is.null(embeddings[[candidate]])) { kind <- candidate; break }
  }
  if (!nzchar(kind)) {
    for (name in reductions) {
      if (tolower(name) == "pca") next
      if (!is.null(embeddings[[name]]) && ncol(embeddings[[name]]) >= 2) { kind <- name; break }
    }
  }
  if (nzchar(kind) && !is.null(embeddings[[kind]]) && ncol(embeddings[[kind]]) >= 2) {
    write_embedding_csv(embeddings[[kind]], file.path(out_dir, "embed.csv"), 2, "dim")
    state$embed_kind <- kind
  }
  vf_values <- as.character(unlist(state$variable_features))
  vf_values <- vf_values[!is.na(vf_values) & nzchar(vf_values)]
  if (length(vf_values) > 0) write_lines_utf8(vf_values, file.path(out_dir, "variable_features.txt"))
  list(state = state, warnings = warnings)
}

write_outputs <- function(result, out_dir) {
  warnings <- character(0)
  san <- sanitize_matrix(result$matrix)
  warnings <- c(warnings, san$warnings)
  mm <- to_csparse(san$matrix)
  mm <- methods::as(Matrix::t(mm), "CsparseMatrix")   # cells × genes
  mm <- to_csparse(mm)
  Matrix::writeMM(mm, file.path(out_dir, "matrix.mtx"))

  write_lines_utf8(result$genes, file.path(out_dir, "features.tsv"))
  write_lines_utf8(result$cells, file.path(out_dir, "barcodes.tsv"))
  warnings <- c(warnings, write_obs_csv(result$obs, result$cells, file.path(out_dir, "obs.csv")))
  warnings
}

# ---------------------------------------------------------------------------
# 6. 轻量 inspect（只报 class/dim，不导出矩阵）
# ---------------------------------------------------------------------------

inspect_object <- function(obj) {
  d <- cell_gene_dim(obj)
  assays <- character(0)
  assay <- ""
  layers <- character(0)
  if (inherits(obj, "Seurat")) {
    assays <- tryCatch(names(obj@assays), error = function(e) character(0))
    assay <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) "") %||% ""
    if (nzchar(assay) && assay %in% assays) layers <- seurat_assay_layers(obj, assay)
  } else if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) {
    assays <- tryCatch(SummarizedExperiment::assayNames(obj), error = function(e) character(0))
    assay <- if (length(assays) > 0) assays[1] else ""
    layers <- assays
  }
  list(class = class_label(obj),
       assay = assay,
       assays = as.list(assays),
       layers = as.list(layers),
       dim = if (is.null(d)) list(NA_integer_, NA_integer_) else list(unname(d["cells"]), unname(d["genes"])),
       n_obs = if (is.null(d)) NA_integer_ else unname(d["cells"]),
       n_vars = if (is.null(d)) NA_integer_ else unname(d["genes"]))
}

# ---------------------------------------------------------------------------
# 7. main
# ---------------------------------------------------------------------------

run_selftest <- function() {
  checks <- list()
  add <- function(id, ok, detail = "") {
    checks[[length(checks) + 1]] <<- list(id = id, ok = isTRUE(ok), detail = as.character(detail))
  }
  # 1) feature 顺序/集合不同的两层：按第一层顺序取并集，值必须对得上
  mat_a <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2,
                                  dimnames = list(c("G1", "G2"), c("c1", "c2"))), sparse = TRUE)
  mat_b <- Matrix::Matrix(matrix(c(30, 10, 20, 40), nrow = 2,
                                  dimnames = list(c("G3", "G1"), c("c1", "c2"))), sparse = TRUE)
  aligned <- align_layer_matrices(list(mat_a, mat_b), c("counts.1", "counts.2"))
  expected <- matrix(0, nrow = 3, ncol = 4,
                     dimnames = list(c("G1", "G2", "G3"), c("c1", "c2", "c1", "c2")))
  expected["G1", ] <- c(1, 3, 10, 40)
  expected["G2", ] <- c(2, 4, 0, 0)
  expected["G3", ] <- c(0, 0, 30, 20)
  add("align_permuted_layers",
      !is.null(aligned$matrix) && isTRUE(all.equal(unname(as.matrix(aligned$matrix)), unname(expected))),
      sprintf("features=%s", paste(aligned$features, collapse = ",")))
  # 2) 并集 feature 数
  add("align_union_features",
      identical(length(aligned$features), 3L) && identical(aligned$union_added, 1L),
      sprintf("n=%d added=%d", length(aligned$features), aligned$union_added))
  # 3) 同层 feature 重名 → 必须放弃而不是错位
  dup <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 1,
                                dimnames = list(c("G1", "G1"), "c1")), sparse = TRUE)
  rejected <- align_layer_matrices(list(dup), "counts.bad")
  add("align_duplicate_features",
      is.null(rejected$matrix) && any(grepl("重复", rejected$warnings)),
      paste(rejected$warnings, collapse = " "))
  # 4) 跨层 barcode 重名 → 重命名且仍唯一
  dup_cells <- align_layer_matrices(list(mat_a, mat_a), c("counts.1", "counts.2"))
  add("align_duplicate_cells",
      !is.null(dup_cells$matrix) && !anyDuplicated(colnames(dup_cells$matrix)) && dup_cells$cells_renamed > 0,
      sprintf("renamed=%d", dup_cells$cells_renamed))
  # 5) 显式 --layer：精确名优先，拆分 layer 用前缀匹配
  pick_exact <- pick_requested_layer(c("counts", "data"), "data")
  pick_prefix <- pick_requested_layer(c("counts.1", "counts.2", "data.1"), "counts")
  add("explicit_layer_arg",
      identical(pick_exact$matched, "data") && length(pick_prefix$matched) == 2L &&
        identical(pick_prefix$strategy, "prefix"),
      sprintf("exact=%s prefix=%s", pick_exact$matched, paste(pick_prefix$matched, collapse = ",")))
  payload <- list(status = "ok", checks = checks)
  cat(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"), "\n", sep = "")
  if (!all(vapply(checks, function(item) isTRUE(item$ok), logical(1)))) quit(save = "no", status = 1L)
  quit(save = "no", status = 0L, runLast = FALSE)
}

main <- function() {
  if (!file.exists(IN)) bridge_error(sprintf("文件不存在：%s", IN))
  loaded <- tryCatch(read_any_r_object(IN),
                     error = function(e) bridge_error(sprintf("读取 R 文件失败：%s", conditionMessage(e))))
  objs <- expand_objects(loaded$objects)
  if (length(objs) == 0) bridge_error(sprintf("文件里没有对象：%s", IN))

  sel <- select_object(objs)
  if (identical(sel$status, "ambiguous")) {
    meta <- list(status = "ambiguous",
                 kind = loaded$kind,
                 candidates = sel$candidates,
                 warnings = list("R 文件里存在多个可用的单细胞对象，请用 object 参数指定"),
                 notes = list())
    write_meta(meta, file.path(OUT, "meta.json"))
    quit(save = "no", status = 0L, runLast = FALSE)
  }

  obj <- sel$object
  if (INSPECT) {
    info <- inspect_object(obj)
    meta <- c(list(status = "ok", kind = loaded$kind, inspect = TRUE,
                   orientation = "cells_genes", warnings = list(), notes = list()),
              info)
    write_meta(meta, file.path(OUT, "meta.json"))
    quit(save = "no", status = 0L, runLast = FALSE)
  }

  result <- finalize_extract(extract_any(obj))
  extra <- tryCatch(export_extra(obj, OUT),
                    error = function(e) list(state = list(), warnings = sprintf("状态导出失败：%s", conditionMessage(e))))
  extra_warnings <- write_outputs(result, OUT)
  extra_warnings <- c(extra_warnings, as.character(extra$warnings %||% character(0)))
  warnings <- unique(c(result$warnings, extra_warnings))
  meta <- list(
    status = "ok",
    kind = loaded$kind,
    class = result$class,
    assay = result$assay %||% "",
    assays = as.list(result$assays %||% character(0)),
    layers = as.list(result$layers %||% character(0)),
    layer_used = result$layer_used %||% "",
    orientation = "cells_genes",
    dim = list(ncol(result$matrix), nrow(result$matrix)),
    n_obs = ncol(result$matrix),
    n_vars = nrow(result$matrix),
    has_counts = isTRUE(result$has_counts),
    assay_requested = result$assay_requested %||% "",
    layer_requested = result$layer_requested %||% "",
    feature_alignment = result$feature_alignment %||% "",
    import_state = extra$state,
    warnings = as.list(warnings),
    notes = as.list(result$notes %||% character(0))
  )
  write_meta(meta, file.path(OUT, "meta.json"))
  quit(save = "no", status = 0L, runLast = FALSE)
}

if (SELFTEST) run_selftest()
main()
