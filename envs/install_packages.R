# ==============================================================================
# 代谢组学复现项目：最小隔离依赖检查与安装脚本
#
# 设计说明:
#   - 采用 "核心基础组 (core) + 项目 Profile 扩展组 (metabonaut / transwell / spamtp)" 结构，
#     避免一次性安装无关 Bioconductor 大包引发依赖冲突。
#   - 特殊包 MsIO (Bioconductor/RforMassSpectrometry 生态新型序列化包) 若在当前镜像未收录，
#     将自动回退至 GitHub 源码安装。
#
# 用法:
#   Rscript envs/install_packages.R --dry-run   # 仅检查当前缺失包，不执行安装
#   Rscript envs/install_packages.R             # 自动安装缺失的 CRAN & Bioconductor 包
#   Rscript envs/install_packages.R --profile metabonaut
#   Rscript envs/install_packages.R --profile transwell
# ==============================================================================

# 解析命令行参数，支持 --dry-run 只读体检及 --profile 指定项目配置
args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

profile <- "core"
if ("--profile" %in% args) {
  idx <- which(args == "--profile")
  if (idx < length(args)) profile <- args[idx + 1]
} else if (any(grepl("^--profile=", args))) {
  profile <- sub("^--profile=", "", args[grep("^--profile=", args)][1])
}

# 核心通用基础包：可视化、表格读取与质谱定量核心容器
cran_core <- c(
  "BiocManager", "tidyverse", "data.table", "ggplot2",
  "patchwork", "pheatmap", "RColorBrewer", "jsonlite", "here", "readxl",
  "remotes", "pander", "vioplot", "gridExtra", "matrixStats", "cowplot"
)

bioc_core <- c(
  "Biobase", "SummarizedExperiment", "BiocParallel", "MsCoreUtils",
  "Spectra", "MsExperiment", "xcms", "MetaboAnnotation",
  "CompoundDb", "notame", "notameStats", "notameViz",
  "alabaster.se", "limma"
)

profile_extras <- switch(
  profile,
  "metabonaut" = list(
    cran = c("ggfortify", "ggVennDiagram", "UpSetR", "missForest", "vioplot"),
    bioc = c("MsBackendMetaboLights", "MsBackendMgf", "MsBackendMassbank", "AnnotationHub", "MetaboCoreUtils")
  ),
  "transwell" = list(
    cran = c("ggforce", "fuzzyjoin", "cowplot", "openxlsx"),
    bioc = c("ComplexHeatmap", "MSnbase", "mzR")
  ),
  "spamtp" = list(
    cran = c("cowplot", "ggnewscale", "ggridges", "FNN", "Seurat"),
    bioc = c("Cardinal", "EBImage", "edgeR", "fgsea", "limma", "MSnbase", "rgoslin", "scater")
  ),
  list(cran = character(0), bioc = character(0))
)

target_cran <- unique(c(cran_core, profile_extras$cran))
target_bioc <- unique(c(bioc_core, profile_extras$bioc))

installed <- rownames(installed.packages())
missing_cran <- setdiff(target_cran, installed)
missing_bioc <- setdiff(target_bioc, installed)

cat("=== 代谢组学环境依赖检查 (Profile:", profile, ") ===\n")
cat("R 版本:", R.version.string, "\n")
cat("CRAN 目标:", length(target_cran), "个 | 已装:", length(target_cran) - length(missing_cran),
    "| 缺失:", length(missing_cran), "\n")
if (length(missing_cran) > 0) cat("  缺失 CRAN:", paste(missing_cran, collapse = ", "), "\n")

cat("Bioc 目标:", length(target_bioc), "个 | 已装:", length(target_bioc) - length(missing_bioc),
    "| 缺失:", length(missing_bioc), "\n")
if (length(missing_bioc) > 0) cat("  缺失 Bioc:", paste(missing_bioc, collapse = ", "), "\n")

if (dry_run) {
  cat("\n[--dry-run] 仅完成检查，未执行任何包安装。\n")
  quit(status = 0)
}

if (length(missing_cran) > 0) {
  cat("\n>>> 安装缺失 CRAN 包...\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (length(missing_bioc) > 0) {
  cat("\n>>> 安装缺失 Bioconductor 包...\n")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

if (profile == "metabonaut" && !("MsIO" %in% rownames(installed.packages()))) {
  cat("\n>>> 安装 MsIO (优先 Bioconductor，若无则 GitHub RforMassSpectrometry/MsIO)...\n")
  tryCatch({
    BiocManager::install("MsIO", update = FALSE, ask = FALSE)
  }, error = function(e) {
    remotes::install_github("RforMassSpectrometry/MsIO", upgrade = "never")
  })
  if (!("MsIO" %in% rownames(installed.packages()))) {
    remotes::install_github("RforMassSpectrometry/MsIO", upgrade = "never")
  }
}

cat("\n=== 依赖安装流程结束 ===\n")
