# ==============================================================================
# 代谢组学复现工程总调度器 (Metabonaut Reproduction Runner)
#
# 设计说明:
#   - 通过 Rscript 命令行 --file= 自动定位项目根目录，无论从仓库任何位置调用均能稳定解析相对路径。
#   - 全部子步骤依次以独立 Rscript 子进程运行，任一环节退出码非 0 即中止流水线并报错。
#   - 全流程运行完毕后自动汇总所有核心软件包版本号并写入 repro_results/environment_versions.tsv 与 logs/sessionInfo.txt。
#
# 使用方式:
#   Rscript reproduction/metabolomics/Metabonaut/repro_scripts/run_all.R              # 默认 full 模式 (01-04)
#   Rscript reproduction/metabolomics/Metabonaut/repro_scripts/run_all.R --mode smoke # 仅运行 00_smoke_test.R
#   Rscript reproduction/metabolomics/Metabonaut/repro_scripts/run_all.R --mode full  # 全流程: 01 -> 02 -> 03 -> 04
# ==============================================================================

suppressPackageStartupMessages({
  library(utils)
})

# 解析 --mode 参数 (支持 "--mode smoke" 与 "--mode=smoke" 两种写法)
args <- commandArgs(trailingOnly = TRUE)
mode <- "full"
if ("--mode" %in% args) {
  idx <- which(args == "--mode")
  if (idx < length(args)) mode <- args[idx + 1]
} else if (any(grepl("^--mode=", args))) {
  mode <- sub("^--mode=", "", args[grep("^--mode=", args)][1])
}

# 动态定位 Metabonaut/ 项目根目录，避免依赖人工预先 setwd()
get_project_root <- function() {
  c_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", c_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
    return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
  }
  getwd()
}

proj_root <- get_project_root()
setwd(proj_root)

cat("===================================================================\n")
cat("  Metabonaut (Metabolomics 2026) 全流程自动化复现调度器\n")
cat("  工作目录:", proj_root, "\n")
cat("  运行模式:", toupper(mode), "\n")
cat("===================================================================\n\n")

scripts_to_run <- if (mode == "smoke") {
  c("00_smoke_test.R")
} else {
  c(
    "01_metabonaut_download_and_qc.R",
    "02_metabonaut_xcms_preprocessing.R",
    "03_metabonaut_norm_limma_annotation.R",
    "04_notame_drift_pqn_selection.R"
  )
}

status_records <- list()

for (s in scripts_to_run) {
  s_path <- file.path("repro_scripts", s)
  cat(sprintf(">>> 正在执行: %s ...\n", s))
  t_start <- proc.time()

  res <- tryCatch({
    system2("Rscript", args = c(s_path), stdout = "", stderr = "")
  }, error = function(e) {
    cat("执行出错:", conditionMessage(e), "\n")
    return(1)
  })

  elapsed <- round((proc.time() - t_start)[3], 1)
  status_str <- if (res == 0) "SUCCESS" else "FAILED"
  cat(sprintf("<<< 完成: %s | 状态: %s | 耗时: %s 秒\n\n", s, status_str, elapsed))
  status_records[[s]] <- data.frame(Script = s, Status = status_str, Elapsed_sec = elapsed)
  if (res != 0) stop("调度流程中止于: ", s)
}

cat("=== 收集并导出代谢组环境软件包版本快照 ===\n")
# 覆盖 00-04 脚本直接加载的所有核心与辅助 CRAN/Bioconductor 软件包
pkgs <- c(
  "R", "BiocManager", "Biobase", "SummarizedExperiment", "alabaster.se",
  "MsExperiment", "MsIO", "MsBackendMetaboLights", "MsBackendMgf",
  "MsBackendMassbank", "xcms", "Spectra", "MetaboCoreUtils",
  "MetaboAnnotation", "CompoundDb", "AnnotationHub", "BiocParallel",
  "limma", "notame", "notameStats", "notameViz", "matrixStats",
  "tidyverse", "data.table", "ggplot2", "patchwork", "pheatmap",
  "ggfortify", "vioplot", "gridExtra", "RColorBrewer", "cowplot"
)

ver_list <- lapply(pkgs, function(p) {
  v <- if (p == "R") {
    paste0(R.version$major, ".", R.version$minor)
  } else {
    tryCatch(as.character(packageVersion(p)), error = function(e) "Not Installed")
  }
  src <- if (p == "R") {
    "Base R"
  } else if (p %in% c("MsIO")) {
    "GitHub / Bioconductor"
  } else if (p %in% c("Biobase", "SummarizedExperiment", "alabaster.se", "MsExperiment",
                      "MsBackendMetaboLights", "MsBackendMgf", "MsBackendMassbank",
                      "xcms", "Spectra", "MetaboCoreUtils", "MetaboAnnotation",
                      "CompoundDb", "AnnotationHub", "BiocParallel", "limma",
                      "notame", "notameStats", "notameViz")) {
    "Bioconductor 3.23"
  } else {
    "CRAN"
  }
  data.frame(Package = p, Version = v, Repository = src, stringsAsFactors = FALSE)
})

ver_df <- do.call(rbind, ver_list)
env_path <- file.path("repro_results", "environment_versions.tsv")
write.table(ver_df, env_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("环境版本快照已写入:", env_path, "\n")
print(ver_df, row.names = FALSE)

sink(file.path("logs", "sessionInfo.txt"))
print(sessionInfo())
sink()
