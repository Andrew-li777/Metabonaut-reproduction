# ==============================================================================
# 02_metabonaut_xcms_preprocessing.R
# Metabonaut Step 2: xcms v4 色谱峰检测、峰精修、QC 锚点对齐、峰分组与缺峰回填
#
# 上游源码对应:
#   upstream_Metabonaut/vignettes/a-end-to-end-untargeted-metabolomics.qmd
#   (Preprocessing / Alignment / Correspondence / Gap-filling 章节)
# 输入文件:
#   - data/processed/01_lcms1_filtered/ (Step 01 导出的 MsExperiment)
#   - data/processed/01_qc_metadata.rds
# 输出产物:
#   - repro_results/04_peak_detection_summary.png (对应原著 Fig 9)
#   - repro_results/05_rt_alignment_diagnostics.png (对应原著 Fig 15)
#   - repro_results/06_gap_filling_qc_correlation.png (对应原著 Fig 23)
#   - repro_results/02_xcms_preprocessing_stats.tsv
#   - data/processed/02_preprocessed_lcms1/ (XcmsExperiment)
#   - data/processed/02_preprocessed_res/ (SummarizedExperiment, 供 Step 03 & Step 04)
# ==============================================================================

suppressPackageStartupMessages({
  library(MsExperiment)
  library(MsIO)
  library(alabaster.se)
  library(Spectra)
  library(xcms)
  library(BiocParallel)
  library(vioplot)
  library(matrixStats)
})

get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
    return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
  }
  getwd()
}

proj_root <- get_project_root()
res_dir <- file.path(proj_root, "repro_results")
proc_dir <- file.path(proj_root, "data", "processed")

if (.Platform$OS.type == "unix") {
  register(MulticoreParam(2))
} else {
  register(SnowParam(2))
}

cat(">>> [1/6] 加载 Step 01 保存的 MsExperiment 对象...\n")
lcms1 <- readMsObject(MsExperiment(), AlabasterParam(path = file.path(proc_dir, "01_lcms1_filtered")))
meta_info <- readRDS(file.path(proc_dir, "01_qc_metadata.rds"))
col_phenotype <- meta_info$col_phenotype
col_sample <- meta_info$col_sample
intern_standard <- meta_info$intern_standard

cat(">>> [2/6] 执行 CentWave 峰检测与裂分峰合并精修...\n")
# CentWave 参数依据 Step 01 内标 EIC 经验测定：HILIC UHPLC 峰宽约 1-8s，高分辨质谱质量精度 15 ppm，
# integrate = 2 基于墨西哥帽小波拟合精确确定峰边界（较 integrate=1 抗拖尾噪声更好）
param_cw <- CentWaveParam(peakwidth = c(1, 8), ppm = 15, integrate = 2)
lcms1 <- findChromPeaks(lcms1, param = param_cw, chunkSize = 2)

# 将同一物质被 CentWave 误拆分的相邻肩峰/裂分峰重新合并（RT 间距 <= 2.5s, m/z <= 0.0015 Da, 谷峰比 >= 0.75）
param_merge <- MergeNeighboringPeaksParam(expandRt = 2.5, expandMz = 0.0015, minProp = 0.75)
lcms1 <- refineChromPeaks(lcms1, param = param_merge, chunkSize = 5)

n_peaks_raw <- nrow(chromPeaks(lcms1))
n_merged <- sum(chromPeakData(lcms1)$merged, na.rm = TRUE)
cat(sprintf("    检测并精修后总色谱峰数: %d (其中合并峰: %d)\n", n_peaks_raw, n_merged))

# 峰检测汇总图 (Fig 9/Peak summary)
pk_list <- split.data.frame(
  chromPeaks(lcms1, columns = c("mzmin", "mzmax", "rtmin", "rtmax")),
  chromPeaks(lcms1, columns = "sample")[, "sample"]
)
pk_list <- lapply(pk_list, function(z) {
  cbind(
    z,
    mz_width_ppm = (z[, "mzmax"] - z[, "mzmin"]) * 1e6 / z[, "mzmax"],
    rt_width = z[, "rtmax"] - z[, "rtmin"]
  )
})

png(file.path(res_dir, "04_peak_detection_summary.png"), width = 1600, height = 1500, res = 200)
par(mfrow = c(3, 1), mar = c(1.5, 4.5, 2, 0.5))
barplot(unlist(lapply(pk_list, nrow)), col = col_sample, ylab = "Peak count",
        main = "Chromatographic Peak Detection Summary (MTBLS8735)", xaxt = "n")
grid()
legend("top", horiz = TRUE, col = col_phenotype, pch = 15, legend = names(col_phenotype), bty = "n")
par(mar = c(1.5, 4.5, 0.5, 0.5))
vioplot(lapply(pk_list, function(z) z[, "mz_width_ppm"]), outline = FALSE,
        ylab = "m/z width [ppm]", xaxt = "n", col = col_sample)
grid()
par(mar = c(4, 4.5, 0.5, 0.5))
vioplot(lapply(pk_list, function(z) z[, "rt_width"]), ylab = "RT width [s]",
        names = sampleData(lcms1)$sample_name, col = col_sample)
grid()
dev.off()

cat(">>> [3/6] 基于 Pooled QC 锚点进行保留时间对齐 (PeakGroups)...\n")
# 先对 QC 样本执行粗分组（minFraction = 0.9 确保锚点在 >=90% 的 QC 中稳定存在），再用 LOESS (span=0.5) 拟合全局 RT 漂移，
# subsetAdjust = "average" 通过相邻 QC 进样序插值推算研究样本的校正偏移
param_init <- PeakDensityParam(
  sampleGroups = sampleData(lcms1)$phenotype == "QC",
  minFraction = 0.9, binSize = 0.01, ppm = 10, bw = 2
)
lcms1 <- groupChromPeaks(lcms1, param = param_init)

subset_qc <- which(sampleData(lcms1)$phenotype == "QC")
param_align <- PeakGroupsParam(
  minFraction = 0.9, extraPeaks = 50, span = 0.5,
  subsetAdjust = "average", subset = subset_qc
)
lcms1 <- adjustRtime(lcms1, param = param_align)

png(file.path(res_dir, "05_rt_alignment_diagnostics.png"), width = 1600, height = 1000, res = 200)
plotAdjustedRtime(lcms1, col = paste0(col_sample, "80"), peakGroupsPch = 1,
                  main = "Retention Time Alignment (QC-subset PeakGroups)")
grid()
legend("topright", col = col_phenotype, legend = names(col_phenotype), lty = 1, bty = "n")
dev.off()

# 将调整后的保留时间永久覆盖回 spectra 与 chromPeaks
lcms1 <- applyAdjustedRtime(lcms1)

cat(">>> [4/6] 执行正式特征分组 (PeakDensity, bw = 1.8)...\n")
# 对齐后带宽 bw 收紧至 1.8s；特征至少在任一生物分组或 QC 组 75% 以上样本中被检测到方形成共识特征
param_grp <- PeakDensityParam(
  sampleGroups = sampleData(lcms1)$phenotype,
  minFraction = 0.75, binSize = 0.01, ppm = 10, bw = 1.8
)
lcms1 <- groupChromPeaks(lcms1, param = param_grp)

na_before <- sum(is.na(featureValues(lcms1))) / length(featureValues(lcms1)) * 100
n_feat_before <- nrow(featureDefinitions(lcms1))
cat(sprintf("    分组后特征数: %d | 缺峰回填前缺失率: %.2f%%\n", n_feat_before, na_before))

cat(">>> [5/6] 执行缺峰回填 (ChromPeakAreaParam) 与缺失率过滤...\n")
# 直接在特征共识 m/z-RT 边界内积分原始质谱信号以填补漏检 NA（避免简单用 0 或全局最小值带来假阳性差异）
lcms1 <- fillChromPeaks(lcms1, param = ChromPeakAreaParam(), chunkSize = 5)
na_after <- sum(is.na(featureValues(lcms1))) / length(featureValues(lcms1)) * 100

# 评估 QC 样本检测峰 vs 回填峰相关性 (Fig 23)
qc_idx <- sampleData(lcms1)$phenotype == "QC"
vals_detect <- featureValues(lcms1, filled = FALSE)[, qc_idx]
vals_filled <- featureValues(lcms1)[, qc_idx]
vals_filled[!is.na(vals_detect)] <- NA
has_filled <- is.na(rowSums(vals_detect))
avg_detect <- rowMeans(vals_detect[has_filled, ], na.rm = TRUE)
avg_filled <- rowMeans(vals_filled[has_filled, ], na.rm = TRUE)
valid_lm <- is.finite(avg_detect) & is.finite(avg_filled) & avg_detect > 0 & avg_filled > 0
fit_gf <- lm(log2(avg_filled[valid_lm]) ~ log2(avg_detect[valid_lm]))

png(file.path(res_dir, "06_gap_filling_qc_correlation.png"), width = 1200, height = 1000, res = 200)
plot(log2(avg_detect[valid_lm]), log2(avg_filled[valid_lm]),
     xlab = expression(log[2]~Detected~QC~Intensity),
     ylab = expression(log[2]~Gap-Filled~QC~Intensity),
     main = sprintf("Gap-Filling Validation (Slope = %.2f, R^2 = %.2f)",
                    coef(fit_gf)[2], summary(fit_gf)$r.squared),
     pch = 21, bg = "#00000025", col = "#00000080")
grid()
abline(0, 1, lty = 2, col = "blue")
abline(fit_gf, col = "red", lwd = 1.5)
legend("topleft", legend = c("y = x", "Linear fit"), col = c("blue", "red"), lty = c(2, 1), bty = "n")
dev.off()

# 过滤条件: QC 样本设 NA 不参与组内检出率考核；特征至少在 CVD 或 CTR 任意一组中有 >60% 样本存在（即缺失 <= 40%）
f <- sampleData(lcms1)$phenotype
f[f == "QC"] <- NA
f <- as.factor(f)
lcms1 <- filterFeatures(lcms1, PercentMissingFilter(f = f, threshold = 40))

n_feat_after <- nrow(featureDefinitions(lcms1))
cat(sprintf("    回填后缺失率: %.2f%% | 缺失过滤后保留特征数: %d\n", na_after, n_feat_after))

cat(">>> [6/6] 导出 SummarizedExperiment 与 XcmsExperiment 供下游使用...\n")
# method = "sum" 对同一特征内多个整合峰采用峰面积求和；同时保存回填前 (raw) 与回填后 (raw_filled) 双 assay
res <- quantify(lcms1, method = "sum", filled = FALSE)
assays(res)$raw_filled <- featureValues(lcms1, method = "sum", filled = TRUE)

lcms_out <- file.path(proc_dir, "02_preprocessed_lcms1")
se_out <- file.path(proc_dir, "02_preprocessed_res")
if (dir.exists(lcms_out)) unlink(lcms_out, recursive = TRUE)
if (dir.exists(se_out)) unlink(se_out, recursive = TRUE)

saveMsObject(lcms1, AlabasterParam(path = lcms_out))
# 清空复杂元数据以免 alabaster.se 序列化过程报 S4 内部环境不支持错
res_export <- res
metadata(res_export) <- list()
saveObject(res_export, se_out)

prep_stats <- data.frame(
  Metric = c("ChromPeaks_After_Refine", "Merged_Peaks", "Features_Grouped",
             "Missing_Pct_Before_GapFill", "Missing_Pct_After_GapFill",
             "GapFill_QC_Slope", "Features_After_MissingFilter"),
  Value = c(n_peaks_raw, n_merged, n_feat_before,
            round(na_before, 2), round(na_after, 2),
            round(coef(fit_gf)[2], 3), n_feat_after)
)
write.table(prep_stats, file.path(res_dir, "02_xcms_preprocessing_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Step 02 完成！预处理统计已保存。\n")
