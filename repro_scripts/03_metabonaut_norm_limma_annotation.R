# ==============================================================================
# 03_metabonaut_norm_limma_annotation.R
# Metabonaut Step 3: 归一化、D-ratio 质控过滤、limma 差异检验与 MS1/MS2 注释
#
# 上游源码对应:
#   upstream_Metabonaut/vignettes/a-end-to-end-untargeted-metabolomics.qmd
#   (Normalization / Filtering / Differential Abundance / Annotation 章节)
# 输入文件:
#   - data/processed/02_preprocessed_lcms1/ (XcmsExperiment)
#   - data/processed/02_preprocessed_res/ (SummarizedExperiment)
#   - data/processed/01_qc_metadata.rds
#   - 远程按需缓存: AnnotationHub AH116166 (MassBank) 与 MTBLS8735 LC-MS/MS .mzML
# 输出产物:
#   - repro_results/07_rla_before_after_norm.png (对应原著 Fig 24)
#   - repro_results/08_study_samples_pca.png (对应原著 Fig 25)
#   - repro_results/09_limma_pvals_and_volcano.png (对应原著 Fig 26-27)
#   - repro_results/10_caffeine_eic_and_ms2_mirror.png (对应原著 Fig 35-37)
#   - repro_results/03_limma_significant_features.tsv
#   - repro_results/03_ms1_massbank_annotations.tsv & 03_ms2_massbank_annotations.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(MsExperiment)
  library(MsIO)
  library(alabaster.se)
  library(MsBackendMetaboLights)
  library(Spectra)
  library(xcms)
  library(MetaboCoreUtils)
  library(limma)
  library(matrixStats)
  library(AnnotationHub)
  library(CompoundDb)
  library(MetaboAnnotation)
  library(BiocParallel)
  library(ggfortify)
  library(gridExtra)
  library(RColorBrewer)
})

set.seed(20261114)

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

cat(">>> [1/5] 加载 Step 02 预处理对象...\n")
# 修复要点 1: readMsObject 必须显式传入 XcmsExperiment() 原型，否则将丢失 chromPeaks 与特征定义
lcms1 <- readMsObject(XcmsExperiment(), AlabasterParam(path = file.path(proc_dir, "02_preprocessed_lcms1")))
res <- readObject(file.path(proc_dir, "02_preprocessed_res"))
# 修复要点 2: alabaster.se 读回的矩阵默认为 DelayedMatrix，直接 sweep("/", ...) 或传给 limma/xcms 过滤会报错，
# 故立即转为标准内存 matrix
for (aname in assayNames(res)) {
  assay(res, aname) <- as.matrix(assay(res, aname))
}
meta_info <- readRDS(file.path(proc_dir, "01_qc_metadata.rds"))
col_phenotype <- meta_info$col_phenotype
col_sample <- col_phenotype[res$phenotype]

cat(">>> [2/5] 缺失值低值均匀填补、中位数缩放与 RLA/PCA 诊断...\n")
# 针对回填后仍残留的 NA (代表真实低丰度或未检出)，按该特征最小观测值的 [min/2, min] 均匀分布随机低值填补
na_unidis <- function(z) {
  na <- is.na(z)
  if (any(na)) {
    mn <- min(z, na.rm = TRUE)
    z[na] <- runif(sum(na), min = mn / 2, max = mn)
  }
  z
}
assays(res)$raw_filled_imputed <- t(apply(assay(res, "raw_filled"), MARGIN = 1, na_unidis))

# 样本内中位数缩放 (Median scaling)：消除样本上样量/雾化效率系统偏差，使其围绕全样本中位数基线对齐
mdns <- apply(assay(res, "raw_filled"), MARGIN = 2, median, na.rm = TRUE)
nf_mdn <- mdns / median(mdns)
assays(res)$norm <- sweep(assay(res, "raw_filled"), MARGIN = 2, nf_mdn, "/")
assays(res)$norm_imputed <- sweep(assay(res, "raw_filled_imputed"), MARGIN = 2, nf_mdn, "/")

# RLA 对比图
png(file.path(res_dir, "07_rla_before_after_norm.png"), width = 1500, height = 1300, res = 200)
par(mfrow = c(2, 1), mar = c(3.5, 4.5, 2.5, 1))
boxplot(rowRla(assay(res, "raw_filled"), group = res$phenotype),
        cex = 0.5, pch = 16, ylab = "RLA", border = col_sample,
        col = paste0(col_sample, "80"), outline = FALSE, xaxt = "n", main = "Raw Data RLA")
grid(nx = NA, ny = NULL)
abline(h = 0, lty = 3)
legend("topright", col = col_phenotype, legend = names(col_phenotype), lty = 1, lwd = 2, horiz = TRUE, bty = "n")

boxplot(rowRla(assay(res, "norm"), group = res$phenotype),
        cex = 0.5, pch = 16, ylab = "RLA", border = col_sample,
        col = paste0(col_sample, "80"), outline = FALSE, xaxt = "n", main = "Median-Normalized RLA")
axis(side = 1, at = seq_len(ncol(res)), labels = res$sample_name)
grid(nx = NA, ny = NULL)
abline(h = 0, lty = 3)
dev.off()

# D-ratio 过滤 (阈值 0.4): D-ratio = MAD(QC) / MAD(Study)，剔除技术波动相对生物变异过大的特征
n_before_dratio <- nrow(res)
filter_dratio <- DratioFilter(
  threshold = 0.4,
  qcIndex = res$phenotype == "QC",
  studyIndex = res$phenotype != "QC",
  mad = TRUE
)
res <- filterFeatures(res, filter = filter_dratio, assay = "norm_imputed")
n_after_dratio <- nrow(res)
cat(sprintf("    D-ratio <= 0.4 过滤: %d -> %d 个高质量特征\n", n_before_dratio, n_after_dratio))

# 剥离 QC 样本用于统计 QC CV，研究样本 (3 CVD vs 3 CTR) 进入下游 PCA 与 limma
res_qc <- res[, res$phenotype == "QC"]
res <- res[, res$phenotype != "QC"]
rowData(res)$qc_cv <- assay(res_qc, "norm") |> rowRsd()

# 过滤后样本 PCA
vals <- assay(res, "norm_imputed") |> t() |> log2() |> scale(center = TRUE, scale = TRUE)
pca_res <- prcomp(vals, scale = FALSE, center = FALSE)
vals_df <- as.data.frame(vals)
vals_df$phenotype <- res$phenotype
vals_df$age <- res$age
p_pca1 <- autoplot(pca_res, data = vals_df, colour = "phenotype", scale = 0) +
  scale_color_manual(values = col_phenotype) + theme_bw() + ggtitle("Study Samples PCA (by Phenotype)")
p_pca2 <- autoplot(pca_res, data = vals_df, colour = "age", scale = 0) +
  theme_bw() + ggtitle("Study Samples PCA (by Age)")
png(file.path(res_dir, "08_study_samples_pca.png"), width = 1600, height = 750, res = 200)
grid.arrange(p_pca1, p_pca2, ncol = 2)
dev.off()

cat(">>> [3/5] 执行 limma 差异丰度检验 (~ phenotype + age)...\n")
# 使用经验贝叶斯线性模型借用全特征方差信息提升小样本 (n=3/组) 检验效能，同时将年龄纳入协变量控制混杂
age <- res$age
phenotype <- factor(res$phenotype, levels = c("CTR", "CVD"))
design <- model.matrix(~ phenotype + age)
fit <- lmFit(log2(assay(res, "norm_imputed")), design = design)
fit <- eBayes(fit)

tmp <- data.frame(
  feature_id = rownames(res),
  coef.CVD = fit$coefficients[, "phenotypeCVD"],
  pvalue.CVD = fit$p.value[, "phenotypeCVD"],
  adjp.CVD = p.adjust(fit$p.value[, "phenotypeCVD"], method = "BH"),
  avg.CVD = rowMeans(log2(assay(res, "norm_imputed")[, res$phenotype == "CVD"])),
  avg.CTR = rowMeans(log2(assay(res, "norm_imputed")[, res$phenotype == "CTR"]))
)
tmp$significant.CVD <- tmp$adjp.CVD < 0.05
rowData(res) <- cbind(rowData(res), tmp)

sig_tab <- as.data.frame(rowData(res)[rowData(res)$significant.CVD,
  c("feature_id", "mzmed", "rtmed", "coef.CVD", "pvalue.CVD", "adjp.CVD", "avg.CTR", "avg.CVD", "qc_cv")])
sig_tab <- sig_tab[order(abs(sig_tab$coef.CVD), decreasing = TRUE), ]
write.table(sig_tab, file.path(res_dir, "03_limma_significant_features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("    limma BH-FDR < 0.05 显著差异特征数: %d\n", nrow(sig_tab)))

png(file.path(res_dir, "09_limma_pvals_and_volcano.png"), width = 1800, height = 650, res = 200)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
hist(rowData(res)$pvalue.CVD, breaks = 40, xlab = "Raw p-value", main = "Raw p-value Distribution")
hist(rowData(res)$adjp.CVD, breaks = 40, xlab = "BH Adjusted p-value", main = "FDR Distribution")
plot(rowData(res)$coef.CVD, -log10(rowData(res)$adjp.CVD),
     xlab = expression(log[2]~Fold~Change~(CVD/CTR)),
     ylab = expression(-log[10]~p[BH]),
     pch = 16, col = "#00000060", main = "Volcano Plot (limma)")
grid()
abline(h = -log10(0.05), col = "blue", lty = 2)
if (any(rowData(res)$significant.CVD)) {
  points(rowData(res)$coef.CVD[rowData(res)$significant.CVD],
         -log10(rowData(res)$adjp.CVD[rowData(res)$significant.CVD]),
         col = "blue", pch = 16, cex = 1.2)
}
dev.off()

cat(">>> [4/5] 基于 MassBank (AnnotationHub) 执行 MS1 精确质量注释...\n")
# 对全部 BH-FDR < 0.05 显著特征按正离子常见加合离子 ([M+H]+, [M+Na]+, [M+H-NH3]+) 与 5 ppm 容差进行 MS1 粗匹配
res_sig <- res[rowData(res)$significant.CVD, ]
ah <- AnnotationHub()
mb_q <- query(ah, "MassBank")
mb_id <- if ("AH116166" %in% names(mb_q)) "AH116166" else tail(names(mb_q), 1)
cat("    使用 AnnotationHub MassBank 记录:", mb_id, "\n")
mb <- ah[[mb_id]]
cmps <- compounds(mb, columns = c("name", "formula", "exactmass", "inchikey"))

param_ms1 <- Mass2MzParam(adducts = c("[M+H]+", "[M+Na]+", "[M+H-NH3]+"), tolerance = 0, ppm = 5)
mtch_ms1 <- matchValues(res_sig, cmps, param = param_ms1, mzColname = "mzmed")
ms1_df <- matchedData(mtch_ms1, c("feature_id", "mzmed", "rtmed", "adduct", "ppm_error",
                                  "target_formula", "target_name", "target_inchikey")) |>
  as.data.frame()
ms1_df <- ms1_df[!is.na(ms1_df$ppm_error), ]
if (nrow(ms1_df) > 0) {
  # 同一特征对应同一化合物多个参考条目时仅保留 ppm 误差最小者
  ms1_df <- split(ms1_df, paste(ms1_df$feature_id, ms1_df$target_inchikey)) |>
    lapply(function(z) z[which.min(z$ppm_error), ]) |>
    do.call(what = rbind)
}
write.table(ms1_df, file.path(res_dir, "03_ms1_massbank_annotations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("    MS1 匹配到 %d 条唯一 (feature, InChIKey) 候选注释\n", nrow(ms1_df)))

cat(">>> [5/5] 拉取 MTBLS8735 LC-MS/MS 数据并执行 MS2 碎片谱图匹配...\n")
# 单独拉取 MTBLS8735 的 MS/MS 实验序列，仅筛选对照组含二级碎裂的运行并裁切至相同 [10, 240] s 色谱区间
param_ms2 <- MetaboLightsParam(
  mtblsId = "MTBLS8735",
  assayName = "a_MTBLS8735_LC-MSMS_positive_hilic_metabolite_profiling.txt",
  filePattern = ".mzML"
)
lcms2 <- readMsObject(MsExperiment(), param_ms2, keepOntology = FALSE, keepProtocol = FALSE, simplify = TRUE)
colnames(sampleData(lcms2)) <- c(
  "sample_name", "derived_spectra_data_file", "metabolite_asssignment_file",
  "source_name", "organism", "blood_sample_type",
  "sample_type", "age", "unit", "phenotype"
)
keep_ms2 <- sampleData(lcms2)$phenotype == "CTR" & grepl("MSMS", sampleData(lcms2)$derived_spectra_data_file)
lcms2 <- lcms2[keep_ms2]
lcms2 <- filterRt(lcms2, c(10, 240))

# MS2 去噪标准步骤：过滤低于基峰 5% 的碎片峰、剔除 >= 母离子 m/z 的未碎裂峰、仅保留含 >1 个碎片的谱图并归一化强度
low_int <- function(x, ...) x > max(x, na.rm = TRUE) * 0.05
ms2_ctr <- spectra(lcms2) |>
  filterMsLevel(2L) |>
  filterIntensity(intensity = low_int) |>
  filterPrecursorPeaks(ppm = 50, mz = ">=")
ms2_ctr <- ms2_ctr[lengths(ms2_ctr) > 1] |> scalePeaks()
ms2_ctr <- setBackend(ms2_ctr, MsBackendMemory()) |> applyProcessing()

# 将 MS1 显著特征的 [rtmin, rtmax] x [mzmin, mzmax] 边界映射回 MS2 谱图集合
target_area <- featureArea(lcms1)[rownames(res_sig), , drop = FALSE]
ms2_ctr_fts <- apply(target_area[, c("rtmin", "rtmax", "mzmin", "mzmax"), drop = FALSE],
                     MARGIN = 1, FUN = filterRanges, object = ms2_ctr,
                     spectraVariables = c("rtime", "precursorMz"))
l_fts <- lengths(ms2_ctr_fts)
ms2_ctr_fts <- concatenateSpectra(ms2_ctr_fts)
ms2_ctr_fts$feature_id <- rep(rownames(res_sig), l_fts)

# 从 MassBank 提取正离子模式参考碎片谱，采用母离子严格匹配 + 归一化点积余弦相似度 >= 0.7 阈值鉴定
ms2_ref <- Spectra(mb) |>
  filterPolarity(1L) |>
  filterIntensity(intensity = low_int) |>
  filterPrecursorPeaks(ppm = 50, mz = ">=")
ms2_ref <- ms2_ref[lengths(ms2_ref) > 1] |> scalePeaks()

# Windows 下 SerialParam 执行光谱对比更稳健，避免多进程 SQLite 句柄序列化开销
register(SerialParam())
prm_cmp <- CompareSpectraParam(ppm = 40, tolerance = 0.05, requirePrecursor = TRUE,
                               THRESHFUN = function(x) which(x >= 0.7))
ms2_mtch <- matchSpectra(ms2_ctr_fts, ms2_ref, param = prm_cmp)
ms2_mtch <- ms2_mtch[whichQuery(ms2_mtch)]
ms2_res <- matchedData(ms2_mtch) |> as.data.frame()

if (nrow(ms2_res) > 0) {
  ms2_res <- split(ms2_res, paste(ms2_res$feature_id, ms2_res$target_inchikey)) |>
    lapply(function(z) z[which.max(z$score), ]) |>
    do.call(what = rbind)
  write.table(ms2_res[, c("feature_id", "precursorMz", "rtime", "target_name", "score", "target_inchikey")],
              file.path(res_dir, "03_ms2_massbank_annotations.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("    MS2 高置信匹配结果:\n")
  print(ms2_res[, c("feature_id", "target_name", "score", "target_inchikey")], row.names = FALSE)

  # 绘制 Top MS2 匹配特征 (Caffeine) EIC + MS2 Mirror 图
  idx_best <- which.max(ms2_res$score)
  best_ft <- ms2_res$feature_id[idx_best]
  q_idx <- which(query(ms2_mtch)$.original_query_index == ms2_res$.original_query_index[idx_best])
  t_idx <- which(target(ms2_mtch)$spectrum_id == ms2_res$target_spectrum_id[idx_best])
  q_spec <- query(ms2_mtch)[q_idx]
  t_spec <- target(ms2_mtch)[t_idx]

  col_all <- col_phenotype[sampleData(lcms1)$phenotype]
  eic_best <- featureChromatograms(lcms1, features = best_ft)

  png(file.path(res_dir, "10_caffeine_eic_and_ms2_mirror.png"), width = 1600, height = 750, res = 200)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  plot(eic_best, col = col_all, peakCol = col_all[chromPeaks(eic_best)[, "sample"]],
       peakBg = paste0(col_all[chromPeaks(eic_best)[, "sample"]], "40"),
       main = sprintf("EIC: %s (%s)", best_ft, ms2_res$target_name[idx_best]))
  abline(v = ms2_res$rtime[idx_best], lty = 2, col = "black")
  legend("topright", col = col_phenotype, legend = names(col_phenotype), lty = 1, bty = "n")

  plotSpectraMirror(q_spec, t_spec,
                    main = sprintf("MS2 Mirror (Cosine Score = %.3f)", ms2_res$score[idx_best]))
  legend("topleft", legend = paste0("Precursor m/z: ", format(precursorMz(q_spec), digits = 6)), bty = "n")
  dev.off()
}

cat("Step 03 完成！\n")
