# ==============================================================================
# 01_metabonaut_download_and_qc.R
# Metabonaut Step 1: 从 MetaboLights (MTBLS8735) 拉取 LC-MS 数据并完成质控可视化
#
# 上游源码对应:
#   upstream_Metabonaut/vignettes/a-end-to-end-untargeted-metabolomics.qmd (Loading data & QC 章节)
# 输入文件:
#   - 远程自动拉取: MetaboLights MTBLS8735 (a_MTBLS8735_LC-MS_positive_hilic_metabolite_profiling.txt)
#   - 本地静态表: data/intern_standard_list.txt (20 种稳定同位素标记内标)
# 输出产物:
#   - repro_results/01_mtbls8735_sample_metadata.tsv
#   - repro_results/01_bpc_raw_and_filtered.png (对应原著 Fig 2-3)
#   - repro_results/02_tic_similarity_heatmap.png (对应原著 Fig 4)
#   - repro_results/03_internal_standards_eic.png (对应原著 Fig 5)
#   - data/processed/01_lcms1_filtered/ & 01_qc_metadata.rds
# ==============================================================================

suppressPackageStartupMessages({
  # 修复要点: 原著 .qmd 隐式借助 BiocStyle/MSnbase 暴露 fData()，独立脚本须显式加载 Biobase
  library(Biobase)
  library(MsExperiment)
  library(MsIO)
  library(MsBackendMetaboLights)
  library(Spectra)
  library(xcms)
  library(BiocParallel)
  library(RColorBrewer)
  library(pheatmap)
})

# 解析 --file= 定位 Metabonaut 项目根目录
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
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

# 跨平台并行后端配置：Windows 缺乏 fork 机制须使用 SnowParam，Unix 可用 MulticoreParam
if (.Platform$OS.type == "unix") {
  register(MulticoreParam(2))
} else {
  register(SnowParam(2))
}

cat(">>> [1/5] 从 MetaboLights (MTBLS8735) 加载 LC-MS1 数据集...\n")
# 通过 MsBackendMetaboLights 按需下载 MTBLS8735 HILIC 正离子模式的 10 个 .mzML 文件并缓存至 BiocFileCache
param <- MetaboLightsParam(
  mtblsId = "MTBLS8735",
  assayName = "a_MTBLS8735_LC-MS_positive_hilic_metabolite_profiling.txt",
  filePattern = ".mzML"
)

lcms1 <- readMsObject(
  MsExperiment(),
  param,
  keepOntology = FALSE,
  keepProtocol = FALSE,
  simplify = TRUE
)

cat(">>> [2/5] 整理样本元数据与分组配色...\n")
# 规范化 ISA-Tab 元数据列名，并将 Pooled QC 样本统一重命名为 POOL1-POOL4 便于区分进样次序
colnames(sampleData(lcms1)) <- c(
  "sample_name", "derived_spectra_data_file", "metabolite_asssignment_file",
  "source_name", "organism", "blood_sample_type",
  "sample_type", "age", "unit", "phenotype"
)
sampleData(lcms1)$phenotype[sampleData(lcms1)$sample_name == "POOL"] <- "QC"
sampleData(lcms1)$sample_name[sampleData(lcms1)$sample_name == "POOL"] <- c("POOL1", "POOL2", "POOL3", "POOL4")
# 记录连续进样顺序 (用于后续 Step 04 notame 批次内信号漂移校正)
sampleData(lcms1)$injection_index <- seq_len(nrow(sampleData(lcms1)))

# 统一全流程配色：QC 灰、CVD 紫、CTR 蓝（与原著保持一致）
col_phenotype <- brewer.pal(9, name = "Set1")[c(9, 5, 4)]
names(col_phenotype) <- c("QC", "CVD", "CTR")
col_sample <- col_phenotype[sampleData(lcms1)$phenotype]

write.table(
  as.data.frame(sampleData(lcms1)[, c("sample_name", "phenotype", "age", "injection_index", "derived_spectra_data_file")]),
  file.path(res_dir, "01_mtbls8735_sample_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat(sprintf("    已加载 %d 个样本 | 总谱图数: %d\n", nrow(sampleData(lcms1)), length(spectra(lcms1))))

cat(">>> [3/5] 绘制 BPC 过滤前后对比图 (保留 RT 10-240s)...\n")
# 提取基峰色谱图 (BPC, 每个扫描点取最大强度)；HILIC 柱前 10s 为死时间盐分洗脱区，240s 后为柱平衡冲洗段
bpc_raw <- chromatogram(lcms1, aggregationFun = "max")
lcms1 <- filterRt(lcms1, c(10, 240))
bpc_filt <- chromatogram(lcms1, aggregationFun = "max")

png(file.path(res_dir, "01_bpc_raw_and_filtered.png"), width = 1800, height = 1400, res = 200)
par(mfrow = c(2, 1), mar = c(4, 4.5, 2.5, 1))
plot(bpc_raw, col = paste0(col_sample, "80"), main = "BPC (Raw MTBLS8735 LC-MS)", lwd = 1.5)
grid()
abline(v = c(10, 240), lty = 2, col = "red")
legend("topright", col = col_phenotype, legend = names(col_phenotype), lty = 1, lwd = 2, horiz = TRUE, bty = "n")

plot(bpc_filt, col = paste0(col_sample, "80"), main = "BPC (After RT Filter: 10-240 s)", lwd = 1.5)
grid()
legend("topright", col = col_phenotype, legend = names(col_phenotype), lty = 1, lwd = 2, horiz = TRUE, bty = "n")
dev.off()

cat(">>> [4/5] 绘制样本间 TIC Pearson 相关性热图...\n")
# 按 2 秒为箱对总离子流色谱 (TIC) 分箱后计算样本间相关性，快速排查进样失败或异常离群样
tic <- chromatogram(lcms1, aggregationFun = "sum") |> bin(binSize = 2)
ticmap <- do.call(cbind, lapply(tic, intensity)) |> cor()
rownames(ticmap) <- colnames(ticmap) <- sampleData(lcms1)$sample_name
ann <- data.frame(phenotype = sampleData(lcms1)[, "phenotype"])
rownames(ann) <- rownames(ticmap)

png(file.path(res_dir, "02_tic_similarity_heatmap.png"), width = 1400, height = 1200, res = 200)
pheatmap(
  ticmap,
  annotation_col = ann,
  annotation_colors = list(phenotype = col_phenotype),
  main = "Sample TIC Similarity (Pearson Correlation, 2s bins)"
)
dev.off()

cat(">>> [5/5] 提取并绘制同位素内标 EIC...\n")
is_path <- file.path(proj_root, "data", "intern_standard_list.txt")
intern_standard <- read.delim(is_path)

# 利用预设的 m/z ± 0.005 Da 与 RT ± 10s 窗口抽取 20 种稳定同位素标记内标的 EIC，
# 目的是在找峰前经验性核查色谱峰宽度（确定后续 CentWave peakwidth = c(1, 8)）及保留时间偏移幅度
eic_is <- chromatogram(
  lcms1,
  rt = as.matrix(intern_standard[, c("rtmin", "rtmax")]),
  mz = as.matrix(intern_standard[, c("mzmin", "mzmax")])
)
fData(eic_is)$mz <- intern_standard$mz
fData(eic_is)$rt <- intern_standard$RT
fData(eic_is)$name <- intern_standard$name
fData(eic_is)$abbreviation <- intern_standard$abbreviation
rownames(fData(eic_is)) <- intern_standard$abbreviation

# 选取峰形评级为 0_perfect / 1_good 的 L-Cystine 与 L-Methionine 进行典型展示
eic_cystine <- eic_is["cystine_13C_15N"]
eic_met <- eic_is["methionine_13C_15N"]

png(file.path(res_dir, "03_internal_standards_eic.png"), width = 1800, height = 850, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
plot(eic_cystine, main = fData(eic_cystine)$name, cex.axis = 0.85, cex.main = 0.9, col = paste0(col_sample, "80"))
grid()
abline(v = fData(eic_cystine)$rt, col = "red", lty = 3)

plot(eic_met, main = fData(eic_met)$name, cex.axis = 0.85, cex.main = 0.9, col = paste0(col_sample, "80"))
grid()
abline(v = fData(eic_met)$rt, col = "red", lty = 3)
legend("topright", col = col_phenotype, legend = names(col_phenotype), lty = 1, bty = "n")
dev.off()

# 使用 Alabaster 跨语言规范保存过滤后的 MsExperiment 对象供 Step 02 加载
obj_dir <- file.path(proc_dir, "01_lcms1_filtered")
if (dir.exists(obj_dir)) unlink(obj_dir, recursive = TRUE)
saveMsObject(lcms1, AlabasterParam(path = obj_dir))
saveRDS(list(col_phenotype = col_phenotype, col_sample = col_sample, intern_standard = intern_standard),
        file.path(proc_dir, "01_qc_metadata.rds"))

cat("Step 01 完成！已保存:", obj_dir, "\n")
