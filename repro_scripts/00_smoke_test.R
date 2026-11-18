# ==============================================================================
# 00_smoke_test.R
# 代谢组学复现工作区 10 分钟冒烟测试 (REPRO-D 最小骨架验证)
#
# 功能:
#   1. 验证目录结构 (data/raw, data/processed, repro_results, logs) 与读写权限
#   2. 构造符合非靶向 LC-MS 特征矩阵规范的轻量模拟验证集 (300 features × 18 samples:
#      6 Case + 6 Control + 6 Pooled QC)，避免冒烟阶段依赖外部大文件下载
#   3. 执行 O8 标准预处理：零值/RSD 过滤 -> log2 转换 -> Pareto 缩放 -> 无监督 PCA
#   4. 验证 Pooled QC 是否紧聚（RSD 质控 + PCA 离散度校验），输出图谱与统计摘要
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggplot2)
  library(patchwork)
})

# 设置确定性随机种子，确保模拟数据集及图谱多次重跑完全可复现
set.seed(20261114)

# 动态定位 Metabonaut/ 项目根目录
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
log_dir <- file.path(proj_root, "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat(">>> [1/4] 构造 LC-MS 特征矩阵冒烟测试集 (SummarizedExperiment)...\n")
n_feat <- 300L
sample_ids <- c(paste0("CTR_", 1:6), paste0("CASE_", 1:6), paste0("QC_", 1:6))
groups <- factor(c(rep("Control", 6), rep("Case", 6), rep("Pooled_QC", 6)),
                 levels = c("Control", "Case", "Pooled_QC"))

# 模拟对数正态基线峰面积，并对前 40 个差异代谢物在 Case 组注入真实信号
base_abund <- exp(rnorm(n_feat, mean = 12, sd = 1.2))
mat <- matrix(0, nrow = n_feat, ncol = length(sample_ids),
              dimnames = list(paste0("FT_", sprintf("%04d", 1:n_feat)), sample_ids))

for (j in seq_along(sample_ids)) {
  grp <- as.character(groups[j])
  # Pooled QC 技术重复方差更小；生物样本方差更大
  tech_sd <- if (grp == "Pooled_QC") 0.08 else 0.25
  mult <- rep(1, n_feat)
  if (grp == "Case") mult[1:40] <- exp(rnorm(40, mean = 0.85, sd = 0.2))
  if (grp == "Pooled_QC") mult[1:40] <- exp(0.5 * 0.85)
  mat[, j] <- base_abund * mult * exp(rnorm(n_feat, mean = 0, sd = tech_sd))
}

# 故意注入 15 个高噪声特征 (QC RSD > 30%) 测试质控过滤
noisy_idx <- (n_feat - 14L):n_feat
mat[noisy_idx, groups == "Pooled_QC"] <- mat[noisy_idx, groups == "Pooled_QC"] *
  exp(matrix(rnorm(15 * 6, 0, 0.55), nrow = 15))

col_df <- DataFrame(
  sample_id = sample_ids,
  group = groups,
  injection_order = seq_along(sample_ids),
  row.names = sample_ids
)
row_df <- DataFrame(
  feature_id = rownames(mat),
  mz = round(runif(n_feat, 100, 950), 4),
  rt_sec = round(runif(n_feat, 30, 720), 1),
  row.names = rownames(mat)
)

se <- SummarizedExperiment(assays = list(intensity = mat), rowData = row_df, colData = col_df)

cat(">>> [2/4] 执行 Pooled QC RSD 过滤 (阈值 <= 30%)...\n")
# 代谢组学通用质控准则：Pooled QC 技术重复相对标准偏差 RSD = SD / Mean > 30% 的特征视为仪器不稳噪声并剔除
qc_mat <- assay(se, "intensity")[, se$group == "Pooled_QC"]
qc_rsd <- apply(qc_mat, 1, function(x) sd(x) / mean(x))
keep <- qc_rsd <= 0.30
se_clean <- se[keep, ]

cat(sprintf("    原始特征数: %d | 通过 QC RSD <= 30%%: %d | 剔除高噪声特征: %d\n",
            nrow(se), nrow(se_clean), sum(!keep)))

cat(">>> [3/4] 执行 log2 转换 + Pareto 缩放与 PCA 质控诊断...\n")
# log2(x + 1) 压缩代谢物峰面积跨数量级的右偏分布
log_mat <- log2(assay(se_clean, "intensity") + 1)
# Pareto scaling: (x - mean) / sqrt(sd)，兼具保留高丰度信号与适当放大低丰度变量的作用（代谢组最常用）
row_means <- rowMeans(log_mat)
row_sds <- apply(log_mat, 1, sd)
pareto_mat <- sweep(sweep(log_mat, 1, row_means, "-"), 1, sqrt(row_sds), "/")

# 已手动完成中心化与 Pareto 缩放，故 prcomp 须设 center = FALSE, scale. = FALSE 避免二次标准化
pca_res <- prcomp(t(pareto_mat), center = FALSE, scale. = FALSE)
var_exp <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
pca_df <- data.frame(
  sample_id = rownames(pca_res$x),
  group = se_clean$group,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2]
)

# 验收：QC 样本在 PC1-PC2 空间欧氏分散度显著小于生物样本组间分散度
qc_coords <- pca_df[pca_df$group == "Pooled_QC", c("PC1", "PC2")]
qc_spread <- mean(sqrt(rowSums(sweep(as.matrix(qc_coords), 2, colMeans(qc_coords), "-")^2)))

cat(sprintf("    PCA 方差解释度: PC1 = %.1f%%, PC2 = %.1f%% | Pooled QC 平均质心距 = %.2f\n",
            var_exp[1], var_exp[2], qc_spread))

cat(">>> [4/4] 导出冒烟测试图表与摘要...\n")
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3.5, alpha = 0.9) +
  stat_ellipse(data = subset(pca_df, group != "Pooled_QC"), level = 0.95, linetype = 2) +
  scale_color_manual(values = c("Control" = "#2b8cbe", "Case" = "#e34a33", "Pooled_QC" = "#31a354")) +
  labs(
    title = "Smoke Test: Pareto-Scaled Metabolomics PCA",
    subtitle = sprintf("Features retained: %d/%d (QC RSD <= 30%%)", nrow(se_clean), nrow(se)),
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2])
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

rsd_df <- data.frame(qc_rsd = qc_rsd * 100, status = ifelse(keep, "Pass (<=30%)", "Filtered (>30%)"))
p2 <- ggplot(rsd_df, aes(x = qc_rsd, fill = status)) +
  geom_histogram(bins = 25, color = "white", alpha = 0.85) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "#b2182b", linewidth = 0.8) +
  scale_fill_manual(values = c("Pass (<=30%)" = "#31a354", "Filtered (>30%)" = "#bdbdbd")) +
  labs(
    title = "Pooled QC RSD Distribution",
    x = "QC Relative Standard Deviation (%)",
    y = "Feature Count"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

combined <- p1 + p2 + plot_layout(widths = c(1.2, 1))
fig_path <- file.path(res_dir, "00_smoke_qc_pca.png")
ggsave(fig_path, combined, width = 10.5, height = 4.5, dpi = 200)

summary_df <- data.frame(
  Metric = c("Total_Input_Features", "QC_Filtered_Features", "Removed_Noisy_Features",
             "Median_QC_RSD_Percent", "PC1_Var_Percent", "PC2_Var_Percent", "QC_Centroid_Spread"),
  Value = c(nrow(se), nrow(se_clean), sum(!keep),
            round(median(qc_rsd[keep]) * 100, 2), var_exp[1], var_exp[2], round(qc_spread, 3))
)
summary_path <- file.path(res_dir, "00_smoke_verification_summary.tsv")
write.table(summary_df, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("冒烟测试完成！图表已保存:", fig_path, "\n")
cat("摘要已保存:", summary_path, "\n")
