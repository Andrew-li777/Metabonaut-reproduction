# ==============================================================================
# 04_notame_drift_pqn_selection.R
# Metabonaut Step 4: 基于 notame 生态的漂移校正、PQN 归一化与特征筛选
#
# 上游源码对应:
#   upstream_Metabonaut/vignettes/notame_normalization_and_feature_selection.qmd
# 输入文件:
#   - 优先加载 Step 02 生成的 data/processed/02_preprocessed_res/ (SummarizedExperiment)
#   - 回退选项: upstream_Metabonaut/inst/extdata/preprocessed_res/
# 输出产物:
#   - repro_results/11_notame_injection_lm.png (进样漂移线性 P 值分布对比)
#   - repro_results/12_notame_boxplots_and_pca.png (样品箱线图与 PCA 对比)
#   - repro_results/13_notame_beeswarm_and_volcano.png (FT0845 蜂群图与火山图)
#   - repro_results/04_notame_flag_report.tsv & 04_notame_top_lm_features.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(alabaster.se)
  library(BiocParallel)
  library(cowplot)
  library(ggplot2)
  library(limma)
  library(notame)
  library(notameViz)
  library(notameStats)
  library(SummarizedExperiment)
})

set.seed(2024)

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

if (.Platform$OS.type == "unix") {
  register(MulticoreParam(2))
} else {
  register(SnowParam(2))
}

# 优先读取本机 Step 02 现算产物；若尚未运行 Step 02，可回退到上游 extdata 预计算对象
local_se_dir <- file.path(proc_dir, "02_preprocessed_res")
upstream_se_dir <- file.path(proj_root, "upstream_Metabonaut", "inst", "extdata", "preprocessed_res")
se_dir <- if (dir.exists(local_se_dir)) local_se_dir else upstream_se_dir

cat(">>> [1/5] 加载 SummarizedExperiment 对象:", se_dir, "\n")
se <- readObject(se_dir)

# 对齐 notame 要求的 colData 标准列名 (Sample_ID, QC, Injection_order) 与 rowData 标准列 (Feature_ID, Split, Flag)
rename_ind <- which(colnames(colData(se)) %in% c("sample_name", "sample_type", "injection_index"))
colnames(colData(se))[rename_ind] <- c("Sample_ID", "QC", "Injection_order")
colnames(se) <- colData(se)$Sample_ID
se$QC[se$QC == "pool"] <- "QC"
se$phenotype <- as.factor(se$phenotype)

rowData(se) <- cbind(Feature_ID = rownames(se), rowData(se))
rowData(se)$Split <- "HILIC_pos"
rowData(se)$Flag <- NA

# 将 alabaster DelayedMatrix 转为内存 matrix 兼容 missForest 随机森林填补与 notame 内部矩阵运算
assay(se, "raw_filled") <- as.matrix(assay(se, "raw_filled"))
assay(se, "raw") <- as.matrix(assay(se, "raw"))

cat(">>> [2/5] 执行检出率过滤、三次样条漂移校正与 PQN 归一化...\n")
# 检出率标记：在 QC 中检出率 >= 70% 或在任意生物组内检出率 >= 80% 方予通过
se <- flag_detection(se, qc_limit = 0.7, group_limit = 0.8,
                     group = "phenotype", assay.type = "raw_filled")
# 基于进样序列与 Pooled QC 样本拟合平滑三次样条曲线，去除随进样时间的信号衰减/漂移
se <- correct_drift(se, assay.type = "raw_filled", name = "drift_norm")
# 概率商归一化 (PQN)：以 Pooled QC 中位数谱为参考消除样本间总浓度/稀释差异
se <- pqn_normalization(se, ref = "qc", method = "median",
                        assay.type = "drift_norm", name = "dil_norm")

cat(">>> [3/5] 绘制 notame 归一化前后诊断对比图...\n")
p_inj <- plot_grid(
  plot_injection_lm(se, assay.type = "raw_filled"),
  plot_injection_lm(se, assay.type = "dil_norm"),
  labels = "AUTO"
)
ggsave(file.path(res_dir, "11_notame_injection_lm.png"), p_inj, width = 14, height = 6.5, dpi = 180)

p_box_pca <- plot_grid(
  plot_sample_boxplots(se, order_by = "Injection_order", fill_by = "phenotype", assay.type = "raw_filled"),
  plot_sample_boxplots(se, order_by = "Injection_order", fill_by = "phenotype", assay.type = "dil_norm"),
  plot_pca(se, color = "Injection_order", shape = "phenotype", assay.type = "raw_filled"),
  plot_pca(se, color = "Injection_order", shape = "phenotype", assay.type = "dil_norm"),
  ncol = 2, labels = "AUTO"
)
ggsave(file.path(res_dir, "12_notame_boxplots_and_pca.png"), p_box_pca, width = 13, height = 9, dpi = 180)

cat(">>> [4/5] 质量标记 (RSD_r < 0.2 & D_ratio_r < 0.4) 与随机森林填补...\n")
# 稳健相对标准偏差 RSD_r < 0.2 且稳健分散比 D_ratio_r < 0.4；未达标特征仅在 Flag 列做标记而先保留供诊断
se <- flag_quality(se, assay.type = "dil_norm", condition = "RSD_r < 0.2 & D_ratio_r < 0.4")
f_report <- flag_report(se)
write.table(f_report, file.path(res_dir, "04_notame_flag_report.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 剔除 QC 样本后用 missForest 对残余缺失值做非线性随机森林填补，再清理带质量标记的特征
set.seed(2024)
pretreated <- impute_rf(drop_qcs(se), assay.type = "dil_norm", name = "imputed")
base <- drop_flagged(pretreated)
assays(base)[names(assays(base)) != "imputed"] <- NULL

cat(sprintf("    质控过滤后保留特征数: %d / %d\n", nrow(base), nrow(se)))

cat(">>> [5/5] 线性模型差异检验 (Feature ~ phenotype + age) 与关键代谢物可视化...\n")
# 对数转换后逐特征拟合普通最小二乘线性模型 (控制 age 协变量) 并做 BH-FDR 多重校验
assay(base, "log2") <- log2(assay(base))
lm_results <- perform_lm(base, assay.type = "log2", formula_char = "Feature ~ phenotype + age")
base <- join_rowData(base, lm_results[, c("Feature_ID", "phenotypeCVD.p.value",
                                          "phenotypeCVD.p.value_FDR", "phenotypeCVD.estimate")])

# 导出 Top 15 最低 p 值特征
ord <- order(rowData(base)$phenotypeCVD.p.value)
top_notame <- as.data.frame(rowData(base)[ord[seq_len(min(15, nrow(base)))],
  c("Feature_ID", "mzmed", "rtmed", "phenotypeCVD.estimate", "phenotypeCVD.p.value", "phenotypeCVD.p.value_FDR")])
write.table(top_notame, file.path(res_dir, "04_notame_top_lm_features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 定位咖啡因或最小 p 值特征 (原著为 FT0845)
target_ft <- if ("FT0845" %in% rownames(base)) "FT0845" else rownames(base)[ord[1]]

p_bee_a <- save_beeswarm_plots(pretreated[target_ft, ], x = "phenotype", color = "phenotype",
                               save = FALSE, assay.type = "raw_filled")[[1]]
p_bee_b <- save_beeswarm_plots(base[target_ft, ], x = "phenotype", color = "phenotype",
                               save = FALSE, assay.type = "imputed")[[1]]
p_volc <- volcano_plot(base, x = "phenotypeCVD.estimate", p = "phenotypeCVD.p.value") +
  geom_text(data = as.data.frame(rowData(base)[target_ft, , drop = FALSE]),
            aes(label = target_ft), vjust = -1.5, color = "#b2182b")

p_final <- plot_grid(p_bee_a, p_bee_b, p_volc, ncol = 3, labels = "AUTO")
ggsave(file.path(res_dir, "13_notame_beeswarm_and_volcano.png"), p_final, width = 15, height = 5, dpi = 180)

cat("Step 04 (notame) 完成！核心特征:", target_ft, "\n")
