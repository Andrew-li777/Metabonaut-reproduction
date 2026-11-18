# `Metabonaut` 中间处理对象目录说明 (`data/processed/`)

本目录存放 `Step 01` 与 `Step 02` 生成的序列化中间对象，采用 Bioconductor **Alabaster 跨语言规范**（HDF5 + JSON，确保大质谱对象在不同 R/Bioconductor 版本间无损反序列化）与原生 `.rds` 存储，供下游 `Step 03`（`limma` + MassBank）与 `Step 04`（`notame` 对比管线）直接加载，避免重复耗时峰检测。

> **Git 策略**：为保持代码库轻量，`data/processed/` 下的 Alabaster 目录与 `.rds` 均受 `.gitignore` 规则隔离，仅追踪本说明书 (`README.md`) 与 `.gitkeep`。运行 `Rscript reproduction/metabolomics/Metabonaut/repro_scripts/run_all.R --mode full` 即可完整再生。

---

## 1. 中间产物清单与内部结构

| 目录 / 文件 | 格式 | 实测大小 | 对象类型与关键维度 | 生产脚本 $\rightarrow$ 消费脚本 | 核心内容与使用注意事项 |
|---|---|---:|---|---|---|
| `01_lcms1_filtered/` | Alabaster (`MsIO`) | `0.40 MB` (`422 KB`) | `MsExperiment` (`10` 样本 / 裁切后 `8,250` 张 MS1 谱图，原始 `17,210` 张裁切至 RT `[10, 240] s`) | `Step 01` $\rightarrow$ `Step 02` | 包含原始 `sampleData`（分组、年龄、进样序 `injection_index`）及 `spectra` 后端索引 |
| `01_qc_metadata.rds` | R `.rds` | `1.65 KB` (`1,686 B`) | `list`（键名：`col_phenotype`, `col_sample`, `intern_standard`） | `Step 01` $\rightarrow$ `Step 02`, `Step 03` | 统一全流程 `QC`/`CVD`/`CTR` 调色板与 20 种内标表 |
| `02_preprocessed_lcms1/` | Alabaster (`MsIO`) | `6.86 MB` | `XcmsExperiment`（含回填后共 `108,454` 个 `chromPeaks`，精修阶段为 `89,185` 峰；`8,722` 过滤后特征） | `Step 02` $\rightarrow$ `Step 03` | **关键坑点**：读取时必须传入 `XcmsExperiment()` 原型（见下），否则丢失 `chromPeaks` 与 `featureDefinitions` |
| `02_preprocessed_res/` | Alabaster (`alabaster.se`) | `1.53 MB` | `SummarizedExperiment`（`8,722` 特征 × `10` 样本；双 assay：`raw`, `raw_filled`） | `Step 02` $\rightarrow$ `Step 03`, `Step 04` | **关键坑点**：`readObject()` 返回 `DelayedMatrix`，必须 `as.matrix()` 后再用于 `sweep`/`limma`/`notame` |

---

## 2. 正确反序列化代码示例

```r
library(MsIO); library(MsExperiment); library(xcms); library(alabaster.se); library(SummarizedExperiment)

# 1. 恢复 XcmsExperiment 对象 (必须使用 XcmsExperiment() 原型，不可用 MsExperiment())
lcms1 <- readMsObject(XcmsExperiment(), AlabasterParam(path = "data/processed/02_preprocessed_lcms1"))

# 2. 恢复 SummarizedExperiment 定量矩阵并从 DelayedMatrix 转为内存 matrix
res <- readObject("data/processed/02_preprocessed_res")
for (aname in assayNames(res)) {
  assay(res, aname) <- as.matrix(assay(res, aname))
}
```
