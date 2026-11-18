# `Metabonaut` 数据架构与数据字典 (`data/`)

本目录遵循 `REPRO-D` 工程规范组织，将外部原始质谱缓存、先验内标参考表、以及 `xcms` 预处理阶段产出的中间对象分层隔离：

- [**`./raw/`**](./raw/README.md)：只读原始数据挂载点（由 `MsBackendMetaboLights` + `BiocFileCache` 自动缓存管理 `MTBLS8735`，参见 [data/raw/README.md](./raw/README.md)）。
- [**`./processed/`**](./processed/README.md)：Alabaster HDF5 + `.rds` 中间产物目录，存储 `Step 01–02` 滤波与峰表构建结果，参见 [data/processed/README.md](./processed/README.md)。
- [**`./intern_standard_list.txt`**](./intern_standard_list.txt)：随仓库直接追踪的 `20` 种稳定同位素标记内标精确质量与保留时间先验表。

---

## 1. 全项目外部与本地数据集登记表

| 数据集 / 文件 | 来源平台 / Accession | 样本设计与生物学背景 | 格式与规模 | 消费脚本 |
|---|---|---|---|---|
| **MTBLS8735 (LC-MS1)** | EMBL-EBI MetaboLights [`MTBLS8735`](https://www.ebi.ac.uk/metabolights/MTBLS8735) | 人心血管疾病 (CVD) 血清/血浆 HILIC 正离子代谢组：3 CVD (`A, D, F`) + 3 CTR (`B, C, E`) + 4 Pooled QC (`POOL1–4`) | `.mzML`（10 样本，17,210 张一级谱图） | `Step 01`–`03` |
| **MTBLS8735 (LC-MS/MS)** | EMBL-EBI MetaboLights [`MTBLS8735`](https://www.ebi.ac.uk/metabolights/MTBLS8735) | CTR 对照样本含入母离子列表（Inclusion List）MS2 采集，含 `CE20`, `CE30`, `CES` 碰撞能量 | `.mzML`（MS2 碎片谱图） | `Step 03` |
| **MassBank 参考谱库** | Bioconductor `AnnotationHub` (`AH116166`) | 公共高分辨小分子精确质量与正离子 MS2 碎片谱图库 | `CompDb` SQLite（由 `AnnotationHub` 缓存） | `Step 03` |
| **同位素内标参考表** | 源自 `upstream_Metabonaut/vignettes/intern_standard_list.txt` | 20 种稳定同位素标记氨基酸/代谢物（如 `cystine_13C_15N`、`methionine_13C_15N`），用于评估仪器 RT 漂移与 `CentWave` 峰宽设定 | TSV（20 行 × 17 列，2.6 KB） | `Step 01`–`02` |
| **冒烟模拟数据集** | 运行期由 `00_smoke_test.R` 内存构造 | 符合真实非靶向统计规律的对数正态峰面积矩阵（300 特征 × 18 样本：6 CTR + 6 Case + 6 QC），前 40 个注入差异信号 | 内存 `SummarizedExperiment` | `00_smoke_test.R` |

---

## 2. `intern_standard_list.txt` 字段字典

该表共 **20 种内标 × 17 列**，核心列说明如下：

- `name` / `abbreviation`：标准名称（如 `cystine_13C_15N` / `L-Cystine (13C6, 99%; 15N2, 99%)`）
- `formula`：中性分子式及单同位素分子量列
- `POS` / `NEG`：正/负离子模式下优先加合离子形式（如 `[M+H]+`, `[M+Na]+`）
- `RT`、`rtmin`、`rtmax`：预期色谱保留时间及 EIC 抽提窗口（s，正负 10s 容差）
- `mz`、`mzmin`、`mzmax`：加合离子质荷比及 EIC 抽提质量窗口（Da，正负 `0.005 Da`）
- `quality_POS` / `quality_NEG`：原实验室标注峰形质量评级（`0_perfect` 至 `3_bad`）
