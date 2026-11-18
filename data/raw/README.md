# `Metabonaut` 原始数据缓存目录说明 (`data/raw/`)

本目录按照 `REPRO-D` 生信工程规范保留为**只读原始数据挂载点**。目前目录内仅存放 `.gitkeep` 与本说明文件，未手动落盘大体积 `.mzML` 原始质谱文件，原因与真实存储架构如下：

---

## 1. 为何本目录不直接存放 `.mzML` 文件？

1. **自动缓存托管机制**：[01_metabonaut_download_and_qc.R](../../repro_scripts/01_metabonaut_download_and_qc.R) 与 [03_metabonaut_norm_limma_annotation.R](../../repro_scripts/03_metabonaut_norm_limma_annotation.R) 使用 Bioconductor 官方 **`MsBackendMetaboLights`** + **`MsIO`**（`MetaboLightsParam(mtblsId = "MTBLS8735", ...)`）按需从 EMBL-EBI MetaboLights API 自动拉取数据，并由 **`BiocFileCache`** 自动持久化在系统用户本地缓存中。
2. **避免 Git 仓库膨胀**：LC-MS1 与 LC-MS/MS 原始 `.mzML` 体积较大且为只读上游公共资源，遵循 `.gitignore`（屏蔽 `**/data/raw/**` 下的质谱二进制与大数据文件，仅追踪本说明文档与 `.gitkeep`）。

---

## 2. 实际托管的 `MTBLS8735` 原始数据集清单

| 采集模式 | Assay 文件名 | 样本设计 | 谱图与用途 |
|---|---|---|---|
| **HILIC 正离子 LC-MS1** | `a_MTBLS8735_LC-MS_positive_hilic_metabolite_profiling.txt` | **10 个人血浆/血清样本**：3 CVD (`A, D, F`) + 3 CTR (`B, C, E`) + 4 Pooled QC (`POOL1–POOL4`) | 共 `17,210` 张一级谱图，用于 Step 01 BPC/TIC/内标 EIC 与 Step 02 `xcms` 峰检测定量 |
| **HILIC 正离子 LC-MS/MS** | `a_MTBLS8735_LC-MSMS_positive_hilic_metabolite_profiling.txt` | 对照组 (`CTR`) 含入母离子列表（Inclusion List）的二级质谱（`CE20`, `CE30`, `CES` 碰撞能量） | 用于 Step 03 对 `limma` 显著特征执行 `matchSpectra()` 余弦相似度 MS2 鉴定（锁定 `FT0845` = Caffeine） |

---

## 3. 查看本机缓存路径

如需在本地确认已下载缓存的 `.mzML` 与 `AnnotationHub` 数据库所在位置，可在 R 中执行：

```r
# 查看 MsBackendMetaboLights 缓存目录
tools::R_user_dir("MsBackendMetaboLights", which = "cache")

# 查看 BiocFileCache / AnnotationHub (MassBank AH116166) 缓存目录
AnnotationHub::getAnnotationHubOption("CACHE")
```
