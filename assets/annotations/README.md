# 内置水稻注释数据快照（assets/annotations/）

本目录是 Rice-skills 内置的**官方水稻注释资源快照**，开箱即用，无需联网下载。
所有文件均来自 RAP-DB / PlantTFDB 的官方发布数据，仅做格式整理，未修改任何内容。

## 文件清单与用途

| 文件 | 内容 | 直接用途 |
|---|---|---|
| `msu_rapdb_map.csv` | 34,133 对 MSU(RAP7) ↔ RAP-DB(IRGSP) 基因 ID 映射（两列：msu,rapdb） | ID 对照查阅 |
| `rapdb_to_msu_map.csv` | 同上，列名 from,to（RAP-DB → MSU） | `adapt_company_data.py --id-map` |
| `msu_to_rapdb_map.csv` | 同上，列名 from,to（MSU → RAP-DB） | 反向转换 |
| `rice_TF_PlantTFDB_msu.csv` | 1,862 个水稻转录因子 / 56 个家族（MSU ID） | `06_tf_analysis.R --tf-annotation` |
| `rice_TF_PlantTFDB.csv` | 同一份 TF 列表的原始列（tf_id,gene_id,family） | 查阅/自定义处理 |
| `rice_TF_families.gmt` | 53 个 TF 家族基因集（GSEA GMT 格式，基因集 ≥2 个基因） | `04_gsea.R --gmt` |
| `raw_RAP-MSU_2018-03-29.txt.gz` | RAP-DB 官方原始转换表（gz 压缩，供复核） | 溯源留存 |
| `raw_Osj_TF_list.txt.gz` | PlantTFDB 官方原始 TF 列表（gz 压缩，供复核） | 溯源留存 |

## 数据来源与溯源

- **基因 ID 转换表**：RAP-DB（Rice Annotation Project Database）官方归档发布
  `RAP-MSU_2018-03-29`（表头 `rapdb / msu7`）。
  主站 `rapdb.dna.affrc.go.jp` 对部分地区/网络的自动化请求返回 HTTP 403
  （浏览器访问正常），本快照经由镜像副本获取，内容与该官方版本一致（34,133 对有效映射，
  官方表中 `None` 行已剔除）。
- **转录因子列表**：PlantTFDB（`planttfdb.gao-lab.org`）官方下载
  `Osj_TF_list.txt.gz`（Oryza sativa Japonica），2026-10-14 直接获取。

## 使用建议与更新

```bash
# 直接使用（无需任何下载）：
Rscript scripts/transcriptomics/06_tf_analysis.R --de-all de_all.csv \
    --tf-annotation assets/annotations/rice_TF_PlantTFDB_msu.csv
Rscript scripts/transcriptomics/04_gsea.R --de-all de_all.csv \
    --gmt assets/annotations/rice_TF_families.gmt
python scripts/utils/adapt_company_data.py --mode transcriptome --input vendor.csv \
    --id-map assets/annotations/rapdb_to_msu_map.csv

# 需要更新到官方最新版本时（默认优先使用本快照，--force 强制联网更新）：
python scripts/utils/fetch_rice_annotations.py --force
#   若网络拦截 TLS：加 --insecure
#   若主站 403：用浏览器下载官方文件后，通过 --rapdb-url / --tf-url 传入本地路径
```

## 许可与引用

RAP-DB 与 PlantTFDB 均为公开发布的科研数据库，其 ID 映射与注释列表为事实性数据，
供学术研究使用。请在论文中按各数据库官方指引引用：

- RAP-DB: <https://rapdb.dna.affrc.go.jp>
- PlantTFDB: <http://planttfdb.gao-lab.org>
