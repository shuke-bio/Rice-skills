# 运行时外部依赖与许可核验清单

> 本仓库**不打包**下列工具；脚本只生成命令或调用本机安装版本。
> 表中"许可"列一律标注为**待核验**：在对外分发、商业部署或打包镜像前，必须打开官方
> LICENSE 原文确认，并把结论写入 `third_party/NOTICE.md`。不要仅依据 GitHub 标签判断。

| 工具 | 用途 | 官方入口 | 许可 |
|---|---|---|---|
| fastp | FASTQ 质控 | https://github.com/OpenGene/fastp | 待核验 |
| HISAT2 | RNA-seq 比对 | https://github.com/DaehwanKimLab/hisat2 | 待核验 |
| samtools | BAM 处理 | https://github.com/samtools/samtools | 待核验 |
| featureCounts (subread) | 基因计数 | https://subread.sourceforge.net/ | 待核验 |
| PLINK | GWAS 质控 | https://www.cog-genomics.org/plink/ | 待核验 |
| rMVP | GWAS（GLM/MLM/FarmCPU） | https://github.com/xiaolei-lab/rMVP | 待核验 |
| GAPIT | GWAS 多模型 | https://github.com/jiabowang/GAPIT | 待核验 |
| Cas-OFFinder | 全基因组脱靶 | http://www.rgenome.net/cas-offinder/ | 待核验 |
| CRISPOR | 脱靶/效率 | http://crispor.tefor.net/ | 待核验 |
| vg | 图泛基因组比对 | https://github.com/vgteam/vg | 待核验 |
| PanGenie | 图基因型 | https://github.com/eblerjana/pangenie | 待核验 |
| minigraph-cactus | 图构建 | https://github.com/ComparativeGenomicsToolkit/cactus | 待核验 |
| odgi | 图操作/可视化 | https://github.com/pangenome/odgi | 待核验 |
| cooler | Hi-C 矩阵 | https://github.com/open2c/cooler | 待核验 |
| HiCExplorer | TAD/loops | https://github.com/deeptools/HiCExplorer | 待核验 |
| Sniffles2 | 长读长 SV | https://github.com/fritzsedlazeck/Sniffles | 待核验 |
| cuteSV | 长读长 SV | https://github.com/tjiangHIT/cuteSV | 待核验 |
| Delly | 短读长 SV | https://github.com/dellytools/delly | 待核验 |
| Manta | 短读长 SV | https://github.com/Illumina/manta | 待核验 |
| SURVIVOR | SV 合并 | https://github.com/fritzsedlazeck/SURVIVOR | 待核验 |
| selscan | iHS/CLR | https://github.com/szpiech/selscan | 待核验 |
| Dsuite | 渐渗/Dsuite | https://github.com/millanek/Dsuite | 待核验 |
| SIRIUS | 代谢物注释 | https://bio.informatik.uni-jena.de/software/sirius/ | 待核验 |
| GNPS | 质谱库匹配 | https://gnps.ucsd.edu/ | 待核验 |
| susieR | fine-mapping | https://github.com/stephenslab/susieR | 待核验 |
| rrBLUP | 基因组选择 | https://github.com/jeffreyendelman/rrBLUP | 待核验 |
| BGLR | BayesB/贝叶斯 GS | https://github.com/variani/BGLR | 待核验 |
| GSVA | ssGSEA/GSVA | https://github.com/rcastelo/GSVA | 待核验 |
| harmonypy | 单细胞整合 | https://github.com/slowkow/harmonypy | 待核验 |
| celltypist | 细胞注释 | https://github.com/Teichlab/celltypist | 待核验 |
| scrublet | 双细胞检测 | https://github.com/swolock/scrublet | 待核验 |
| cell2location | 空间反卷积 | https://github.com/BayraktarLab/cell2location | 待核验 |
| COMMOT | 空间通讯 | https://github.com/zcang/COMMOT | 待核验 |
| BE-Hive | 碱基编辑效率 | https://github.com/maxwshen/be-hive | 待核验 |
| PRIDICT | 先导编辑效率 | https://github.com/uzh-dqbm-cmi/PRIDICT | 待核验 |

## 使用原则

1. **不打包**：本仓库不复制上述工具的源码或二进制；
2. **运行时安装**：由用户自行 `conda`/`pip`/`apt` 安装；
3. **GPL 隔离**：GPL 系工具只通过命令行调用，不并入 MIT 主体；
4. **核验后再分发**：如需制作 Docker 镜像对外发布，逐一核对 LICENSE 并写入 NOTICE；
5. **商业使用**：部分工具对商业用途有限制，需单独确认。
