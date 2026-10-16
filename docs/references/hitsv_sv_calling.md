# HitSV 跨测序平台结构变异检测——方法笔记

> 依据 Nat Biotechnol 2026 论文（Li G … Wang Y）与其开源工具公开文档**用自己的语言重写**的方法笔记；
> 未复制论文正文、图或任何第三方源码。工具为 GPL-3.0，**运行时安装、不打包**（见 §5）。

## 1. 适用场景

- 群体中**只有部分材料有基因组组装**，其余只有重测序 reads——没法整体走 minimap2+SYRI
  组装比对路线时，用 read-based 路线补。
- 需要 SV **基因分型**（0/1、1/1）而不只是"检没检到"的场景：等位基因频率、PAV 矩阵、
  SV-GWAS。分型信息是下游关联分析的基础，只检出不分型的 callset 用不起来。
- 水稻具体场景：3K-RG 等公共数据以短读长为主；NBS-LRR 抗病基因簇、反转录转座子富集区、
  着丝粒周边的大片段 SV/CNV 与耐逆、株型、品质性状关联密切；长读重测序的种质可用
  HiFi/ONT 模式拿到碱基分辨率断点。
- 超低深度（<5×）三代数据的 SV 初筛（工具宣称支持，基因型可信度需另行评估）。

## 2. 方法要点（语言重写）

**问题本质**：SV 检出来难，不是因为信号不存在，而是 reads 比对到参考时——重复背景让
断点邻域的比对不可信，短读长跨不过重复区、长读长准确但成本高，单一信号（split/clip/
discordant/depth）各有盲区，靠人工阈值拼规则的工具在重复区要么漏、要么假阳性。

HitSV 的两个核心设计（按论文与工具文档重写）：

1. **重复感知 + 信号密度感知的断点识别**：不孤立地看单个比对信号，而是把断点邻域的
   重复程度与支持信号的密度联合评估——重复区里稀疏的软剪切更可能是比对噪声；信号富集
   且邻域重复度可控的位置才进入候选。本质是用"信号 × 上下文"替代固定阈值。
2. **单倍型解析的局部组装**：在候选 SV 区域把 reads 局部重新组装成 contig，再与参考
   比对，从而碱基级重建插入序列、精修断点两侧坐标——这是"base-resolution"的来源，
   也让插入型 SV 的序列内容可解读（而非只报长度）。

**基因分型与分相**：基于局部组装的等位证据分型；未定相杂合位点默认随机分相
（`-f 1`，可用 `-f 0` 输出未分相结果）。

**平台与物种适配**（水稻迁移的关键参数）：

- LRS 误差模式预设 `-p`：`ASM`（组装结果）/ `HIFI`、`ONT_Q20`（错误率 <1%）/
  `ERR_PRONE`（1–6%）。按实际数据质量选，不要默认。
- contig→参考比对灵敏度预设 `-P`：`asm5`（<1% 分歧）/ `asm10`（约 1–5%，**水稻种内
  默认这个**）/ `asm20`（>5%）。论文明确：非人物种用更宽松的 asm10/asm20
  （其植物基准为拟南芥、大豆）。

**公开基准效果**（摘要口径，随覆盖度变化）：HG002 GIAB 上 genotype F1 相对 SOTA 提升
长读长 12–68%、短读长 3–36%、混合 13%；1KGP Phase 4 队列多调出 31.5% 的 SV；大型长读
队列中鉴定出串联重复–可移动元件复合阵列（TR-MEI）这类多等位复杂 SV——提示**重复结构
本身是 SV 的重要来源**，水稻这类富含 TE 的基因组同理。

## 3. SV 检测路线决策（数据形式 → 路线）

| 数据形式 | 推荐路线 | 工具/脚本 |
|---|---|---|
| 多数材料有基因组组装（FASTA） | 组装比对 | 本包 `02_sv_detect.sh --mode assembly`（minimap2+SYRI） |
| 已建泛基因组图（minigraph-cactus） | 图分型 | 直接用 cactus 产出 VCF 中的 SV 基因型（`01_pangenome_graph.sh`） |
| 只有短读长（Illumina/BGI-T7）大队列 | read-based SRS | HitSV SRS 模式（备选 Manta / Delly） |
| ONT/HiFi 单样本或小群体 | read-based LRS | HitSV（备选 Sniffles2 / pbsv / cuteSV） |
| 低深度三代 + 高通量二代 | read-based 混合 | HitSV hybrid 模式 |
| 超低深度三代（<5×） | HitSV | 基因型可信度需专门评估 |
| 育种群体找因果 SV | call → GT 矩阵 → 关联 | 本包 `02_sv_detect.sh` → `03_pav_analysis.R` / `04_pan_gwas.R` |

SV 与 PAV 的关系：SYRI 组装比对直接给组装间的结构差异（基因存在/缺失 = PAV）；
read-based 从 reads 信号找缺失/插入，不需组装。两者互补，理想做法是取并集后交叉验证，
或用长读数据验证短读 callset 的断点。

## 4. 与 Rice-skills 的衔接路径

| 上游输出 | 下游模块 |
|---|---|
| `02_sv_detect.sh --mode read-based` 产出的 VCF（含 GT） | `bcftools query -f '[%GT\t]\n'` 导 0/1 矩阵 → `03_pav_analysis.R` / `04_pan_gwas.R` |
| SV 打断的基因列表 | 模块 5 注释管线主要面向小变异；SV 建议用 bedtools intersect 基因模型 + 手动核查 |
| SV 位点验证 | 设计 PCR/长读长特异验证引物（模块 6 设计模块） |
| SV 与 SNP 联合 | 与模块 7 SNP-GWAS 命中比较，SV 常捕获 SNP 遗漏的因果变异 |

## 5. 落地选择

- **HitSV（主推 read-based 路线）**：GitHub Releases 提供预编译静态二进制（免编译），
  或源码 `make`。用法封装在 `02_sv_detect.sh --mode read-based`（含 `--demo` 合成数据
  自检）。**许可提示：GPL-3.0**——本包 MIT 主体不打包其源码；用户自行安装二进制不构成
  再分发，若改装其源码并再分发需同样以 GPL-3.0 开源。
- **备选工具**（短读长 Manta / Delly；长读长 Sniffles2 / pbsv / cuteSV）：短读大队列、
  或 HitSV 安装受阻时使用；各自的许可需在安装前自行核验。
- **组装路线（零新增依赖）**：材料有组装就用本包默认 minimap2+SYRI，成本最低。
- **报告规范**：caller 名+版本、参考基因组版本与 preset（asm5/10/20）、覆盖度、
  基因型方法、验证方式（PCR/长读长/IGV 人工检查）——审稿常问，缺一不可。

## 6. 参考（文献，非源码）

- Li G, Liu Y, Zhang A, Lin X, Hu Z, Liu B, Qian L, Wang Y. Structural variant detection
  across sequencing technologies. *Nat Biotechnol* (2026).
  <https://doi.org/10.1038/s41587-026-03332-9>
- 工具与文档：<https://github.com/hitbc/HitSV>（GPL-3.0；含 HG002 demo 数据与多进程脚本）
- 预印本（2026-02，CC BY 4.0）：<https://www.researchsquare.com/article/rs-8913458>
- 官方 call results（HG002 / 1KGP 等）：<https://github.com/hitbc/HitSV_call_results>
- 水稻 SV/PAV 背景：Wang W, et al. Genomic variation in 3,010 diverse accessions of Asian
  cultivated rice. *Nature* (2018) 555:398–403.
- SYRI：Goel M, et al. SyRI: finding genomic rearrangements and local sequence differences
  from whole-genome assemblies. *Bioinformatics* (2019) 35:5066–5068.
- Sniffles：Sedlazeck FJ, et al. Accurate detection of complex structural variations using
  single-molecule sequencing. *Nat Methods* (2018) 15:461–468.
- Manta：Chen X, et al. Manta: rapid detection of structural variants and indels for
  germline and cancer sequencing applications. *Bioinformatics* (2016) 32:1220–1222.
