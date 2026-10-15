# 水稻 CRISPR 敲除实验参考（含引物与靶点设计）

> **来源声明**：本文整理自公开方法文档（刘耀光课题组公开协议；Ma et al., Mol Plant 2015；
> Ma & Liu, Curr Protoc Mol Biol 2016），条目按本项目的结构重新撰写，供实验设计参考；
> 执行前请与课题组最新版官方文档核对。本包**不包含**任何未授权的第三方文件。

---

## 1. 靶点设计与选择

### 1.1 设计平台（在线工具）

| 工具 | 网址 | 用途 |
|---|---|---|
| **CRISPR-GE 靶点设计** | http://skl.scau.edu.cn/targetdesign/ | 水稻靶点设计（最常用） |
| **CRISPR-P 2.0** | http://crispr.hzau.edu.cn/CRISPR2/ | 靶点设计 + 全基因组脱靶分析 |
| Cas-OFFinder | http://www.rgenome.net/cas-offinder/ | 脱靶位点全基因组搜索 |
| OligoEvaluator / mfold | oligoevaluator.com · mfold.rna.albany.edu | 靶点二级结构预测 |

### 1.2 靶点选择七原则（要点）

1. 编辑区域**优先 CDS**；正/负链打靶效率基本相当，皆可设计；
2. GC 含量尽量不低于 40%，**50–75% 打靶效率较高**；
3. 靶点内（5'-N20NGG-3' 方向）**不得有连续 ≥4 个 T**（Pol III 转录终止信号）；
4. 靶点与 sgRNA 骨架**连续配对 ≥7 bp 会抑制结合**（注意 RNA 的 U-G 配对），须避免；
5. **双靶点**构建时两靶点位置至少相距 **100 bp**，且尽量取基因中上游；
6. 一个靶点敲除多个同源基因时，选同源基因中**碱基完全相同的区域**；
7. 优先 On-score 高、脱靶位点少、Off-score 低；做 BLAST 特异性检查时，
   避免与同源序列差异少于 5 个碱基的靶点（切点附近与 PAM 有 2 个碱基差异即可特异）。

> 本包 `scripts/design/01_sgrna_design.py` 已将上述原则（除全基因组脱靶）全部代码化为
> 本地离线筛选，并输出推荐的**双靶点组合**与构建引物；全基因组脱靶请回站 CRISPR-P / Cas-OFFinder 复核。

### 1.3 受体材料验证（不可省略）

设计引物 PCR 扩增"靶点 ±100 bp"区域（本包 `02_primer_design.py --mode validation` 一键生成），
送 Sanger 测序，确认受体植株中**靶点与 NGG 无突变**；若靶点有 SNP → 用受体的序列替换靶点，
或换靶点；若 NGG 丢失 → 必须换靶点。

## 2. pYLCRISPR/Cas9 双靶点载体构建（流程概览）

```
靶点选择(2个, ≥100bp间隔) → 合成特异引物 → PCR1 重叠延伸(2 个 sgRNA 表达盒片段)
→ PCR2 加 BsaI 接头 → BsaI 酶切-连接入 pYLCRISPR/Cas9 → 转化 DH5α → 菌液 PCR
→ 测序验证 → 转化农杆菌 → 遗传转化 → T0 植株鉴定 → T1 无 T 植株筛选
```

**关键引物体系**（完整表见 `scripts/design/01_sgrna_design.py` 输出的
`construction_primers.csv`，含 shared 通用引物区块）：

| 引物 | 组成 | 说明 |
|---|---|---|
| gRT1+ / gRT2+ | 靶点序列（N19/N20）+ `gttttagagctagaaat` | 靶点 + sgRNA 骨架起始 |
| U6aT1- | 靶点方向**互补序列** + `Cggcagccaagccagca` | 靶点 1 + U6a 端 |
| U6bT2- | 靶点方向**互补序列** + `Caacacaagcggcagc` | 靶点 2 + U6b 端 |
| U-F / gR-R | 见输出表 | PCR1 首轮扩增 |
| Pps-GGL / Pgs-GG2 / Pps-GG2 / Pgs-GGR | 见输出表 | PCR2 加 BsaI 位点 |
| PB-L / PB-R | 见输出表 | 连接产物扩增检测 |
| SP-L1 / SP-R | `GCGGTGTCATCTATGTTACTAG` / `CCGACATAGATGCAATAACTTC` | 菌液 PCR 与测序验证 |

**引物规则**：G 开头的靶点取**不含该 G 的 N19**；非 G 开头取 **N20**。

**电泳判断**：单子叶表达盒 PCR 产物 LacZ-OsU6a-gRNA ≈ **831 bp**（BsaI 切后 801），
OsU6b-gRNA ≈ **515 bp**（切后 485）；连接产物扩增条带 ≈ **1500 bp**。

**其他要点**：
- pYLCRISPR/Cas9 菌种（TOP10F'）卡那 25 μg/mL；sgRNA 载体菌种（DH10B）氨苄 50 μg/mL；
- BsaI 的识别位点（GGTCTC）与切割位点隔 1 bp，引物在切割位点附近错配碱基可形成不同黏性末端，
  使多个产物有序连接；
- 转化涂板含 IPTG/X-gal；ccdB 空载致死，但 Lacz 启动子过短常白斑为主，**白色菌落须逐个验证**。

## 3. 水稻 DNA 提取

### 3.1 CTAB 法（常规，质量好）

1. 2–3 cm 叶片 + 钢珠液氮速冻，磨样机 50 Hz 60 s；加 1 mL 65 ℃ 预热的 2×CTAB，65 ℃ 水浴 30 min（摇匀 2–3 次），12000 rpm 离心 3 min；
2. 取 0.85 mL 上清 + 0.85 mL 氯仿，充分摇匀 3–5 min，静置 5 min，12000 rpm 离心 10 min；
3. 取 0.7 mL 上清 + 2/3 体积异丙醇，-20 ℃ 静置 ≥1 h，12000 rpm 4 ℃ 离心 10 min；
4. 70% 乙醇洗涤 2 次，晾干，30–200 μL ddH₂O 溶解 + 1 μL RNase 37 ℃ 1 h；
5. 0.8% 琼脂糖电泳质检，-20 ℃ 保存。

**试剂配方**（2×CTAB 100 mL）：CTAB 2 g，PVP 1 g，1 M Tris-HCl(pH 8.0) 10 mL，
0.5 M EDTA(pH 8.0) 4 mL，5 M NaCl 28 mL，β-巯基乙醇 0.1 mL，定容。

### 3.2 TPS 法（大规模分子标记筛选，T1 无 T 植株初筛推荐）

适用于**基因编辑材料 T1 世代 T-DNA 筛选**：可插秧前完成鉴定、集中取样、省地省工。

1. 1–2 cm 叶片 + 200 μL TPS 提取液，磨样 50 Hz 180 s（剩余苗按顺序水培备用）；
2. 65 ℃ 温浴 30 min（每 10 min 摇匀），10000 rpm 离心 10 min；
3. 取 150 μL 上清入预装 150 μL 无水乙醇的 1.5 mL 管，-20 ℃ ≥30 min，10000 rpm 4 ℃ 10 min；
4. 70% 乙醇洗 2 次，去残液，开盖晾干 ≥15 min，20 μL ddH₂O 溶解（30 ℃ 助溶）；

**PCR（分子标记检测）**：9.7 μL 体系（5 μL Mix + 各 0.2 μL 引物 + 4.3 μL H₂O）+ 0.3 μL DNA；
程序 95 ℃ 5 min → 36×(95 ℃ 15 s，55 ℃ 15 s，72 ℃ 100 s)；1% 琼脂糖 120 V 18 min。

**TPS 配方**（1 L）：1 M Tris-HCl(pH 8.0) 100 mL + 0.5 M EDTA(pH 8.0) 20 mL + KCl 74.5 g，定容。

## 4. 编辑体鉴定与突变解码

| 环节 | 方法 | 工具/引物 |
|---|---|---|
| T0 代（筛**有** T 植株） | PCR 检测 T-DNA（如 SP-L1/SP-R、kanaF/R、Hyg±等） | 见构建引物表 |
| T1 代（筛**无** T 植株） | TPS 法 + 分子标记大规模初筛 | 基因型特异标记 |
| 靶点突变确认 | 靶点区域 PCR + **Sanger 测序** | BioEdit 比对；双峰杂合 → **DSDecodeM** 解码（http://skl.scau.edu.cn/dsdecode/） |
| 规模化纯合体初筛 | 二代测序方案（如 HiDecode） | 输入：双端测序 reads + 参考序列 FASTA + 靶点序列（tab 分隔）；结果需 Sanger 复核 |

## 5. 与本包工具的衔接

| 步骤 | 本包工具 |
|---|---|
| 靶点设计与双靶点配对 | `scripts/design/01_sgrna_design.py`（七原则代码化 + 构建引物自动生成） |
| 受体靶点验证引物（靶点 ±100 bp） | `scripts/design/02_primer_design.py --mode validation` |
| qPCR 表达量验证（跨内含子） | `scripts/design/02_primer_design.py --mode qpcr` |
| 通用/克隆引物 | `scripts/design/02_primer_design.py --mode generic` |

> 使用说明与参数详见 [design_usage.md](design_usage.md)。
> 在线平台（skl.scau.edu.cn / CRISPR-P）用于全基因组脱靶复核与二次确认——**离线设计与在线复核结合使用**。
