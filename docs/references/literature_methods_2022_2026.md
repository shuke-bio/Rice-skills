# 水稻生信顶刊方法索引（2022–2026）：S/A1/A2/A3/B 全量蒸馏

> 配套文件：`落地性分级_SAB.csv`（逐条仓库/许可/理由）、
> `水稻生信顶刊方法调研_2022-2026.md`（检索方法与统计）。
> 本文件是**可直接使用的蒸馏笔记**：每个方法给出「做什么 / 输入输出 / 落地方式 / 引用」，
> 论文原文未复制，方法细节以自有语言重述。

**分级口径**（按可执行程度）：S=代码+许可+契约清晰、可出三件套；
A1=可封装但重依赖；A2=无仓库、蒸馏笔记；A3=湿实验协议/设计模板；B=资源登记。

---

## 一、S 级（5 项，脚本已落地）

### S1. 水稻泛基因组分析流程（Genome Biology 2023）
- **做什么**：从多个组装出发，做共线性/泛基因集/PAV/SV 的群体级泛基因组分析。
- **输入→输出**：多份基因组 FASTA + GFF → 泛基因矩阵、PAV 表、SV 集合。
- **落地**：`scripts/pangenome/05_oryza_pangenome_pipeline.sh`（Snakemake 阶段运行器，
  MIT，上游 `LengFeng00/oryza-pangenome-pipeline`）；依赖 snakemake 与各阶段工具。
- **注意**：上游是研究流程集合，不是打包工具；按 `--stage` 逐段运行。

### S2. 3VmrMLM 族（已集成，v3.3.1）
- 多位点 QTN / QEI / 上位性 / 二分类 / 多亲本，见
  `docs/references/3vmrmlm_multilocus_qtn.md` 等四篇笔记与 `scripts/gwas/09-11`。

### S3. MeNet 多环境基因组预测（Plant Communications 2025）
- **做什么**：混合效应深度网络；双嵌入（遗传相关随机效应 + 变异累积固定效应），
  自适应调整两者贡献，支持跨环境迁移与可解释性归因。
- **输入→输出**：基因型矩阵 + 表型（多环境）→ 预测值、模块贡献度。
- **落地**：`scripts/phenomics/06_menet_gs.py`（学术专用许可，运行时安装，不打包）；
  水稻 12 性状×3 环境验证，34/36 评测最优。
- **引用**：Li Y, Ren S, Li J, Lee J, Wan J, Gan X. Plant Commun 2025.
  doi:10.1016/j.xplc.2025.101620

### S4. ELMERF 耐盐表型（Plant Phenomics 2026）
- **做什么**：MMSegmentation 语义分割，分离绿苗/黄苗/根，提取器官级表型，
  用于耐盐评价与遗传定位。
- **输入→输出**：水培 RGB 图像 → 分割掩膜、绿/黄/根面积分数等性状表。
- **落地**：`scripts/phenomics/07_elmerf_salt.py`（Apache-2.0，不打包，需作者权重）。

### S5. Chloro-Count（New Phytologist 2025）
- **做什么**：分割/计数/测量叶绿体与维管束鞘细胞，量化叶绿体占有率。
- **输入→输出**：共聚焦 .lif/图像栈 → 分割图、计数与测量表。
- **落地**：`scripts/phenomics/08_chlorocount.py`（BSD-3-Clause，模型需从上游下载）。

### S6. AnnoSINE（Plant Physiology 2022）
- **做什么**：植物基因组 SINE 注释；同源/结构/混合三模式，产出非冗余 SINE 库与全基因组注释。
- **输入→输出**：基因组 FASTA → `Seed_SINE.fa` + `.out` 注释。
- **落地**：`scripts/pangenome/06_annosine.sh`（MIT，需 blastn/TRF/node）。
- **注意**：论文以水稻 SINE 库为基准；水稻基因组实测前先在小染色体上试跑。

---

## 二、A1 级（4 项：可封装但重依赖）

### A1-1. PDLLMs（Molecular Plant 2025）
- 植物 DNA 大语言模型族（含水稻）；启动子/调控元件/变异效应预测。
- 代码+权重：`zhangtaolab/Plant_DNA_LLMs`，**CC BY-NC-SA 4.0**（禁商用、相同方式共享）。
- 落地：推理封装（GPU），权重运行时下载；文档需写明非商业限制。

### A1-2. PlantCaduceus（PNAS 2025）
- 16 个被子植物基因组（**含水水稻**）预训练的 DNA LM；跨物种剪接位点/翻译起始/
  变异有害性预测；零样本整合 GWAS 定位因果变异。
- 仓库：`plantcad/plantcad-dev`（**无 LICENSE 文件**，先核再谈封装）。
- 落地：仅笔记 + 待许可明确后再推理封装。

### A1-3. RicePilaf（GigaScience 2024）
- Post-GWAS/QTL 仪表板：整合泛基因组、共表达、调控、表观、本体、通路与文本挖掘。
- 仓库：`bioinfodlsu/rice-pilaf`，**MIT**；Shiny/Docker 部署较重。
- 落地：本地启动指引 + 输出解读笔记。

### A1-4. SeedLLM·Rice（Molecular Plant 2025）
- 水稻知识图谱 + 大语言模型；面向育种知识问答与基因-性状关联检索。
- 未找到公开代码；落地为方法笔记（架构、语料、评测方式）。

---

## 三、A2 级（39 项：蒸馏笔记）

> 每条格式：**方法（期刊 年）** — 要点 / 输入输出 / 建议。

### 泛基因组与群体遗传
1. **图泛基因组 GWAS**（Nat Commun 2025）— 把 graph pangenome 作为参考做 meta-GWAS，
   提升稀有等位与 SV 的关联功效；输入 VCF+graph，输出位点-性状关联。建议：先构建
   graph（pggb/odgi），再接关联。
2. **TE 插入多态性泛基因组**（Nat Commun 2025）— 以 TE 插入存在/缺失为标记做关联，
   解析耐冷；输入 TE 注释+群体，输出 TE-PAV 关联。
3. **泛基因组倒位索引**（Nat Commun 2023）— 系统鉴定倒位并建索引，分析亚群结构；
   输入组装，输出倒位区间与基因型。
4. **syntelog 泛基因组**（Genome Biology 2023）— 以 syntelog 定义同源基因，
   重建驯化/去驯化历史；输入多组装+注释。
5. **Meta-QTL**（Rice 2026）— 整合多研究 QTL 区间，缩小区间并配 InDel 标记；
   输入文献 QTL 表，输出一致性区间。
6. **TWAS**（Nat Commun 2023）— 用表达预测模型把 GWAS 信号映射到基因表达，
   定位穗型因果基因；输入 GWAS+表达参考面板。
7. **regulome 建图与建模**（Nat Commun 2024）— 整合 ATAC/组蛋白/DNase 构建调控图谱，
   建模复杂性状；输入多表观数据。
8. **单倍型关联**（Rice 2023）— 以基因为单位做单倍型-表型关联，比单 SNP 更稳健；
   输入基因型+表型。

### 单细胞、空间与表观
9. **空间转录组流程**（Mol Plant 2025）— 根系干旱适应的空间表达分析；
   输入空间矩阵+图像，输出区域特异基因。
10. **scATAC 流程**（Nat Plants 2025）— 10.4 万核、126 细胞态的染色质可及性图谱；
    输入 fragments，输出 peak-by-cell 与细胞态。
11. **COPILOT 管线 + 根图谱整合**（Nature 2025）— kallisto|bustools 预处理 + 参考图谱
    标签迁移；输入 FASTQ，输出细胞类型比例。
12. **snRNA 流程**（JIPB 2026）— 单核转录组解析种子早期发育；输入核 RNA counts。
13. **3D 表观调控抽穗期**（Genome Biology 2026）— Hi-C + 表观 + 表达整合，
    解析 cis/trans 调控；输入 Hi-C 矩阵。
14. **H3K9ac 标记的 3D 互作**（Genome Biology 2024）— 干旱下超级增强子样启动子
    介导的染色质互作；输入 ChIP-seq + Hi-C。
15. **3D+表观+转录整合（生物钟）**（NAR 2023）— 三维修饰与节律基因表达的联合分析。

### 表型组与基因组预测
16. **空间效应 GP**（Rice 2025）— 在全国历史数据中加入空间效应，提升抽穗期预测；
    输入历史表型+基因型+地点。
17. **多代多点 GS**（Rice 2023）— 优化循环选择的多代模型；输入系谱+基因型。
18. **稀疏表型 GS**（Rice 2023）— 少量表型 + 单倍型模型；输入稀疏表型。
19. **SV+SNP 深度学习预测**（Plant Methods 2024）— CNN/MLP 对比 BayesC/RKHS，
    二元性状优势明显。
20. **高光谱蛋白 GWAS**（Plant Phenomics 2024）— 高光谱预测籽粒蛋白并做 GWAS，
    低成本替代化学测定。
21. **BloomSight**（Plant Phenomics 2026）— 高频开花动态表型，解析开花时间遗传。
22. **极端表型 QTL 策略**（Plant Commun 2022）— 用极端个体组合提高 QTL 检出效率。
23. **高光谱→代谢物 + mGWAS**（Genome Biology 2025）— 1848 高光谱指数预测 887 代谢物，
    再做 mGWAS 与共定位；输入高光谱+代谢组+基因型。
24. **挥发组 GWAS**（Mol Plant 2024）— 挥发物作分子表型定位香气基因。
25. **广靶挥发组学**（Mol Plant 2022）— 挥发物广靶检测方法学（含水稻香气）。
26. **遗留数据 GWAS**（Plant Physiol 2023）— 整合历史 GWAS 结果提升 QTL 可信度。
27. **GPS 数据融合**（Plant Commun 2025）— 基因组+表型选择的数据/特征/结果三种融合，
    四作物验证。
28. **AI 辅助野生稻抗病种质挖掘**（Rice Science 2025）— 图像/表型+基因组筛选抗病种质。
29. **GWAS+单倍型耐碱**（Rice Science 2025）— 联合关联与单倍型缩小候选区间。

### AI 与预测工具
30. **RiSpy**（BIB 2026）— 基于特征选择的基因编辑株系指纹识别，用于编辑事件质控。
31. **RiceSNP-ABST**（BIB 2025）— 深度学习预测非生物胁迫相关 SNP（水稻专用）。
32. **RiceSNP-BST**（BIB 2024）— 同上，生物胁迫版。
33. **MINI-EX**（Mol Plant 2022）— 植物单细胞基因调控网络推断（转录因子-靶基因）。
34. **FLASH pipeline**（Mol Plant 2022）— 阵列 CRISPR 文库构建与水稻 RLK 基因功能发现。
35. **GreenPhos**（Mol Plant 2024）— 植物磷酸化蛋白质组深度定量。
36. **TE-based markers**（Plant Commun 2022）— 转座子标记鉴定农艺性状等位变异。
37. **TE 标记泛基因组**（Plant Commun 2022）— 泛基因组尺度 TE 标记与等位变异。
38. **SeedLLM·Rice**（Mol Plant 2025）— 见 A1-4（此处登记为方法条目）。
39. **AI 辅助表型-基因型联合筛选**（Rice Science 2025）— 见第 28 条。

---

## 四、A3 级（5 项：协议/设计模板）

### A3-1. PrimeRoot 大片段定点插入（Nature Biotechnology 2023）
- **原理**：先导编辑引入位点特异性重组酶识别位点，再由重组酶介导大片段整合。
- **设计字段**：目标位点、PAM、pegRNA（含 RTT/PBS）、重组酶供体（attB/attP）、
  筛选标记；水稻已验证 11.1 kb 插入抗病基因。
- **模板**：`pegRNA: spacer(20nt)+sgRNA scaffold+RTT(10-20nt)+PBS(8-12nt)`；
  供体两端加对应重组位点。

### A3-2. uORF 梯度下调（Nature Biotechnology 2023）
- **原理**：在 5'UTR 新创 uORF 或延长内源 uORF，梯度抑制主 ORF 翻译（2.5–84.9%）。
- **设计字段**：uORF 起始密码子强度（Kozak）、与主 ORF 的距离、uORF 长度、
  终止密码子位置；水稻 OsDLT/OsTB1/OsTCP19 已验证。
- **模板**：按目标抑制强度选择 uORF 类型（新创/延长）与数量。

### A3-3. uORF 调控 DEP1/GIF1（Rice 2024）
- **原理**：编辑 DEP1/GIF1 的 uORF 提升产量相关性状；设计要点同上。

### A3-4. 编辑新体系（2024–2026）
- **PE6c**（JIPB 2024）：增强先导编辑效率；**Cas12i3**（Plant Physiol/Genome Biol）：
  小体积核酸酶；**IscB/TnpBmax**（Plant Commun/JIPB）：超紧凑编辑工具。
- **设计字段**：PAM 兼容性、编辑窗口、蛋白体积（AAV 适配）、水稻密码子优化。

### A3-5. InPACT（Plant Physiology 2026）
- **原理**：in-planta 农杆菌介导谷类转化，绕过组织培养。
- **字段**：注射时期、菌株/载体、表面活性剂、共培养条件；水稻/小麦/玉米验证。

---

## 五、B 级资源登记（16 项）

| 资源 | 期刊·年 | 类型 | 用途 |
|---|---|---|---|
| 野生-栽培稻泛基因组 | Nature 2025 | 数据 | 145 组装，3.87 Gb 新序列 |
| PanOryza 泛基因目录 | Genome Res 2026 | 数据 | 泛基因 catalog |
| 泛基因组基因分型芯片 | Plant J 2023 | 设计 | 育种芯片标记选择 |
| 超级泛基因组着丝粒图谱 | JIPB 2024 | 数据 | 着丝粒结构 |
| 非洲野生稻 T2T | GigaScience 2025 | 数据 | O. longistaminata 参考 |
| 水稻单细胞多组学图谱 | Nature 2025 | 数据 | 细胞类型参考 |
| 时空生命周期图谱 | Cell 2026 | 数据 | 空间发育参考 |
| 早期种子单核图谱 | JIPB 2026 | 数据 | 种子发育 |
| RiceSEG | Plant Phenomics 2025 | 数据 | 分割基准 |
| Rice Gene Index | Mol Plant 2023 | 数据库 | 泛基因组基因索引 |
| RAP-DB | NAR 2025 | 数据库 | 注释与组学资源 |
| RicePilaf | GigaScience 2024 | 平台 | post-GWAS 仪表板（MIT） |
| ricePSP | Genome Biology 2025 | 数据库 | 相分离蛋白 |
| ROADv2 | Rice 2024 | 数据库 | 表达谱图谱 |
| Oryza CLIMtools | Plant Commun 2024 | 平台 | 基因组-环境关联 |
| Rice3kGS / RiceTFtarget | PC 2025 / PP 2023 | 平台 | GS / TF 靶基因 |

---

## 六、许可与使用边界

- **可打包/复制**（保留 LICENSE+NOTICE）：oryza-pangenome-pipeline（MIT）、
  ELMERF（Apache-2.0）、Chloro-Count（BSD-3）、AnnoSINE（MIT）、RicePilaf（MIT）。
- **仅运行时使用，不打包**：MeNet（学术专用）、PDLLMs（CC BY-NC-SA，禁商用）、
  3VmrMLM 族（GPL-3）。
- **无许可，仅笔记/网页指引**：RiceSNP-ABST/BST、PlantCaduceus（待核）。
- 所有 A2/A3 条目只做方法蒸馏，未复制论文原文或图件。

---

## 七、落地状态与实测说明（2026-10-19）

本节记录 5 个 S 级 wrapper 的**真实**执行状态，避免把"代码已就位"误读为"结果已跑通"。

| 脚本 | 上游 | demo 实测 | 缺失依赖 |
|---|---|---|---|
| `pangenome/05_oryza_pangenome_pipeline.sh` | Oryza 泛基因组（MIT） | **SKIP** | snakemake、minimap2、samtools |
| `pangenome/06_annosine.sh` | AnnoSINE（MIT） | **SKIP** | blastn、TRF（node 已具备） |
| `phenomics/06_menet_gs.py` | MeNet（学术专用） | **SKIP** | PyTorch |
| `phenomics/07_elmerf_salt.py` | ELMERF（Apache-2.0） | **SKIP** | PyTorch、mmsegmentation、作者权重 |
| `phenomics/08_chlorocount.py` | Chloro-Count（BSD-3） | **SKIP** | PyTorch、作者权重 |

说明：
1. 5 个 demo 都已实际执行，均生成合成 fixture 与依赖报告，并在缺依赖时打印安装命令后
   以 SKIP 退出（退出码 0），**没有伪装成功**；
2. `05` 的上游是**研究流程集合**（Snakefile + 分阶段脚本），不是打包工具，故封装为
   `--stage` 阶段运行器；真实运行需按上游 `envs/` 建环境；
3. 真实跑通的最小环境：conda + PyTorch(CPU) + mmsegmentation/mmcv（ELMERF）、
   BLAST+/TRF（AnnoSINE）、snakemake + 比对/组装工具（泛基因组流程）；
4. 模型权重（ELMERF、Chloro-Count）由原作者分发，本仓库不转存。
