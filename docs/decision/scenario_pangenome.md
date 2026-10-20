# 场景卡：泛基因组与结构变异

## 1. 适用与前提

- 数据：多个染色体级组装（建议 ≥5，理想 ≥10）+ 统一注释；
- 前提：**注释流程与基因 ID 映射统一**，否则 PAV 矩阵不可比；
- 目标：核心/可变/私有基因、SV/PAV、Pan-GWAS、属级比较。

## 2. 先判断"能不能做泛基因组"

| 检查 | 不通过时的替代 |
|---|---|
| 组装连续性（N50/ BUSCO） | 碎片化 → 先补组装，或改用 SV 调用 |
| 注释来源一致 | 不一致 → 统一重注释或只做 SV |
| 参考基因组明确 | 未定 → 先确定 IRGSP-1.0 或相应组装 |
| 样本代表性 | 近缘重复 → 增加多样性材料 |

## 3. 三档方案

### 最小可行
1. `pangenome/01_pangenome_graph.sh`（minigraph/cactus/PGGB）；
2. 统一注释与基因 ID 映射；
3. 组装质量评估（BUSCO/连续性）。

### 标准
1. `pangenome/02_sv_detect.sh`（minimap2 + SYRI）；
2. `pangenome/03_pav_analysis.R`（核心/可变/私有）；
3. `pangenome/04_pan_gwas.R`（Pan-GWAS）。

### 进阶
1. `pangenome/05_oryza_pangenome_pipeline.sh`（属级 Snakemake 阶段运行器）；
2. `pangenome/06_annosine.sh`（SINE 注释）；
3. 与表达/表型联合。

## 4. 文献依据（白名单）

- 野生-栽培稻泛基因组（Nature 2025，145 个染色体级组装）——泛基因集与 PAV 的设计范式；
- 图泛基因组 GWAS（Nature Communications 2025）——把 graph 作为参考做关联；
- TE 插入多态性泛基因组（Nature Communications 2025）——TE-PAV 作为标记；
- PanOryza 泛基因目录（Genome Research 2026）——基因目录与 ID 映射；
- 属级泛基因组流程（Genome Biology 2023，MIT）——阶段化流程的实现参考。

## 5. 失败模式

| 现象 | 原因 | 处理 |
|---|---|---|
| PAV 矩阵大量缺失 | 注释不一致 | 统一注释/ID 映射后重算 |
| SV 数量异常多 | 组装错误/嵌合 | 过滤低支持 SV，多种调用交叉验证 |
| 核心基因数过少 | 组装碎片化 | 提高组装质量或降低阈值（需说明） |
| Pan-GWAS 无信号 | 群体小/PAV 稀疏 | 扩大群体，改用 SV 基因型 |

## 6. 脚本映射

`scripts/pangenome/01–06` + `scripts/gwas/03,09`（关联与多位点）。
