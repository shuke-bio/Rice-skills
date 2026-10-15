# 启动子 motif 与 TF 匹配模块使用说明（`scripts/tf/`）

> 模块定位：从基因列表提取启动子 → k-mer 富集找候选 motif → 匹配已知 TF PWM
> → 得到"可能调控这批基因的转录因子"名单。与模块 1 的 TF 分析（表达差异角度）、
> `docs/research_workflow_tf.md`（湿实验路线图）互补。

## 输入准备

1. **基因组 FASTA + GFF3**（用于按基因位置提取启动子序列；水稻推荐 RAP-DB/IRGSP-1.0 或 MSU7 注释）
2. **目标基因列表**（如 DEG、WGCNA 模块基因、GWAS 候选基因的 `gene_id` 列）
3. （第 2 步用）**PWM 库**：水稻 TF motif 集合，MEME 格式——从 PlantPAN 4.0 或
   Cis-BP 导出（见 `docs/rice_databases.md`/`docs/tools_navigation.md`）

## 第 1 步：启动子 k-mer 富集 `01_promoter_kmer.py`

```bash
python scripts/tf/01_promoter_kmer.py \
    --genome genome.fa --gff annotation.gff3 \
    --genes candidate_genes.csv \
    --up 2000 --k 6 \
    --outdir results/tf/01_promoter_kmer
```

- 提取每个基因 TSS 上游 `--up`（默认 2 kb）启动子；链感知、双链 k-mer 合并为规范方向
- 背景 = **全部注释基因自动扣除候选**（`--max-bg` 控制抽样上限，默认 3000），无需另备背景文件
- Fisher 精确检验比较候选 vs 背景 k-mer 频率，BH 校正；输出 `kmer_enrichment.csv`
  （kmer/n_cand/freq/fold/p/fdr）+ top k-mer 条形图 + `promoters.fasta`
- `--demo` 自检：注入的 motif 排名 #1（FDR 3e-9）

## 第 2 步：motif → 已知 TF 匹配 `02_motif_tf_match.py`

```bash
python scripts/tf/02_motif_tf_match.py \
    --kmers results/tf/01_promoter_kmer/kmer_enrichment.csv \
    --pwm data/motifs/oryza_plantpan.meme \
    --norm-cut 0.6 --outdir results/tf/02_motif_match
```

- 解析 MEME PWM 库（PlantPAN/Cis-BP 导出格式）；对每条 top k-mer 做双链滑窗
  log-odds 匹配（0 阶均匀背景，可 `--bg` 改），输出归一化分数
- 输出 `motif_matches.csv`（每条 k-mer 最佳命中）+ `tf_ranking.csv`（TF 按命中
  k-mer 数及最佳分排序——**候选调控 TF 名单**）
- 注意：这是快速筛查匹配（非 TOMTOM E-value 检验）；发表级结论请用
  TOMTOM/PlantPAN 复验，或结合模块 1 的差异表达 TF 名单交叉

## 结果解读与衔接

| 输出 | 下游用法 |
|---|---|
| `tf_ranking.csv` 候选 TF | 与模块 1 DEG 交集 → "表达变化且 motif 富集"的双重证据 TF |
| `motif_matches.csv` | 单基因层面的启动子变异解释（配合 `qtl/04` 的启动子区 SNP） |
| top k-mer 列表 | 与 `qtl/04_qtl_region_snps.py` 的启动子区差异 SNP 对照，做"变异破坏 motif"论证 |
| 候选 TF 名单 | → `docs/research_workflow_tf.md`（湿实验验证路线：表达/亚细胞/互作） |

## 常见问题

- **PWM 库从哪来**：PlantPAN 4.0（水稻启动子/TF 面板）或 Cis-BP 下载后为 MEME 格式
  直接可用；文件不进本仓库（版权），请自行从官方获取
- **k 值怎么选**：5–7 起步；短 k 假阳性高、长 k 覆盖低，可同时跑 5/6/7 比较
- **背景基因集**：建议用"检测到的全部基因"而非全基因组注释，减少表达偏向
