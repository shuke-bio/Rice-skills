#!/usr/bin/env python3
###############################################################################
# 01_make_plan.py -- 组学分析方案生成器（决策层）
#
# 作用：读取"数据契约"（JSON），先做可行性校验，再输出三档方案
#       （最小可行 / 标准 / 进阶），每档映射到本仓库真实脚本，并给出
#       前提、失败模式、回退建议与方法依据（仅白名单期刊）。
#
# 用法：
#   python scripts/plan/01_make_plan.py --input my_data.json --outdir results/plan
#   python scripts/plan/01_make_plan.py --demo
#
# 支持 assay：bulk_rna | metabolome | resequencing | rna_metabolome |
#             single_cell | spatial | proteome | wgbs | atac | pangenome
# 支持 goal ：differential_expression | marker_discovery | gwas | prediction |
#             integration | qtl
#
# 输出：plan.md（人读）、plan.json（机读）、commands.sh（执行清单）
###############################################################################
import argparse, json, os, re, sys
from datetime import date

ASSAYS = {"bulk_rna", "metabolome", "resequencing", "rna_metabolome",
          "single_cell", "spatial", "proteome", "wgbs", "atac", "pangenome"}
GOALS = {"differential_expression", "marker_discovery", "gwas", "prediction",
         "integration", "qtl"}

WHITELIST_JOURNALS = {
    "Nature", "Science", "Cell", "Nature Biotechnology", "Nature Methods",
    "Nature Genetics", "Nature Plants", "Nature Communications", "Science Advances",
    "PNAS", "Molecular Plant", "Plant Cell", "Plant Communications", "Plant Physiology",
    "JIPB", "New Phytologist", "Plant Journal", "Briefings in Bioinformatics",
    "Nucleic Acids Research", "Genome Biology", "Genome Research", "Bioinformatics",
    "GigaScience", "Plant Phenomics", "Crop Journal", "Rice",
}

def normalize(d):
    missing = [k for k in ("assay", "goal", "design", "n_per_group") if k not in d]
    if missing:
        sys.exit(f"缺少必填字段: {', '.join(missing)}")
    if d["assay"] not in ASSAYS:
        sys.exit(f"assay 必须是 {sorted(ASSAYS)} 之一")
    if d["goal"] not in GOALS:
        sys.exit(f"goal 必须是 {sorted(GOALS)} 之一")
    d.setdefault("replicates", True)
    d.setdefault("batch", False)
    d.setdefault("phenotype", "none")
    d.setdefault("environments", 1)
    d.setdefault("factors", [])
    d.setdefault("reference", "")
    d.setdefault("quantification", "untargeted")
    d.setdefault("paths", {})
    return d

def load(path):
    with open(path, encoding="utf-8") as f:
        return normalize(json.load(f))

# --------------------------------------------------------------------------- #
# 校验：返回 (blockers, warnings)
# --------------------------------------------------------------------------- #
def validate(d):
    b, w = [], []
    n = int(d["n_per_group"])
    assay, goal = d["assay"], d["goal"]
    if assay in ("bulk_rna", "rna_metabolome") and goal == "differential_expression":
        if n < 2 or not d["replicates"]:
            b.append(f"每组 n={n} 且 replicates={d['replicates']}：无法估计组内方差，"
                     "DESeq2/edgeR 的 P 值不可解释。")
        elif n < 3:
            w.append(f"每组 n={n}：可做推断但功效低；结论定位为“候选”，需独立验证。")
        if d["batch"] and len(d["factors"]) == 1:
            w.append("存在批次且处理是唯一因素：若批次与分组完全混杂，效应不可分离。")
    if assay in ("single_cell", "spatial"):
        if n < 2:
            b.append(f"样本数 n={n}：单细胞/空间无法做组间比较（至少 2 个生物学样本）。")
        w.append("细胞层面 n 很大不等于生物学重复；组间比较以样本为单位，需先做批次整合。")
        w.append("线粒体比例/基因数阈值因组织而异，需按组织设定并报告。")
    if assay in ("wgbs", "atac"):
        if n < 2:
            b.append(f"样本数 n={n}：差异甲基化/差异可及性需要生物学重复。")
        w.append("WGBS 需报告覆盖度与亚硫酸盐转化率；ATAC 需报告 FRiP 与 TSS 富集。")
    if assay == "metabolome":
        if d.get("quantification") == "untargeted":
            w.append("非靶向代谢组为半定量数据：组间比较用相对变化，不做绝对浓度推断。")
        w.append("代谢物为组成型数据：相关分析前建议 CLR/对数比，避免伪相关。")
    if assay == "proteome":
        w.append("蛋白组缺失值常见（MNAR）：需报告缺失率与插补策略，缺失>50% 的蛋白不建议做差异。")
        if n < 3:
            w.append(f"每组 n={n}：蛋白组技术变异较大，建议 n≥3。")
    if assay == "pangenome":
        if n < 5:
            w.append(f"组装数 n={n}：泛基因组分析建议 ≥5 个高质量组装，否则核心/可变基因估计不稳。")
        w.append("需统一注释流程与基因 ID 映射，否则 PAV 矩阵不可比。")
    if goal == "gwas":
        if n < 100:
            w.append(f"群体规模 n={n}：GWAS 功效有限（建议 ≥200–300），可先做候选区间/单倍型。")
        if d.get("phenotype") == "binary" and n < 200:
            w.append("二分类性状对样本量更敏感，n<200 时假阳性风险高。")
        if not d["reference"]:
            w.append("未填写参考基因组/注释版本：关联结果无法与注释对齐。")
    if goal == "qtl":
        if n < 150:
            w.append(f"群体/家系 n={n}：QTL 定位建议 ≥150 个家系或个体。")
        if d["design"] not in ("ril_magic_nam", "case_control", "natural_population"):
            w.append(f"design={d['design']}：QTL 通常需要家系/多亲本设计，请确认。")
    if goal == "prediction":
        if n < 150:
            w.append(f"基因组预测 n={n}：交叉验证方差大，报告重复 CV 的均值±SD。")
    if goal == "integration" and assay != "rna_metabolome":
        w.append("整合分析至少需要两个组学层；当前 assay 只声明了一层。")
    if assay == "rna_metabolome" and n < 20:
        w.append(f"多组学整合 n={n}：MOFA/DIABLO 不稳定，建议先做相关/共表达。")
    return b, w

# --------------------------------------------------------------------------- #
# 方案模板
# --------------------------------------------------------------------------- #
def tier(name, intent, steps, outputs, risk):
    return {"tier": name, "intent": intent, "steps": steps,
            "outputs": outputs, "risk": risk}

def plan_rnaseq(d):
    t1 = tier("最小可行", "先确认数据质量与整体结构，不急于做推断",
        ["01_qc_filter.py：低表达过滤、样本相关性、PCA",
         "02_deseq2.R：仅在 n>=2 且有重复时给出差异表（否则只做描述性汇总）",
         "07_plot_transcriptome.R：PCA/相关性/火山图初版"],
        ["表达矩阵质控报告", "PCA 与样本聚类图", "差异表（若可用）"],
        "无重复时不得报告 P 值；批次未建模时差异可能来自技术因素")
    t2 = tier("标准", "可解释的差异表达 + 通路层面结论",
        ["01_qc_filter.py → 02_deseq2.R（设计公式含批次）",
         "03_enrichment.R（GO/KEGG）→ 04_gsea.R（排序基因集）",
         "06_tf_analysis.R（差异 TF 与靶基因）→ 07_plot_transcriptome.R",
         "敏感性分析：edgeR 与 DESeq2 结果取交集"],
        ["差异基因表", "富集/GSEA 图", "TF-靶基因表", "出版级图表"],
        "富集结果依赖背景基因集与注释版本；需在方法中写明")
    t3 = tier("进阶", "调控网络与多组学联合（有条件时）",
        ["05_wgcna.R（n>=15 时构建模块与性状关联）",
         "integration/01_mofa.R 或 02_diablo.R（与代谢组匹配样本）",
         "integration/05_network_integration.R（基因-代谢物网络）",
         "qtl/02_colocalization.R（有基因型时做 eQTL/mQTL 共定位）"],
        ["模块-性状表", "多组学因子/权重", "候选基因-代谢物网络"],
        "样本不匹配或 n 不足时，整合结果只能作探索性")
    return [t1, t2, t3]

def plan_gwas(d):
    t1 = tier("最小可行", "先看群体结构与数据质量，再决定模型",
        ["01_plink_qc.sh：MAF/缺失/哈代温伯格过滤",
         "02_plink_to_matrix.R：导出基因型矩阵",
         "05_plot_gwas.R：QQ/曼哈顿图模板（先验证绘图链路）"],
        ["质控后 PLINK 文件", "群体结构初判（PCA）"],
        "未做 PCA/亲缘校正的关联结果不可信")
    t2 = tier("标准", "单位点关联 + 候选区间 + 单倍型",
        ["03_association.R（GLM/MLM，含 PCA 协变量）",
         "04_candidate_loci.R（lead SNP ±LD 区间）",
         "07_haplotype.R（候选基因单倍型-表型）",
         "08_2d_association.R（成对性状/上位性初筛）"],
        ["显著位点表", "候选基因列表", "单倍型箱线图"],
        "LD 衰减与群体分层会改变区间大小；需报告多重检验校正")
    t3 = tier("进阶", "多位点、环境互作与特殊设计",
        ["09_multilocus_qtn.R（3VmrMLM：多位点/QEI/上位性）",
         "10_case_control_qtn.R（二分类，Linux/WSL）",
         "11_multiparent_qtl.R（MAGIC/NAM/ROAM）",
         "06_meta_gwas.R（多环境/多群体合并）"],
        ["QTN 表（含环境互作）", "二分类关联结果", "Meta 分析结果"],
        "多位点方法对 LD 与群体结构敏感")
    return [t1, t2, t3]

def plan_multiomics(d):
    t1 = tier("最小可行", "先确认两层数据可匹配、可比较",
        ["metabolomics/01_preprocess.py：缺失填充、归一化、QC",
         "transcriptomics/01_qc_filter.py：表达矩阵质控",
         "按样本 ID 取交集，输出匹配表"],
        ["匹配样本清单", "两层各自 PCA"],
        "样本不匹配是最常见失败原因")
    t2 = tier("标准", "有监督/无监督整合，给出跨层关联",
        ["metabolomics/02_oplsda.R（组间判别与 VIP）",
         "integration/02_diablo.R（需分组标签）",
         "integration/01_mofa.R（无监督因子，n>=20 更稳）",
         "基因-代谢物 Spearman 相关 + BH 校正"],
        ["VIP 代谢物表", "整合因子/权重", "跨层相关表"],
        "组成型数据相关分析易产生假相关；需 CLR/对数比")
    t3 = tier("进阶", "网络与遗传定位",
        ["integration/05_network_integration.R（模块-代谢物网络）",
         "qtl/02_colocalization.R 或 03_coloc_bayesian.R（mQTL/eQTL 共定位）",
         "gwas/07_haplotype.R（候选位点单倍型）"],
        ["跨层网络图", "共定位证据（PP4）", "候选基因-代谢物对"],
        "共定位需要足够样本量与 LD 信息；PP4 高不等于因果")
    return [t1, t2, t3]

def plan_singlecell(d):
    t1 = tier("最小可行", "先做样本级 QC 与整合，再谈细胞类型",
        ["single_cell/01_scanpy_pipeline.py：QC（线粒体/基因数）、归一化、HVG、PCA",
         "按样本评估细胞数与中位基因数（报告每个样本，不合并后再筛）",
         "批次整合（Harmony/scVI 由脚本参数选择）"],
        ["每样本 QC 指标表", "整合后的 UMAP", "细胞类型初注释"],
        "以样本为生物学重复；细胞数不等于 n；整合过度会抹掉真实组间差异")
    t2 = tier("标准", "细胞类型组成与差异表达",
        ["聚类 + marker 注释（参考图谱标签迁移）",
         "细胞类型比例差异（以样本为单位，非细胞为单位）",
         "伪 bulk 差异表达（避免单细胞检验的伪重复）"],
        ["细胞类型比例表", "伪 bulk 差异表", "marker 热图"],
        "比例差异受采样与消化偏差影响，需报告捕获效率")
    t3 = tier("进阶", "轨迹/调控/空间映射",
        ["轨迹推断与分支基因（需明确起点）",
         "SCENIC 类调控网络（资源消耗大）",
         "空间/参考映射（若有空间或图谱数据）"],
        ["轨迹图", "调控模块", "空间映射图"],
        "轨迹方向需生物学依据；推断结果不能直接当因果")
    return [t1, t2, t3]

def plan_spatial(d):
    t1 = tier("最小可行", "先对齐图像与表达矩阵",
        ["spatial/01_spatial_pipeline.py：数据读入、质控、归一化",
         "核对 spot 坐标与组织图像配准",
         "空间域初聚类"],
        ["质控后矩阵", "空间域图"],
        "图像与坐标错位会让所有空间结论失效")
    t2 = tier("标准", "空间域 + 区域特异表达",
        ["空间域聚类与 marker",
         "区域间差异表达（以样本/切片为单位）",
         "与单细胞参考做标签迁移"],
        ["空间域 marker 表", "区域差异表", "映射图"],
        "切片间批次效应常见；不可把 spot 当独立重复")
    t3 = tier("进阶", "细胞邻域与多切片整合",
        ["邻域富集/细胞通讯（需细胞类型）",
         "多切片整合与差异区域比较",
         "与 bulk/单细胞联合验证"],
        ["邻域网络", "跨切片比较", "联合验证结果"],
        "邻域分析对分割精度敏感；结论需组织学验证")
    return [t1, t2, t3]

def plan_metabolome(d):
    t1 = tier("最小可行", "先处理缺失与归一化，确认整体结构",
        ["metabolomics/01_preprocess.py：缺失过滤/插补、归一化、QC",
         "PCA 与样本相关性（找离群）",
         "报告检测到的代谢物数与缺失率"],
        ["预处理矩阵", "PCA/相关性图", "QC 报告"],
        "缺失机制（MNAR）会影响插补；需报告策略")
    t2 = tier("标准", "组间差异与通路",
        ["metabolomics/02_oplsda.R（OPLS-DA、VIP、置换检验）",
         "metabolomics/03_pathway.R（通路富集）",
         "metabolomics/05_biomarker.R（ROC/随机森林，需独立验证）"],
        ["差异代谢物表", "通路图", "标志物候选"],
        "OPLS-DA 易过拟合；必须报告置换检验与交叉验证")
    t3 = tier("进阶", "网络、分类与遗传定位",
        ["metabolomics/06_metabolite_network.R（共表达网络）",
         "metabolomics/04_metabolite_classification.R（注释分类）",
         "qtl/02_colocalization.R（有基因型时做 mQTL）"],
        ["代谢网络", "分类注释", "mQTL 候选"],
        "网络边不代表因果；mQTL 需足够样本量")
    return [t1, t2, t3]

def plan_qtl(d):
    t1 = tier("最小可行", "先确认家系结构与标记质量",
        ["qtl/01_qtl_mapping.R：遗传图谱/连锁分析（家系设计）",
         "01_plink_qc.sh：标记质控（若为芯片/重测序）",
         "表型分布与区组检查"],
        ["连锁图", "QTL 区间初表"],
        "家系结构错误会导致连锁相判断错误")
    t2 = tier("标准", "QTL 区间 + 候选基因",
        ["qtl/01_qtl_mapping.R（LOD 阈值/置换检验）",
         "04_qtl_region_snps.py（区间 SNP 提取与注释）",
         "05_bsa_scan.py（若为 BSA-seq：ΔSNP-index/G′）"],
        ["QTL 区间表", "区间基因注释", "BSA 候选区间"],
        "区间大小受重组与标记密度限制；BSA 需足够混池个体")
    t3 = tier("进阶", "共定位与候选验证",
        ["qtl/02_colocalization.R / 03_coloc_bayesian.R（与 GWAS/表达共定位）",
         "gwas/07_haplotype.R（候选基因单倍型）",
         "design/01_sgrna_design.py（编辑验证靶点）"],
        ["共定位证据", "单倍型表型", "编辑靶点清单"],
        "共定位不等于因果；需功能验证")
    return [t1, t2, t3]

def plan_gs(d):
    t1 = tier("最小可行", "先做表型与基因型的匹配检查",
        ["phenomics/01_phenotype_blup.R（BLUP/遗传力）",
         "05_power_analysis.R（功效评估）",
         "基因型质控与样本 ID 对齐"],
        ["BLUP 表型", "遗传力估计", "样本匹配表"],
        "表型与基因型样本错配是最常见错误")
    t2 = tier("标准", "GBLUP 与交叉验证",
        ["phenomics/03_gs_gblup.R（VanRaden 亲缘 + 岭回归 + k 折 CV）",
         "多环境时：phenomics/04_ammi_gge.R（稳定性）",
         "报告重复 CV 的均值±SD 与预测-观测相关"],
        ["预测准确率", "稳定性参数", "选择指数"],
        "CV 划分必须按家系/环境分组，随机划分会高估准确率")
    t3 = tier("进阶", "深度学习与多环境预测",
        ["phenomics/06_menet_gs.py（MeNet 多环境预测，需 PyTorch）",
         "多环境/多年份联合建模",
         "与 GBLUP 基线对比"],
        ["多环境预测表", "模型对比", "候选品系排序"],
        "深度模型小样本易过拟合；必须与 GBLUP 基线比较")
    return [t1, t2, t3]

def plan_pangenome(d):
    t1 = tier("最小可行", "先统一组装与注释",
        ["pangenome/01_pangenome_graph.sh（minigraph/cactus/PGGB）",
         "统一注释流程与基因 ID 映射",
         "组装质量评估（BUSCO/连续性）"],
        ["图泛基因组", "基因 ID 映射表", "组装质量表"],
        "注释不统一会让 PAV 矩阵不可比")
    t2 = tier("标准", "SV/PAV 与泛基因组关联",
        ["pangenome/02_sv_detect.sh（minimap2+SYRI）",
         "pangenome/03_pav_analysis.R（核心/可变/私有）",
         "pangenome/04_pan_gwas.R（Pan-GWAS）"],
        ["SV 集合", "PAV 矩阵", "Pan-GWAS 结果"],
        "SV 假阳性率高；需多种调用策略交叉验证")
    t3 = tier("进阶", "属级流程与注释扩展",
        ["pangenome/05_oryza_pangenome_pipeline.sh（Snakemake 阶段运行器）",
         "pangenome/06_annosine.sh（SINE 注释）",
         "与表达/表型联合"],
        ["属级泛基因组", "SINE 注释", "联合分析"],
        "上游为研究流程集合；需按 envs 建环境并逐阶段运行")
    return [t1, t2, t3]

def plan_epigenome(d):
    t1 = tier("最小可行", "先看覆盖度与重复一致性",
        ["epigenomics/01_methylation_analysis.R（WGBS）或 02_atac_chip_analysis.R（ATAC/ChIP）",
         "样本相关性/覆盖度检查",
         "批次与技术协变量整理"],
        ["覆盖度/FRiP 表", "样本相关性", "QC 报告"],
        "覆盖度不足的区域不能做差异分析")
    t2 = tier("标准", "差异区域与注释",
        ["差异甲基化区域（DMR）或差异可及区域（DAR）",
         "区域注释（启动子/增强子/基因体）",
         "与表达联合（若有 RNA-seq）"],
        ["DMR/DAR 表", "区域注释", "表达关联"],
        "区域定义与平滑参数会影响结果；需报告阈值")
    t3 = tier("进阶", "三维基因组与多组学",
        ["Hi-C/染色质互作（若有）",
         "与 GWAS/QTL 共定位",
         "多组学整合"],
        ["互作网络", "共定位", "整合因子"],
        "Hi-C 分辨率与文库质量决定可用性")
    return [t1, t2, t3]

def plan_proteome(d):
    t1 = tier("最小可行", "先处理缺失与归一化",
        ["proteomics/01_preprocess.py：缺失过滤/插补、归一化",
         "样本相关性/PCA",
         "报告缺失率与插补策略"],
        ["预处理矩阵", "PCA", "缺失报告"],
        "MNAR 缺失会让插补偏倚；缺失>50% 不做差异")
    t2 = tier("标准", "差异蛋白与功能",
        ["proteomics/02_diff_protein.R（limma）",
         "proteomics/03_plot_proteome.R（火山/热图）",
         "proteomics/04_uniprot_mapping.py（蛋白→基因）",
         "proteomics/05_interpro_domains.py（结构域富集）"],
        ["差异蛋白表", "火山/热图", "结构域富集"],
        "蛋白丰度不等于 mRNA 丰度；结论限定在蛋白层")
    t3 = tier("进阶", "翻译后修饰与多组学",
        ["磷酸化等 PTM 分析（需 PTM 数据）",
         "与转录组联合",
         "通路/网络整合"],
        ["PTM 位点表", "联合分析", "网络"],
        "PTM 数据需专门的富集与归一化")
    return [t1, t2, t3]

# --------------------------------------------------------------------------- #
# 文献依据（仅白名单期刊）
# --------------------------------------------------------------------------- #
SOURCES = {
 "01_qc_filter.py": ["Love et al. Genome Biology 2014 (DESeq2, doi:10.1186/s13059-014-0550-8)"],
 "02_deseq2.R": ["Love et al. Genome Biology 2014 (doi:10.1186/s13059-014-0550-8)",
                 "Robinson et al. Bioinformatics 2010 (edgeR, doi:10.1093/bioinformatics/btp616)",
                 "Law et al. Genome Biology 2014 (limma-voom, doi:10.1186/gb-2014-15-2-r29)"],
 "03_enrichment.R": ["Subramanian et al. PNAS 2005 (GSEA framework, doi:10.1073/pnas.0506580102)"],
 "04_gsea.R": ["Subramanian et al. PNAS 2005 (doi:10.1073/pnas.0506580102)"],
 "01_plink_qc.sh": ["Chang et al. GigaScience 2015 (PLINK 1.9, doi:10.1186/s13742-015-0047-8)"],
 "09_multilocus_qtn.R": ["Li et al. Molecular Plant 2022 (3VmrMLM, doi:10.1016/j.molp.2022.02.012)",
                          "Plant Communications 2025 (Fast3VmrMLM, doi:10.1016/j.xplc.2025.101385)"],
 "02_colocalization.R": ["Genome Biology 2025 (rice mGWAS + colocalization, doi:10.1186/s13059-025-03513-w)"],
 "03_coloc_bayesian.R": ["Genome Biology 2025 (rice mGWAS + colocalization, doi:10.1186/s13059-025-03513-w)"],
 "01_mofa.R": ["Argelaguet et al. Molecular Systems Biology 2018 (MOFA, doi:10.15252/msb.20178157)"],
 "02_diablo.R": ["Singh et al. Bioinformatics 2019 (DIABLO, doi:10.1093/bioinformatics/btz158)"],
 "06_menet_gs.py": ["Li et al. Plant Communications 2025 (MeNet, doi:10.1016/j.xplc.2025.101620)"],
 "metabolomics/01_preprocess.py": ["Genome Biology 2025 (rice hyperspectral-metabolite mGWAS, doi:10.1186/s13059-025-03513-w)"],
}

def sources_for(steps):
    out = []
    for s in steps:
        for key, refs in SOURCES.items():
            if key in s:
                for r in refs:
                    if r not in out:
                        out.append(r)
    return out


# --------------------------------------------------------------------------- #
# 图形方案：每档给出"图名 / 类型 / 生成脚本 / 规格"
# 规格默认 CNS 风格：单栏 85 mm、双栏 170 mm、300 dpi、PDF+SVG+TIFF
# --------------------------------------------------------------------------- #
FIG_SPEC = "85/170 mm, 300 dpi, PDF+SVG+TIFF"
FIGURES = {
 "rnaseq": {
  "最小可行": [("PCA 散点", "散点", "01_qc_filter.py"), ("样本相关性热图", "热图", "01_qc_filter.py")],
  "标准": [("PCA 散点", "散点", "01_qc_filter.py"), ("样本相关性热图", "热图", "01_qc_filter.py"),
           ("差异火山图", "火山", "07_plot_transcriptome.R"), ("差异基因热图", "热图", "07_plot_transcriptome.R"),
           ("GO/KEGG 富集气泡图", "气泡", "03_enrichment.R")],
  "进阶": [("差异火山图", "火山", "07_plot_transcriptome.R"), ("GSEA 富集曲线", "曲线", "04_gsea.R"),
           ("WGCNA 模块-性状热图", "热图", "05_wgcna.R"), ("基因-代谢物网络", "网络", "05_network_integration.R")],
 },
 "gwas": {
  "最小可行": [("QQ 图", "QQ", "05_plot_gwas.R"), ("曼哈顿图", "曼哈顿", "05_plot_gwas.R")],
  "标准": [("QQ 图", "QQ", "05_plot_gwas.R"), ("曼哈顿图", "曼哈顿", "05_plot_gwas.R"),
           ("候选区间 LD 热图", "热图", "04_candidate_loci.R"), ("单倍型-表型箱线图", "箱线", "07_haplotype.R")],
  "进阶": [("多位点 QTN 曼哈顿图", "曼哈顿", "09_multilocus_qtn.R"), ("QEI 环境互作热图", "热图", "09_multilocus_qtn.R"),
           ("共定位 PP4 图", "柱状", "02_colocalization.R")],
 },
 "multiomics": {
  "最小可行": [("各层 PCA", "散点", "metabolomics/01_preprocess.py"), ("样本匹配示意", "示意", "01_make_plan.py")],
  "标准": [("OPLS-DA 得分图", "散点", "metabolomics/02_oplsda.R"), ("VIP 条形图", "条形", "metabolomics/02_oplsda.R"),
           ("跨层相关热图", "热图", "02_diablo.R")],
  "进阶": [("MOFA 因子图", "散点", "01_mofa.R"), ("基因-代谢物网络", "网络", "05_network_integration.R"),
           ("共定位 PP4 图", "柱状", "02_colocalization.R")],
 },
 "singlecell": {
  "最小可行": [("每样本 QC 小提琴图", "小提琴", "01_scanpy_pipeline.py"), ("整合 UMAP", "UMAP", "01_scanpy_pipeline.py")],
  "标准": [("UMAP + 细胞类型", "UMAP", "01_scanpy_pipeline.py"), ("marker 点图", "点图", "01_scanpy_pipeline.py"),
           ("细胞比例堆叠图", "堆叠", "01_scanpy_pipeline.py")],
  "进阶": [("轨迹图", "轨迹", "01_scanpy_pipeline.py"), ("调控模块热图", "热图", "01_scanpy_pipeline.py"),
           ("空间映射图", "空间", "01_spatial_pipeline.py")],
 },
 "spatial": {
  "最小可行": [("空间域图", "空间", "01_spatial_pipeline.py"), ("QC 散点", "散点", "01_spatial_pipeline.py")],
  "标准": [("空间域图", "空间", "01_spatial_pipeline.py"), ("区域 marker 空间图", "空间", "01_spatial_pipeline.py"),
           ("区域差异热图", "热图", "01_spatial_pipeline.py")],
  "进阶": [("邻域网络图", "网络", "01_spatial_pipeline.py"), ("多切片比较图", "空间", "01_spatial_pipeline.py")],
 },
 "metabolome": {
  "最小可行": [("PCA 散点", "散点", "01_preprocess.py"), ("缺失率条形图", "条形", "01_preprocess.py")],
  "标准": [("OPLS-DA 得分图", "散点", "02_oplsda.R"), ("VIP 条形图", "条形", "02_oplsda.R"),
           ("通路气泡图", "气泡", "03_pathway.R"), ("ROC 曲线", "曲线", "05_biomarker.R")],
  "进阶": [("代谢共表达网络", "网络", "06_metabolite_network.R"), ("mQTL 曼哈顿图", "曼哈顿", "02_colocalization.R")],
 },
 "qtl": {
  "最小可行": [("LOD 曲线", "曲线", "01_qtl_mapping.R"), ("表型分布图", "直方", "01_qtl_mapping.R")],
  "标准": [("LOD 曲线", "曲线", "01_qtl_mapping.R"), ("QTL 区间图", "区间", "01_qtl_mapping.R"),
           ("BSA ΔSNP-index 图", "曲线", "05_bsa_scan.py")],
  "进阶": [("共定位图", "柱状", "02_colocalization.R"), ("单倍型箱线图", "箱线", "07_haplotype.R")],
 },
 "gs": {
  "最小可行": [("表型分布", "直方", "01_phenotype_blup.R"), ("遗传力条形图", "条形", "01_phenotype_blup.R")],
  "标准": [("预测-观测散点", "散点", "03_gs_gblup.R"), ("CV 准确率箱线图", "箱线", "03_gs_gblup.R"),
           ("AMMI 双标图", "双标", "04_ammi_gge.R")],
  "进阶": [("MeNet 预测散点", "散点", "06_menet_gs.py"), ("模型对比图", "条形", "06_menet_gs.py")],
 },
 "pangenome": {
  "最小可行": [("泛基因组增长曲线", "曲线", "01_pangenome_graph.sh"), ("组装质量条形图", "条形", "01_pangenome_graph.sh")],
  "标准": [("PAV 热图", "热图", "03_pav_analysis.R"), ("核心/可变基因饼图", "饼图", "03_pav_analysis.R"),
           ("Pan-GWAS 曼哈顿图", "曼哈顿", "04_pan_gwas.R")],
  "进阶": [("SV Circos 图", "Circos", "02_sv_detect.sh"), ("SINE 分布图", "分布", "06_annosine.sh")],
 },
 "epigenome": {
  "最小可行": [("样本相关性热图", "热图", "01_methylation_analysis.R"), ("覆盖度分布图", "分布", "01_methylation_analysis.R")],
  "标准": [("DMR/DAR 火山图", "火山", "01_methylation_analysis.R"), ("区域注释饼图", "饼图", "01_methylation_analysis.R"),
           ("metagene 轮廓图", "曲线", "01_methylation_analysis.R")],
  "进阶": [("Hi-C 互作热图", "热图", "02_atac_chip_analysis.R"), ("共定位图", "柱状", "02_colocalization.R")],
 },
 "proteome": {
  "最小可行": [("PCA 散点", "散点", "01_preprocess.py"), ("缺失率条形图", "条形", "01_preprocess.py")],
  "标准": [("差异蛋白火山图", "火山", "02_diff_protein.R"), ("差异蛋白热图", "热图", "03_plot_proteome.R"),
           ("结构域富集条形图", "条形", "05_interpro_domains.py")],
  "进阶": [("PTM 位点图", "分布", "02_diff_protein.R"), ("蛋白-mRNA 一致性散点", "散点", "01_preprocess.py")],
 },
}

# --------------------------------------------------------------------------- #
# 路由与渲染
# --------------------------------------------------------------------------- #
def build(d):
    b, w = validate(d)
    assay, goal, design = d["assay"], d["goal"], d["design"]
    if goal == "qtl" or design == "ril_magic_nam":
        tiers = plan_qtl(d)
    elif goal == "prediction":
        tiers = plan_gs(d)
    elif assay == "single_cell":
        tiers = plan_singlecell(d)
    elif assay == "spatial":
        tiers = plan_spatial(d)
    elif assay == "pangenome":
        tiers = plan_pangenome(d)
    elif assay in ("wgbs", "atac"):
        tiers = plan_epigenome(d)
    elif assay == "proteome":
        tiers = plan_proteome(d)
    elif assay == "metabolome" and goal in ("marker_discovery", "differential_expression"):
        tiers = plan_metabolome(d)
    elif goal in ("differential_expression", "marker_discovery") and assay in ("bulk_rna", "rna_metabolome"):
        tiers = plan_rnaseq(d)
    elif goal == "gwas" or assay == "resequencing":
        tiers = plan_gwas(d)
    elif goal == "integration" or assay == "rna_metabolome":
        tiers = plan_multiomics(d)
    else:
        tiers = plan_rnaseq(d)
    key = ("qtl" if goal == "qtl" or design == "ril_magic_nam" else
           "gs" if goal == "prediction" else
           "singlecell" if assay == "single_cell" else
           "spatial" if assay == "spatial" else
           "pangenome" if assay == "pangenome" else
           "epigenome" if assay in ("wgbs", "atac") else
           "proteome" if assay == "proteome" else
           "metabolome" if assay == "metabolome" and goal in ("marker_discovery", "differential_expression") else
           "rnaseq" if goal in ("differential_expression", "marker_discovery") and assay in ("bulk_rna", "rna_metabolome") else
           "gwas" if goal == "gwas" or assay == "resequencing" else
           "multiomics")
    for t in tiers:
        t["figures"] = [{"name": n, "type": ty, "script": sc, "spec": FIG_SPEC}
                        for n, ty, sc in FIGURES.get(key, {}).get(t["tier"], [])]
    return {"input": d, "blockers": b, "warnings": w,
            "tiers": tiers, "generated": str(date.today())}

def render_md(r):
    d = r["input"]
    L = ["# 分析方案（自动生成）", "",
         f"- 日期：{r['generated']}",
         f"- 数据类型：`{d['assay']}`｜目标：`{d['goal']}`｜设计：`{d['design']}`",
         f"- 样本量：{d['n_per_group']}｜重复：{d['replicates']}｜批次：{d['batch']}"
         f"｜环境数：{d.get('environments', 1)}", ""]
    if r["blockers"]:
        L += ["## ⛔ 阻断项（先解决再分析）", ""]
        L += [f"- {x}" for x in r["blockers"]] + [""]
    if r["warnings"]:
        L += ["## ⚠️ 风险提示", ""]
        L += [f"- {x}" for x in r["warnings"]] + [""]
    for t in r["tiers"]:
        L += [f"## {t['tier']}：{t['intent']}", "", "**步骤**"]
        L += [f"1. {s}" for s in t["steps"]]
        L += ["", "**预期产出**：" + "；".join(t["outputs"])]
        if t.get("figures"):
            L += ["", "**图形方案**"]
            for f in t["figures"]:
                L.append(f"- {f['name']}（{f['type']}，{f['spec']}）— `{f['script']}`")
        L += ["", f"**风险与边界**：{t['risk']}"]
        refs = sources_for(t["steps"])
        if refs:
            L += ["", "**方法依据**"]
            L += [f"- {x}" for x in refs]
        L.append("")
    L += ["---", "说明：本方案由 `scripts/plan/01_make_plan.py` 依据数据契约生成；"
          "所有脚本均在本仓库 `scripts/` 下，可先跑 `--demo` 验证环境。"]
    return "\n".join(L)

def _script_index():
    idx = {}
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if fn.endswith((".R", ".py", ".sh")):
                rel = os.path.relpath(os.path.join(dp, fn), os.path.dirname(root)).replace("\\", "/")
                idx.setdefault(fn, []).append(rel)
    return idx

def _resolve(token, idx):
    if "/" in token:
        cand = "scripts/" + token.split("scripts/")[-1]
        return cand if os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", cand)) else None
    hits = idx.get(token, [])
    return hits[0] if len(hits) == 1 else None

def render_commands(r):
    idx = _script_index()
    L = ["#!/usr/bin/env bash", "# 由 01_make_plan.py 生成：先跑 --demo 验证环境，再按 docs/*_usage.md 传真实参数",
         "set -euo pipefail", ""]
    seen = []
    for t in r["tiers"]:
        L.append(f"# --- {t['tier']} ---")
        for s in t["steps"]:
            for token in re.findall(r"[A-Za-z0-9_./-]+\.(?:R|py|sh)", s):
                token = token.lstrip("./")
                if token in seen:
                    continue
                seen.append(token)
                path = _resolve(token, idx)
                if path:
                    L.append(f"# bash {path} --demo")
                else:
                    L.append(f"# bash scripts/**/{token} --demo   # 名称不唯一，按模块选择")
        L.append("")
    L.append("# 真实参数见 docs/installation.md 与 docs/*_usage.md")
    return "\n".join(L)

def run(d, outdir):
    r = build(d)
    os.makedirs(outdir, exist_ok=True)
    open(os.path.join(outdir, "plan.md"), "w", encoding="utf-8").write(render_md(r))
    open(os.path.join(outdir, "plan.json"), "w", encoding="utf-8").write(
        json.dumps(r, ensure_ascii=False, indent=2))
    open(os.path.join(outdir, "commands.sh"), "w", encoding="utf-8").write(render_commands(r))
    return r

DEMOS = {
 "demo_rnaseq_ok": {"assay": "bulk_rna", "goal": "differential_expression",
    "design": "treatment_control", "n_per_group": 3, "replicates": True,
    "factors": ["treatment"], "batch": False, "reference": "IRGSP-1.0"},
 "demo_rnaseq_n1": {"assay": "bulk_rna", "goal": "differential_expression",
    "design": "treatment_control", "n_per_group": 1, "replicates": False,
    "factors": ["treatment"], "batch": False, "reference": "IRGSP-1.0"},
 "demo_gwas": {"assay": "resequencing", "goal": "gwas", "design": "natural_population",
    "n_per_group": 300, "replicates": False, "phenotype": "binary",
    "environments": 2, "reference": "IRGSP-1.0"},
 "demo_multiomics": {"assay": "rna_metabolome", "goal": "integration",
    "design": "treatment_control", "n_per_group": 60, "replicates": True,
    "factors": ["treatment"], "batch": True, "reference": "IRGSP-1.0"},
 "demo_singlecell": {"assay": "single_cell", "goal": "differential_expression",
    "design": "treatment_control", "n_per_group": 4, "replicates": True,
    "factors": ["treatment"], "batch": True, "reference": "IRGSP-1.0"},
 "demo_spatial": {"assay": "spatial", "goal": "differential_expression",
    "design": "treatment_control", "n_per_group": 2, "replicates": True,
    "factors": ["treatment"], "batch": False, "reference": "IRGSP-1.0"},
 "demo_metabolome": {"assay": "metabolome", "goal": "marker_discovery",
    "design": "treatment_control", "n_per_group": 6, "replicates": True,
    "factors": ["treatment"], "batch": False, "quantification": "untargeted"},
 "demo_qtl": {"assay": "resequencing", "goal": "qtl", "design": "ril_magic_nam",
    "n_per_group": 200, "replicates": False, "phenotype": "continuous",
    "reference": "IRGSP-1.0"},
 "demo_gs": {"assay": "resequencing", "goal": "prediction",
    "design": "natural_population", "n_per_group": 400, "replicates": False,
    "phenotype": "continuous", "environments": 2, "reference": "IRGSP-1.0"},
 "demo_pangenome": {"assay": "pangenome", "goal": "gwas", "design": "natural_population",
    "n_per_group": 12, "replicates": False, "phenotype": "continuous",
    "reference": "IRGSP-1.0"},
 "demo_proteome": {"assay": "proteome", "goal": "differential_expression",
    "design": "treatment_control", "n_per_group": 4, "replicates": True,
    "factors": ["treatment"], "batch": False},
 "demo_wgbs": {"assay": "wgbs", "goal": "differential_expression",
    "design": "treatment_control", "n_per_group": 4, "replicates": True,
    "factors": ["treatment"], "batch": False, "reference": "IRGSP-1.0"},
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input"); ap.add_argument("--outdir", default="results/plan")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        base = os.path.join(a.outdir, "demo"); ok = 0
        for name, d in DEMOS.items():
            r = run(normalize(dict(d)), os.path.join(base, name))
            print(f"[demo] {name}: tiers={len(r['tiers'])} "
                  f"blockers={len(r['blockers'])} warnings={len(r['warnings'])}")
            if name in ("demo_rnaseq_n1", "demo_spatial_n1"):
                assert r["blockers"], f"{name} 应触发阻断项"
            ok += 1
        r = build(normalize({"assay": "single_cell", "goal": "differential_expression",
                             "design": "treatment_control", "n_per_group": 1,
                             "replicates": False}))
        assert r["blockers"], "single-cell n=1 必须阻断"
        print(f"DEMO OK ({ok}/{len(DEMOS)}) + guardrail -> {base}")
        return
    if not a.input:
        sys.exit("请提供 --input data.json（或使用 --demo）")
    r = run(load(a.input), a.outdir)
    print(f"blockers={len(r['blockers'])} warnings={len(r['warnings'])} -> {a.outdir}/plan.md")

if __name__ == "__main__":
    main()
