# 决策层（Decision Layer）

> 回答的不是"这个脚本怎么跑"，而是"**这批数据该做什么、先做什么、不能做什么**"。

## 目录

| 文件 | 用途 |
|---|---|
| `data_contract.md` | 10 个必答问题 + 三条红线 + JSON 模板 |
| `scenario_rnaseq.md` | bulk RNA-seq：三种设计、三档方案、失败模式 |
| `scenario_gwas.md` | 自然群体/家系/多亲本：模型选择、三档方案、失败模式 |
| `scenario_multiomics.md` | 转录组×代谢组：可整合性判断、三档方案、失败模式 |
| `scenario_singlecell.md` | 单细胞与空间：样本级重复、整合、比例/伪 bulk |
| `scenario_metabolome.md` | 代谢组：靶向/非靶向、组成型数据、标志物 |
| `scenario_qtl.md` | QTL/BSA：家系 vs 混池、区间到候选 |
| `scenario_gs.md` | 表型组与 GS：BLUP、GBLUP、MeNet、交叉验证 |
| `scenario_pangenome.md` | 泛基因组：图构建、SV/PAV、Pan-GWAS |
| `scenario_epigenome.md` | WGBS/ATAC/ChIP：覆盖度、差异区域、3D |
| `scenario_proteome.md` | 蛋白组：缺失机制、差异蛋白、PTM |
| `router.md` | 用户说法 → 该问什么 → 用哪张卡 |
| `en/` | 10 张英文场景卡（与中文版一一对应） |
| `dependencies.md` | 依赖体检与“同意后安装”流程（含安全边界） |
| `figure_plan.md` | 分析+图形方案确认单、图形规格、批准与执行记录 |
| `template.json` | 数据契约模板（含路径字段） |
| `../../scripts/plan/01_make_plan.py` | 读 JSON → 输出 `plan.md` / `plan.json` / `commands.sh` |
| `../../scripts/plan/02_check_data.py` | 数据体检：缺失/文库/PCA/批次混杂/多组学重叠 |
| `../../scripts/plan/03_validate_cards.py` | 校验方法卡字段、脚本存在性、期刊白名单 |
| `../../scripts/plan/04_check_deps.py` | 依赖体检：R/Python 包与外部工具；同意后才安装 |
| `../../scripts/plan/05_run_plan.py` | 方案确认闸门：无 approval.json 拒绝执行 |

## 使用顺序

1. 按 `data_contract.md` 填一份 JSON（可用 `template.json` 复制）；
2. 有数据文件时先体检：`python scripts/plan/02_check_data.py --matrix ... --metadata ...`；
3. 生成方案：`python scripts/plan/01_make_plan.py --input my_data.json --outdir results/plan`；
4. 先看 **阻断项**——有阻断项就不要往下跑；
5. 在"最小可行 / 标准 / 进阶"里选一档；
6. **确认闸门**：`python scripts/plan/05_run_plan.py --plan results/plan/plan.json --tier minimal`
   把分析步骤与图形方案给用户看；确认后 `--approve`，再 `--execute`；
7. 依赖体检：`python scripts/plan/04_check_deps.py --plan results/plan/plan.json --tier minimal`；
8. 把缺失清单交给用户确认，同意后加 `--install --yes` 安装并复检；
9. 按 `commands.sh` 逐个 `--demo` 验证环境，再上真实数据。

## 文献准入标准

- 只引用 **A1–A3 白名单期刊**：CNS 及子刊、PNAS、Molecular Plant、Plant Cell、
  Plant Communications、Plant Physiology、JIPB、New Phytologist、Plant Journal、
  Briefings in Bioinformatics、NAR、Genome Biology、Genome Research、Bioinformatics、
  GigaScience、Plant Phenomics、The Crop Journal、Rice；
- 方法学原始文献若不在白名单（如部分工具的开源论文），**只保留工具名**，不列为
  "方法依据"，避免把非白名单期刊当作设计依据；
- 每张场景卡的"高分论文的设计启发"只取水稻相关、且来自上述期刊的工作。

## 与其它目录的分工

- `docs/references/literature_methods_2022_2026.md`：方法库（70 条，S/A/B 分级）；
- `docs/decision/`：决策层（能不能做、先做哪一步）；
- `docs/*_usage.md`：执行层（脚本参数与输出）。
