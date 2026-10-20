# 安装指南（Installation）

本文说明如何搭建 Rice-skills 的运行环境：conda 一键安装（推荐）、手工安装、以及常见故障排查。

## 0. 三条安装路线（先选一条）

| 路线 | 适用 | 命令 / 操作 |
|---|---|---|
| A. skills CLI | Codex / Claude Code / OpenClaw 等支持 skills 的 agent | `npx skills add shuke-bio/Rice-skills --skill rice-skills --yes --copy` |
| B. Git 克隆 | 需要本地脚本库、要改代码 | `git clone https://github.com/shuke-bio/Rice-skills.git` |
| C. 代理指向 / ZIP | WorkBuddy、DeepSeek Harness 等 | 把仓库链接发给 agent，说“从这个仓库安装水稻技能”；或下载 ZIP 解压到技能目录 |

### 0.1 WorkBuddy（Windows）

1. 从 workbuddy.cn 下载 Windows 版并安装、登录；
2. 把克隆/解压后的 `Rice-skills` 文件夹放入 WorkBuddy 的技能目录，或在对话中
   直接说“从 github.com/shuke-bio/Rice-skills 安装 rice 技能”；
3. 安装 R 4.4+ 与 Python 3.11+（安装时勾选加入 PATH），再按第 2 节创建 conda 环境；
4. 验证：`Rscript scripts/gwas/09_multilocus_qtn.R --demo`。

### 0.2 DeepSeek Harness（dsh）

1. 安装 Node.js 22.19+（或 24 LTS），在终端运行 `npx @deepseek-ai/dsh web`，
   浏览器打开 `http://127.0.0.1:3080`；
2. 在 platform.deepseek.com 创建 API Key，进入 设置 → Models 粘贴并选择模型；
3. 选择工作区（如 `D:\work\rice`），在对话中让它把本仓库安装为技能；
4. 验证：请它运行 `Rscript scripts/phenomics/06_menet_gs.py --demo`（缺依赖会给出安装提示）。

> 说明：dsh 仍为开发者预览版，命令与模型名以官方文档为准；API Key 保存在
> `~/.dsh/.credentials.yaml`，不要提交到仓库。

## 1. 前置要求


| 组件 | 版本 | 说明 |
|---|---|---|
| conda（mamba/miniconda/anaconda 均可） | 最新 | 推荐用 [miniforge](https://github.com/conda-forge/miniforge) |
| R | ≥ 4.3 | 由 conda 自动安装，无需单独装 |
| Python | ≥ 3.10 | 由 conda 自动安装 |
| 磁盘 | ≥ 8 GB | Bioconductor 包体积较大 |
| 网络 | 可访问 Bioconductor / CRAN / KEGG | 富集分析需联网访问 KEGG 与 GO |

> Windows 用户请在 **Git Bash** 中执行所有 `bash` 命令（本项目的 `.sh` 脚本为 POSIX shell）。
> 若已在系统安装了 R，可用 `--rscript "C:/Program Files/R/R-4.6.0/bin/Rscript.exe"` 指定路径。

## 2. 一键安装（推荐）

```bash
cd Rice-skills
conda env create -f environment.yml     # 约 20-40 分钟
conda activate rice-skills
```

`environment.yml` 已包含：Python 3.11、R 4.3、R 端核心包（ggplot2、WGCNA、clusterProfiler、DESeq2、MOFA2、mixOmics 等）以及 `requirements.txt` 中的 Python 包。

### 2.1 ropls 的补充安装

`ropls`（OPLS-DA）在部分渠道不可用，环境创建后补装：

```bash
conda activate rice-skills
Rscript -e 'install.packages("ropls", repos = "https://cloud.r-project.org")'
```

### 2.2 验证安装

```bash
bash tests/test_pipeline.sh
# 若默认 python/Rscript 不在 PATH：
bash tests/test_pipeline.sh --python /path/to/python --rscript /path/to/Rscript
```

该脚本会检查：仓库结构、Python 语法、R 语法（20 个文件）、R/Python 工具层冒烟测试、以及示例数据端到端跑通。全部 PASS = 环境就绪。

## 3. 手工安装（不想用 conda）

### 3.1 Python 依赖

```bash
pip install -r requirements.txt
```

### 3.2 R 依赖（Bioconductor）

```bash
Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'

Rscript -e 'BiocManager::install(c(
  "DESeq2", "edgeR", "clusterProfiler", "org.Osativa.eg.db", "ReactomePA",
  "fgsea", "WGCNA", "MOFA2", "mixOmics", "ComplexHeatmap", "MetaboAnalystR",
  "imputeLCMD"
))'

Rscript -e 'install.packages(c(
  "ropls", "ggplot2", "ggrepel", "igraph", "ggraph", "pheatmap", "jsonlite",
  "UpSetR", "ggalluvial", "circlize", "patchwork", "pROC", "randomForest",
  "dplyr", "readr", "tidyr", "stringr", "ragg"
), repos = "https://cloud.r-project.org")'
```

> `org.Osativa.eg.db` 是水稻（Japonica）注释库。若研究其他物种，替换为目标物种的
> `org.*.eg.db`，并把脚本的 `--org` 参数改为对应 KEGG 代码（如玉米 `zma`、拟南芥 `ath`）。

## 4. 可选依赖与降级行为

本项目在可选包缺失时**不会崩溃**，而是跳过对应步骤并给出提示：

| 缺失包 | 受影响步骤 | 行为 |
|---|---|---|
| `patchwork` | 07 组合图 | 各面板单独保存，不合成 |
| `ComplexHeatmap` | 热图类 | 自动降级到 `pheatmap` |
| `ggraph` | 网络图 | 降级到 igraph 基础图形（同尺寸三格式输出） |
| `UpSetR` / `ggalluvial` / `circlize` | 整合图 | 跳过对应面板，日志提示 |
| `MOFA2` / `mixOmics` / `ropls` / `DESeq2` | 核心步骤 | `need_pkgs()` 友好报错并给出安装命令 |

## 5. 故障排查

**① `BiocManager::install` 编译失败（尤其 WGCNA/preprocessCore）**
- Windows：安装 [Rtools](https://cran.r-project.org/bin/windows/Rtools/) 并加入 PATH。
- Linux：装 `build-essential libcurl4-openssl-dev libxml2-dev libssl-dev`。
- 优先用 conda 安装的二进制包（`environment.yml`）避免编译。

**② `MOFA2` 安装失败（HDF5 相关）**
- 优先用 conda：`conda install -c bioconda bioconductor-mofa2`。
- 失败时跳过整合模块 01（其余步骤不受影响）。

**③ R 启动时出现 `Setting LC_COLLATE=C.UTF-8 failed`**
- 仅是区域设置警告，不影响分析。Windows 可在 R 中设置 `Sys.setlocale("LC_ALL", "English_United States.1252")`。

**④ 出图报 `invalid font type` 或字体缺失**
- 主题层默认使用 Helvetica（PDF 基础字体，任何平台可嵌入）；如需 Arial，请在系统中注册为 Type 1 字体后由 `plot_theme.R` 自动采用。

**⑤ 富集分析无结果 / 通路为空**
- 检查 `org.Osativa.eg.db` 是否安装；
- 基因 ID 是否被正确识别：`LOC_Os*` 默认按 `SYMBOL` 处理，换个来源的 ID 请加 `--id-type ENTREZID`；
- KEGG 富集需要联网，公司/校园网络可能需代理。

**⑥ 样本不匹配报错**
- 矩阵列名必须与 `metadata.csv` 的 `sample_id` 完全一致（大小写、空格敏感）。脚本会打印两侧样例 ID 便于核对。

**⑦ 内存不足（WGCNA / MOFA）**
- 减小 `--top-n`（默认 10000）；增大 `--min-module-size`；WGCNA 已使用 `blockwiseModules` 分块计算。
- 2 万基因 × 100 样本建议 ≥16 GB RAM。

## 6. 目录权限

所有脚本把结果写到 `--outdir`（默认 `results/`），不会修改仓库源码。`results/` 已在 `.gitignore` 中忽略，可安全删除重跑。
