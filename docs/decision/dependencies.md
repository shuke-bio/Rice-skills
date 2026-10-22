# 依赖体检与安装（Dependency Preflight）

> 解决的问题：用户提出分析需求时，本机缺 R 包 / Python 包 / 外部工具。
> 原则：**先体检、明示清单、用户同意、再安装、装完复检、最后执行**；默认不做任何改动。

## 1. 固定流程（agent 必须按此顺序）

1. 生成方案：`python scripts/plan/01_make_plan.py --input data.json --outdir results/plan`
2. 用户选档（最小可行 / 标准 / 进阶）
3. **依赖体检**：
   ```bash
   python scripts/plan/04_check_deps.py --plan results/plan/plan.json --tier minimal
   ```
4. 把 `deps_report.md` 的缺失清单与安装命令展示给用户，**明确询问是否同意安装**；
5. 用户同意后才执行：
   ```bash
   python scripts/plan/04_check_deps.py --plan results/plan/plan.json --tier minimal --install --yes
   ```
6. 安装后脚本自动复检；仍缺失则如实报告并停止；
7. 复检通过后再运行分析（可加 `--run` 以 `--demo` 逐个验证）。

## 2. 会安装什么、不会安装什么

| 类型 | 是否自动安装 | 说明 |
|---|---|---|
| Python 包 | ✅（用户同意后） | `python -m pip install`，装到当前解释器 |
| R 包（CRAN） | ✅（用户同意后） | `install.packages()` |
| R 包（Bioconductor） | ✅（用户同意后） | 先装 BiocManager，再 `BiocManager::install()` |
| 系统工具（PLINK/BWA/BLAST/TRF/Snakemake…） | ❌ 只给指引 | 需 conda/官网安装，涉及系统权限，不自动执行 |
| 模型权重（ELMERF/Chloro-Count） | ❌ 只给链接 | 由原作者分发，仓库不转存 |

## 3. 安全边界

- **默认只检测**：不加 `--install` 不会执行任何安装；
- `--install` 仍需 `--yes` 或交互确认；非交互环境读不到确认时**退出码 2**，不做改动；
- 不调用 `sudo`/管理员权限，不修改系统 PATH；
- 安装日志写入 `results/deps/install.log`，便于排错；
- 安装失败不会被吞掉：脚本返回非零并列出仍缺失的包。

## 4. 输出文件

| 文件 | 内容 |
|---|---|
| `deps_report.md` | 缺失清单、安装方式、可复制的安装命令 |
| `deps.json` | 结构化结果（脚本列表、缺失项、命令），供 agent 读取 |
| `install.log` | 仅安装时生成，完整 stdout/stderr |

## 5. 依赖是怎么识别的

- **Python**：扫描 `import X` / `from X import`，过滤标准库与仓库内本地模块；
- **R**：扫描 `library()` / `require()` / `requireNamespace()` / `need_pkgs(c(...))` / `pkg::fun`；
  纯小写且带下划线的名字（如 `org_pkg`）视为变量，不当作包；
- **Shell**：匹配 PLINK、BWA、samtools、GATK、minimap2、SyRI、Snakemake、BLAST+、TRF 等工具名。

## 6. 常见问题

| 现象 | 原因 | 处理 |
|---|---|---|
| 装了包仍提示缺失 | 安装到了另一个 R/Python | 用 `--python` / `--rscript` 指定实际解释器 |
| Bioconductor 包版本不匹配 | BiocManager 版本与 R 不匹配 | 升级 BiocManager，或按 Bioc 官方对应版本安装 |
| R 包编译失败（Windows） | 缺 Rtools | 安装 Rtools 并重启终端 |
| 公司网络无法访问 CRAN/Bioc | 代理限制 | 配置镜像（`options(repos=...)`）或离线安装包 |
| 可选包也被列出 | 脚本中为条件加载 | 可跳过；不影响主流程 |

## 7. 示例（本机实测）

`scripts/transcriptomics/02_deseq2.R` 在本机检测结果：

```
[R ] DESeq2
[R ] SummarizedExperiment
```

安装命令（同意后执行）：

```bash
Rscript -e if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org");BiocManager::install(c('DESeq2','SummarizedExperiment'),ask=FALSE,update=FALSE)
```
