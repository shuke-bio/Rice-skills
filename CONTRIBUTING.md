# Contributing to Rice-skills

Contributions are welcome — bug reports, module improvements, new analysis
scripts, documentation fixes, and translations.

贡献欢迎：问题报告、模块改进、新分析脚本、文档修正与翻译。

## Before you start

- **Open an issue first** for anything larger than a typo: describe the
  module/scenario, expected behavior, and your data format.
- 较大的改动请先开 issue 讨论，描述模块/场景、期望行为与数据格式。

## Ground rules

1. **Reproducibility first.** Every analysis script must ship with a
   `--demo` (or equivalent) self-test on bundled demo data, and must not require
   network access to pass it.
2. **License hygiene.** Do not copy code, data, or assets from repositories
   without a clear open-source license (MIT / BSD / Apache-2.0 / CC BY 4.0 are
   acceptable with attribution + NOTICE). GPL-family components may only be
   bundled as an isolated directory preserving their own license. When in
   doubt, open an issue first.
3. **Figure style.** Use the shared theme (`config/nature_palette.json` +
   `scripts/utils/plot_theme.*`) — do not hardcode palettes or fonts.
4. **Keep the module contract.** Each script: CLI flags via `--help`,
   `run_info.txt` provenance file in its output directory, no silent
   write-outs outside `--outdir`.

## Pull request checklist

- [ ] `--demo` self-test passes locally
- [ ] New/changed script documented in the relevant `docs/*_usage.md`
- [ ] Counters updated (README / SKILL.md) if the script inventory changed
- [ ] No personal paths, credentials, or third-party files without license
- [ ] Commit messages describe intent (English preferred)

## Reporting bugs

Include: the exact command line, the module, `run_info.txt`, the error output,
and (if possible) the smallest demo dataset reproducing the issue.
Attaching the `--demo` command that fails is the fastest path to a fix.

Other questions: open a [GitHub issue](https://github.com/shuke-bio/Rice-skills/issues),
For anything that should not be public, open an issue without sensitive details and we will follow up.
其他问题：开 [GitHub issue](https://github.com/shuke-bio/Rice-skills/issues)，
不便公开的内容，请开 Issue 并避免贴敏感信息，我们会另行跟进。
