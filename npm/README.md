# rice-skills (npm installer)

Installer for [Rice-skills](https://github.com/shuke-bio/Rice-skills) — an open-source
rice multi-omics and breeding analysis toolkit (18 modules, 91 R/Python scripts, MIT).

The npm package contains **only the installer**; analysis scripts are downloaded from the
pinned GitHub tag, so the package stays small and every install is version-explicit.

## Usage

```bash
# install the latest packaged version into ./Rice-skills
npx rice-skills install

# choose a target directory (e.g. your agent's skills folder)
npx rice-skills install --dir ~/.workbuddy/skills/Rice-skills

# install a specific tag into a non-empty directory
npx rice-skills install --dir ./Rice-skills --ref v3.3.9 --force

# verify an installation
npx rice-skills check --dir ./Rice-skills
```

## What it does

1. Downloads `https://codeload.github.com/shuke-bio/Rice-skills/tar.gz/refs/tags/<ref>`
   (or reads a local tarball with `--from`);
2. Extracts it with the top-level folder stripped;
3. Checks that `SKILL.md` exists;
4. Prints the next steps (demo commands, docs).

It does not install R/Python packages. Use the repository's dependency checker after
installing:

```bash
python scripts/plan/04_check_deps.py --plan results/plan/plan.json --tier minimal
```

## Notes

- Node.js >= 18.
- No telemetry. The only network request is the GitHub download.
- `--from <file>` supports offline installs and CI.
- npm download counts are the reason this package exists: GitHub does not expose
  install counts for `npx skills add`.

## License

MIT — see [LICENSE](LICENSE). The installed repository contains third-party components
with their own licenses; see `third_party/NOTICE.md` after installation.
