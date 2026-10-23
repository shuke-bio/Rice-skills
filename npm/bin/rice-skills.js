#!/usr/bin/env node
/**
 * rice-skills — installer for the Rice-skills repository.
 *
 * The npm package contains only this installer; analysis scripts are downloaded
 * from the pinned GitHub tag so the tarball stays small and versions stay explicit.
 */
"use strict";

const fs = require("fs");
const os = require("os");
const path = require("path");
const https = require("https");
const { spawnSync } = require("child_process");
const tar = require("tar");

const PKG = require("../package.json");
const REPO = "shuke-bio/Rice-skills";
const DEFAULT_REF = `v${PKG.version}`;

const AGENTS_TEXT = `# Project conventions

- Rice omics analyses in this workspace use \`.workbuddy/skills/Rice-skills\`.
- Before starting, read \`SKILL.md\` and \`docs/decision/README.md\` from that folder.
- Always present an analysis + figure plan first, and wait for user approval before running.
- Check dependencies first and list anything missing; install only after explicit consent.
- Demo data are program-generated and must not be used for biological conclusions.
`;

function writeAgents() {
  const file = path.join(process.cwd(), "AGENTS.md");
  if (fs.existsSync(file)) {
    console.log("AGENTS.md already exists — left unchanged.");
    return;
  }
  fs.writeFileSync(file, AGENTS_TEXT, "utf8");
  console.log("Created AGENTS.md (project-level skill binding).");
}

function help() {
  console.log(`rice-skills ${PKG.version}

Usage:
  rice-skills install [--dir <path> | --project] [--ref <tag>] [--force] [--from <tarball>] [--write-agents]
  rice-skills check   [--dir <path> | --project]
  rice-skills version
  rice-skills help

Options:
  --dir <path>    Target directory (default: ./Rice-skills)
  --project       Install into ./.workbuddy/skills/Rice-skills (project-level skill)
  --ref <tag>     GitHub tag to install (default: v${PKG.version})
  --force         Install into a non-empty directory
  --from <file>   Install from a local .tar.gz instead of downloading
  --write-agents  Create ./AGENTS.md binding this project to the installed skill
                  (never overwrites an existing AGENTS.md)

Examples:
  npx rice-skills install
  npx rice-skills install --dir ~/.workbuddy/skills/Rice-skills
  npx rice-skills check --dir ./Rice-skills
`);
}

function parse(argv) {
  const out = { _: [] };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a.startsWith("--")) {
      const key = a.slice(2);
      const next = argv[i + 1];
      if (next && !next.startsWith("--")) { out[key] = next; i++; }
      else out[key] = true;
    } else out._.push(a);
  }
  return out;
}

function download(url, dest) {
  return new Promise((resolve, reject) => {
    const file = fs.createWriteStream(dest);
    const get = (u, depth = 0) => {
      if (depth > 5) return reject(new Error("Too many redirects"));
      https.get(u, { headers: { "User-Agent": `rice-skills/${PKG.version}` } }, (res) => {
        if (res.statusCode >= 300 && res.statusCode < 400 && res.headers.location) {
          res.resume();
          return get(res.headers.location, depth + 1);
        }
        if (res.statusCode !== 200) {
          res.resume();
          return reject(new Error(`Download failed: HTTP ${res.statusCode} (${u})`));
        }
        const total = Number(res.headers["content-length"] || 0);
        let got = 0, last = 0;
        res.on("data", (c) => {
          got += c.length;
          if (total && Date.now() - last > 500) {
            last = Date.now();
            process.stdout.write(`\r  downloading ${(got / 1048576).toFixed(1)} / ${(total / 1048576).toFixed(1)} MB`);
          }
        });
        res.pipe(file);
        file.on("finish", () => { file.close(); if (total) process.stdout.write("\n"); resolve(); });
      }).on("error", reject);
    };
    get(url);
  });
}

async function install(opts) {
  if (opts.project && opts.dir) {
    console.error("Use either --project or --dir, not both.");
    process.exit(1);
  }
  const target = opts.project
    ? path.join(process.cwd(), ".workbuddy", "skills", "Rice-skills")
    : path.resolve(opts.dir || "Rice-skills");
  const ref = String(opts.ref || DEFAULT_REF);
  const force = Boolean(opts.force);

  if (fs.existsSync(target) && fs.readdirSync(target).length && !force) {
    console.error(`Target exists and is not empty: ${target}\nUse --force to overwrite, or --dir to choose another path.`);
    process.exit(1);
  }
  fs.mkdirSync(target, { recursive: true });

  let archive = opts.from;
  let temp = false;
  if (!archive) {
    archive = path.join(os.tmpdir(), `rice-skills-${ref}-${Date.now()}.tar.gz`);
    temp = true;
    const url = `https://codeload.github.com/${REPO}/tar.gz/refs/tags/${ref}`;
    console.log(`Downloading ${REPO} ${ref} ...`);
    await download(url, archive);
  } else if (!fs.existsSync(archive)) {
    console.error(`Archive not found: ${archive}`);
    process.exit(1);
  }

  console.log(`Extracting to ${target} ...`);
  await tar.x({ file: archive, cwd: target, strip: 1 });
  if (temp) fs.rmSync(archive, { force: true });

  const skill = path.join(target, "SKILL.md");
  if (!fs.existsSync(skill)) {
    console.error("Extraction finished but SKILL.md is missing — the archive layout may have changed.");
    process.exit(1);
  }
  console.log(`\nInstalled Rice-skills ${ref} -> ${target}`);
  console.log("Next: see docs/installation.md, then run a demo:");
  console.log("  Rscript scripts/gwas/09_multilocus_qtn.R --demo");
  console.log("  python scripts/design/03_qpcr_analysis.py --demo");
  if (opts["write-agents"]) writeAgents();
}

function check(opts) {
  const dir = opts.project
    ? path.join(process.cwd(), ".workbuddy", "skills", "Rice-skills")
    : path.resolve(opts.dir || "Rice-skills");
  const required = ["SKILL.md", "README.md", "config/nature_palette.json", "scripts"];
  let ok = true;
  for (const f of required) {
    const exists = fs.existsSync(path.join(dir, f));
    console.log(`${exists ? "OK  " : "MISS"} ${f}`);
    if (!exists) ok = false;
  }
  const scriptsDir = path.join(dir, "scripts");
  if (fs.existsSync(scriptsDir)) {
    let n = 0;
    const walk = (d) => {
      for (const e of fs.readdirSync(d, { withFileTypes: true })) {
        if (e.name === "__pycache__") continue;
        const p = path.join(d, e.name);
        if (e.isDirectory()) walk(p);
        else if (/\.(R|py|sh)$/.test(e.name)) n++;
      }
    };
    walk(scriptsDir);
    console.log(`scripts: ${n}`);
  }
  const r = spawnSync("Rscript", ["--version"], { encoding: "utf8" });
  console.log(`Rscript: ${r.status === 0 ? (r.stdout || r.stderr).split("\n")[0].trim() : "not found (optional for Python-only use)"}`);
  const py = spawnSync("python", ["--version"], { encoding: "utf8" });
  console.log(`python : ${py.status === 0 ? (py.stdout || py.stderr).trim() : "not found"}`);
  process.exit(ok ? 0 : 1);
}

async function main() {
  const opts = parse(process.argv.slice(2));
  const cmd = opts._[0] || "help";
  if (cmd === "install") return install(opts);
  if (cmd === "check") return check(opts);
  if (cmd === "version" || opts.version) { console.log(PKG.version); return; }
  help();
}

main().catch((e) => { console.error(e.message || e); process.exit(1); });
