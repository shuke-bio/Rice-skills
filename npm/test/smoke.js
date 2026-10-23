#!/usr/bin/env node
/** Smoke test: help/version + offline install from a local tarball. */
const { spawnSync } = require("child_process");
const fs = require("fs");
const os = require("os");
const path = require("path");

const bin = path.join(__dirname, "..", "bin", "rice-skills.js");
const run = (...args) => spawnSync(process.execPath, [bin, ...args], { encoding: "utf8" });
const runIn = (cwd, ...args) => spawnSync(process.execPath, [bin, ...args], { encoding: "utf8", cwd });

let r = run("version");
if (r.status !== 0 || !/^\d+\.\d+\.\d+/.test(r.stdout.trim())) throw new Error("version failed");
console.log("version OK:", r.stdout.trim());

r = run("help");
if (r.status !== 0 || !r.stdout.includes("rice-skills install")) throw new Error("help failed");
console.log("help OK");

const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "rs-test-"));
const out = path.join(tmp, "pkg");
fs.mkdirSync(path.join(out, "SKILL.md", ".."), { recursive: true });
fs.writeFileSync(path.join(out, "SKILL.md"), "# test\n");
const archive = path.join(tmp, "fixture.tar.gz");
const tarLib = require("tar");
tarLib.c({ gzip: true, file: archive, cwd: path.join(out, "..") }, ["pkg"]).then(() => {
const dest = path.join(tmp, "install");
r = run("install", "--from", archive, "--dir", dest);
if (r.status !== 0) throw new Error("install failed: " + r.stderr);
if (!fs.existsSync(path.join(dest, "SKILL.md"))) throw new Error("SKILL.md missing after install");
console.log("offline install OK");

const proj = fs.mkdtempSync(path.join(os.tmpdir(), "rs-proj-"));
r = runIn(proj, "install", "--project", "--from", archive, "--write-agents");
if (r.status !== 0) throw new Error("project install failed: " + r.stderr);
if (!fs.existsSync(path.join(proj, ".workbuddy", "skills", "Rice-skills", "SKILL.md")))
  throw new Error("project skill not installed");
if (!fs.existsSync(path.join(proj, "AGENTS.md"))) throw new Error("AGENTS.md not created");
console.log("project install + AGENTS OK");

fs.rmSync(tmp, { recursive: true, force: true });
fs.rmSync(proj, { recursive: true, force: true });
}).catch((e) => { throw e; });
