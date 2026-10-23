# 发布 npm 安装器（逐步操作）

> 目标：把 `npm/` 里的安装器发布为 `rice-skills`，从而获得 npm 的下载量统计。
> 当前状态：包名 `rice-skills` 可用；本机未登录 npm（`npm whoami` → ENEEDAUTH）。

## 0. 前置检查

```bash
node -v          # >= 18
npm -v
npm view rice-skills version   # E404 表示名字可用；若返回版本号说明已被占用
```

## 1. 注册 npm 账号

1. 打开 https://www.npmjs.com/signup
2. 填 username / email / password
3. 到邮箱点验证链接
4. 建议立即开启 2FA（Account → Two-Factor Authentication）

## 2. 登录

```bash
npm login
```

按提示输入 username、password、email；若开了 2FA，会要求一次性验证码（OTP）。

验证：

```bash
npm whoami        # 应输出你的用户名
```

## 3. 发布

```bash
cd /path/to/Rice-skills/npm
npm publish --access public
```

成功后会输出：

```
+ rice-skills@3.3.9
```

> 若提示 `E402 You do not have permission to publish "rice-skills"`：名字已被抢注。
> 改 `package.json` 的 `"name": "@shuke-bio/rice-skills"`，再执行
> `npm publish --access public`。代码无需改动。

## 4. 验证

```bash
npm view rice-skills version          # 应返回 3.3.9
npx rice-skills@3.3.9 check --dir /tmp/rs-check
```

在空目录完整安装一次：

```bash
mkdir -p /tmp/rs-test && cd /tmp/rs-test
npx rice-skills install --dir ./Rice-skills
ls Rice-skills/SKILL.md
```

## 5. 自动化：Release 时自动发布

1. npmjs.com → Access Tokens → Generate New Token
   - 类型选 **Automation**（绕过 2FA，适合 CI）；
   - 复制 token（只显示一次）。
2. GitHub 仓库 → Settings → Secrets and variables → Actions → New repository secret
   - Name：`NPM_TOKEN`
   - Value：粘贴 token
3. 之后每次发布 GitHub Release，`.github/workflows/npm-publish.yml` 会自动执行
   `npm ci && npm test && npm publish`。

> 注意：`npm/package.json` 的 `version` 必须与 Release tag 一致（如 tag `v3.4.0` →
> package version `3.4.0`），否则 npm 会拒绝重复版本。

## 6. 每次发版的版本同步

```bash
# 仓库根版本更新后
cd npm
npm version 3.4.0 --no-git-tag-version   # 同步 package.json / lock
npm publish --access public
```

## 7. 常见错误

| 错误 | 原因 | 处理 |
|---|---|---|
| `ENEEDAUTH` | 未登录 | `npm login` |
| `E402` | 包名被占用 | 改用 `@shuke-bio/rice-skills` |
| `E403` | 2FA/权限 | 用 Automation token，或输入 OTP |
| `EPUBLISHCONFLICT` | 版本已存在 | 提升版本号后重发 |
| `npx` 仍装旧版 | npx 缓存 | `npx rice-skills@latest ...` 或清缓存 |

## 8. 发布后要做的两件事

1. README 顶部加下载量徽章：
   ```markdown
   [![npm](https://img.shields.io/npm/dw/rice-skills)](https://www.npmjs.com/package/rice-skills)
   ```
2. 在安装章节补一条：
   ```bash
   npx rice-skills install --dir ~/.workbuddy/skills/Rice-skills
   ```
