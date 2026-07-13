---
title: 为 GitHub 组织私有仓库配置 Deploy Key
tags: default
---

我有一个 GitHub 组织，组织下托管了一个私有的流程仓库。代码和文档都在上面，服务器需要能随时拉取最新版本，偶尔还要把自动生成的结果报告推送回去。问题来了：服务器没有"GitHub 账号"，它只是一台机器——给它配 SSH Key 的时候，这个 Key 挂在谁名下？

挂在个人账号下固然能跑，但有两个问题：一是权限太宽——我的个人 Key 能访问所有我有权限的仓库，放在服务器上不合适；二是如果哪天我离开实验室、账号权限被调整，服务器就断连了。正确的做法是给仓库配一把 **Deploy Key**——一把只属于这个仓库的独立钥匙。

## 三种方案对比

在 GitHub 生态里，让服务器访问私有仓库有三种主流方案：

| 方案 | 权限范围 | 依赖个人账号 | 适合场景 |
| :--: | :------: | :--------: | :------: |
| 个人 SSH Key | 账号下所有仓库 | 是 | 个人开发机 |
| Deploy Key | 单个仓库 | 否 | 服务器拉取/推送 |
| Machine User | 组织授权仓库 | 否（独立账号） | 多仓库 CI/CD |

对于"一台服务器 + 一个私有仓库"的场景，Deploy Key 是最简洁的选择：权限精准、不绑账号、配置简单。

> **Deploy Key 的局限**：一把 Deploy Key 只能绑定一个仓库。如果你有多个仓库需要同一台服务器访问，要么为每个仓库生成独立的 Key，要么考虑 Machine User 方案。但 Machine User 会消耗一个 GitHub 席位，对个人/小团队来说不太划算。

## 第一步：在服务器上生成 SSH Key

登录服务器，生成一把专用于这个仓库的密钥对：

```bash
# 生成 ED25519 密钥（推荐，比 RSA 更短更快）
ssh-keygen -t ed25519 -C "server-deploy-<repo-name>" -f ~/.ssh/id_ed25519_<repo-name>
```

`-C` 是注释，方便后续区分这把 Key 的用途。`-f` 指定密钥文件名——**不要覆盖默认的 `id_ed25519`**，那是你个人用的 Key。

生成后会得到两个文件：

| 文件 | 用途 | 是否上传 |
| :--: | :--: | :-----: |
| `id_ed25519_<repo-name>` | 私钥 | 留在服务器，绝不外传 |
| `id_ed25519_<repo-name>.pub` | 公钥 | 上传到 GitHub |

```bash
# 查看公钥内容，待会要用
cat ~/.ssh/id_ed25519_<repo-name>.pub
```

输出类似 `ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAI... server-deploy-<repo-name>`，复制整行。

> **注意**：生成密钥时不要设置 passphrase，否则服务器每次拉取都需要输入密码，自动化脚本就跑不起来了。如果担心私钥泄露，可以在服务器层面做文件权限控制（`chmod 600`），而不是依赖 passphrase。

## 第二步：在 GitHub 仓库添加 Deploy Key

这一步需要仓库管理员权限。组织管理员天然拥有组织中所有仓库的管理权限，所以你可以直接操作。

1. 打开 GitHub，进入目标仓库的 **Settings**（注意是仓库设置，不是组织设置）
2. 左侧菜单找到 **Deploy keys**（在 Security 分组下）
3. 点击 **Add deploy key**
4. 填写：
   - **Title**：给 Key 起个名字，比如 `Lab Server 01 - Analysis Pipeline`，方便日后识别
   - **Key**：粘贴上一步复制的公钥完整内容
   - **Allow write access**：勾选与否取决于需求

| 需求 | 是否勾选 Write Access |
| :--: | :------------------: |
| 只需 `git pull` / `git clone` | 不勾选 |
| 需要 `git push`（如推送自动生成的结果） | 勾选 ✔️ |

5. 点击 **Add key** 保存

添加成功后，这把 Key 会出现在仓库的 Deploy keys 列表中。如果以后服务器退役或密钥泄露，直接在这里删除即可，不影响其他仓库和账号。

> **组织仓库的 Deploy Key 管理**：作为组织管理员，你可以在组织所有仓库上添加 Deploy Key。但要注意——Deploy Key 是仓库级别的，不是组织级别的。如果你在组织下新建了仓库，需要单独为它配置 Deploy Key。

## 第三步：配置 SSH 使用指定密钥

服务器上可能已经有多个 SSH Key（个人的、其他项目的），需要让 Git 在访问这个仓库时使用正确的 Key。编辑 `~/.ssh/config`：

```bash
# 编辑 SSH 配置（不存在则新建）
vim ~/.ssh/config
```

添加以下内容：

```
# GitHub Deploy Key: <repo-name>
Host github.com-<repo-name>
    HostName github.com
    User git
    IdentityFile ~/.ssh/id_ed25519_<repo-name>
    IdentitiesOnly yes
```

`Host` 是一个别名，`IdentitiesOnly yes` 阻止 SSH 尝试其他密钥（避免触发 GitHub 的认证频率限制）。配置完成后，仓库的 remote URL 需要做对应修改：

```bash
# 克隆新仓库时，把 github.com 替换为别名
git clone git@github.com-<repo-name>:<org>/<repo>.git

# 已有仓库修改 remote
git remote set-url origin git@github.com-<repo-name>:<org>/<repo>.git
```

关键是 `git@github.com` → `git@github.com-<repo-name>`，让 SSH 根据 Host 别名匹配到对应的密钥文件。

## 第四步：验证连接

```bash
# 测试 SSH 认证
ssh -T git@github.com-<repo-name>
```

如果配置正确，会看到：

```
Hi <org>/<repo>! You've successfully authenticated, but GitHub does not provide shell access.
```

注意这里的输出是 `<org>/<repo>` 而不是你的用户名——这说明 GitHub 识别出这是一个 Deploy Key 认证，而非个人账号认证。权限范围被限制在这个仓库内。

接下来试一次拉取：

```bash
git fetch origin
```

能正常拉取就算配置成功。

## 多仓库场景的扩展

如果后续需要同一台服务器访问组织的多个私有仓库，有两种做法：

**方案 A：每仓库一把 Key，SSH Config 加多个 Host**

```bash
# 为每个仓库生成独立 Key
ssh-keygen -t ed25519 -C "deploy-repo-a" -f ~/.ssh/id_ed25519_repo-a
ssh-keygen -t ed25519 -C "deploy-repo-b" -f ~/.ssh/id_ed25519_repo-b

# SSH config 中分别为每个 Host 配置
Host github.com-repo-a
    HostName github.com
    User git
    IdentityFile ~/.ssh/id_ed25519_repo-a
    IdentitiesOnly yes

Host github.com-repo-b
    HostName github.com
    User git
    IdentityFile ~/.ssh/id_ed25519_repo-b
    IdentitiesOnly yes
```

优点是一把 Key 一个仓库，互不干扰，撤销某台机器的某个仓库权限只需删除对应的 Key。缺点是 Key 数量随仓库数增长。

**方案 B：共享 Key + 多仓库注册**

同一把公钥可以添加到多个仓库的 Deploy Key 中。在 GitHub 上给每个目标仓库都添加同一把公钥即可。优点是一把 Key 走天下，缺点是撤销权限时需要逐个仓库删除。

> 我个人倾向于方案 A——Key 文件多几个不碍事，但权限边界清晰，哪天出问题排查起来不费劲。

## 安全注意事项

- **私钥文件权限**：`chmod 600 ~/.ssh/id_ed25519_*`，禁止其他用户读取
- **不要复用个人 Key**：服务器上始终使用独立的 Deploy Key，不要用个人的 SSH Key
- **定期轮换**：如果服务器的 Key 曾有外传记录，或团队成员离职，及时删除旧 Key 并生成新 Key
- **Write Access 审慎开启**：只在确实需要推送时才勾选写权限。对于只读拉取（如部署脚本），一律不勾选。这样即使服务器被攻破，攻击者也推不了代码

## 结论

给组织私有仓库配置 Deploy Key 是"服务器拉代码"场景的标准解法。核心思路就三步：

1. **服务器上生成独立密钥对**——不绑定个人账号
2. **在仓库设置中添加 Deploy Key**——权限只作用于这个仓库
3. **SSH Config 做 Host 别名匹配**——让 Git 自动使用正确的 Key

整个过程五分钟搞定，换来的是权限隔离、账号无关、可随时撤销的安全保障。如果你的服务器还在用个人 SSH Key 拉私有仓库，今天就可以换掉。