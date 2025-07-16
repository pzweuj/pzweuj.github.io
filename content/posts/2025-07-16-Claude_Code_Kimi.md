---
title: 更新Claude Code到工具库并接入Kimi K2
tags: default
---

此前在做前端项目时，一直在白嫖Cursor和Augment。一直听说Claude Code好用，但是Claude 4的价格非常高，同时公司又不愿意采购，因此一直没有采用。

Kimi的K2出来后，可以接入到Claude Code里，同时价格是Claude的五分之一，编程能力据说媲美Gemini 2.5Pro，可以尝试一下。

## 对比其他

根据网上信息，对比Cursor、Augment、Cline、Roo Code等，Claude Code的优势是结果更精准，更能一步达到目的。缺点是只有命令行的交互。

使用的方式是需先制定计划，让Claude Code自行读取计划并执行任务。


## 前置准备

安装Claude Code。

```bash
npm install -g @anthropic-ai/claude-code
```

到[月之暗面申请一个API Key](https://platform.moonshot.cn/console/api-keys)。新用户默认赠送15元人民币。

在Windows下还需要[安装Git](https://git-scm.com/downloads/win)。


## 配置环境变量

配置环境变量即可使用。

在Linux下：

```bash
export ANTHROPIC_AUTH_TOKEN=sk-xxxxxxxxxxxxxxxxxxxxx
export ANTHROPIC_BASE_URL=https://api.moonshot.cn/anthropic
```

在Windows下（PowerShell）：
```powershell
$env:ANTHROPIC_AUTH_TOKEN="sk-xxxxxxxxxxxxxxxxxxxxx"
$env:ANTHROPIC_BASE_URL="https://api.moonshot.cn/anthropic"
```

如果需要永久固化，在Linux可修改 `.bashrc`。

在Windows下，注意设置后，在**新的终端**中才会生效：
```powershell
setx ANTHROPIC_AUTH_TOKEN sk-xxxxxxxxxxxxxxxxxxxxx
setx ANTHROPIC_BASE_URL https://api.moonshot.cn/anthropic
```

## 使用

在终端输入 claude 来使用 Claude Code。

```bash
claude
```

## 其他

在未认证前，Kimi是无法充值的。在使用时还会被限流，并返回429错误。如果[需要提高速率](https://platform.moonshot.cn/docs/pricing/limits#%E9%99%90%E9%80%9F%E6%A6%82%E5%BF%B5%E8%A7%A3%E9%87%8A)，需要进行认证，并充值最小的金额50元（但仍然未能达到生产力的速率）。
