---
title: Claude Code接入智谱GLM4.5
tags: default
---

实测将智谱GLM4.5接入Claude Code比Kimi K2和Qwen3都要好。


## 前置准备

安装Claude Code。

```bash
npm install -g @anthropic-ai/claude-code
```

申请一个智谱BigModel平台的账户，可以用我的[AFF链接](https://www.bigmodel.cn/invite?icode=StKTbhPWugGzeS%2BG8KgW%2Bf2gad6AKpjZefIo3dVEQyA%3D)，这样我们都可以获得GLM-4.5-Air的2000万tokens。

然后在[API Key](https://bigmodel.cn/usercenter/proj-mgmt/apikeys)中创建一个新的Key用于接入Claude Code。

在Windows下还需要[安装Git](https://git-scm.com/downloads/win)。


## 配置环境变量

配置环境变量即可使用。

在Linux下：

```bash
export ANTHROPIC_AUTH_TOKEN=xxxxxxxxxxxxxxxxxxxxx
export ANTHROPIC_BASE_URL=https://open.bigmodel.cn/api/anthropic
```

在Windows下（PowerShell）：
```powershell
$env:ANTHROPIC_AUTH_TOKEN="xxxxxxxxxxxxxxxxxxxxx"
$env:ANTHROPIC_BASE_URL="https://open.bigmodel.cn/api/anthropic"
```

如果需要永久固化，在Linux可修改 `.bashrc`。

在Windows下，注意设置后，在**新的终端**中才会生效：
```powershell
setx ANTHROPIC_AUTH_TOKEN xxxxxxxxxxxxxxxxxxxxx
setx ANTHROPIC_BASE_URL https://open.bigmodel.cn/api/anthropic
```

## 使用

在终端输入 claude 来使用 Claude Code。

```bash
claude
```

## 注意

BigModel平台新用户赠送的GLM4.5额度是200万tokens，同时会赠送200万的不限模型tokens，即我们可以白嫖的tokens数只有400万，只能做一个小项目。同时V0等级用户的GLM-4.5并发限制是20。当前GLM-4.5的定价是70元/1000万tokens，限期3个月（当前在打1折，6.9元出售）。

## 氛围编程

我现在的项目建立流程基本是这样的：

1. 如果是前端，使用V0构建基础的框架；
2. 用Kiro实现基本的功能，一直薅到每天的用量上限；
3. Kiro达到上限后，切换到Claude Code + GLM继续。
