---
title: Gemini-CLI安装和使用
tags: software
---

众所周知，Gemini在各种意义上，是现在（2025年6月）最强的模型，[Gemini-CLI](https://github.com/google-gemini/gemini-cli)是google开源的以命令行的方式使用Gemini的软件。

## 安装

我准备在我的VPS里部署。下面的命令均是在`root`下操作，请注意。

首先需要安装nodeJs v18（注意，实际需要'>=20.18.1'）以上的版本。我的VPS是debian，由于debian的apt源里的node版本现在是v18.19.0，因此需要先更新源。

```bash
curl -fsSL https://deb.nodesource.com/setup_24.x -o nodesource_setup.sh
bash nodesource_setup.sh
apt-get install -y nodejs
```

确认版本，这样安装的node是v24版本了。

```bash
node -v
```

接下来安装Gemini-CLI。

```bash
npm install -g @google/gemini-cli
```

## 使用

认证方式分为下面三种：
- Login with Google
- Gemini API Key
- Vertex AI

我使用第二种Gemini API Key的模式。首先，需要有[Google AI Studio](https://aistudio.google.com/)的API Key。

可以将Key放到环境变量中，因为是自己的机器，所以直接放了

```bash
echo 'export GEMINI_API_KEY="YOUR_GEMINI_API_KEY"' >> ~/.bashrc
source ~/.bashrc
```

输入gemini就成功启动gemini-cli了，初次使用，需要设定主题，然后再进行认证。

```bash
gemini
```

连续两次Ctrl + C可以退出。

## FAILED_PRECONDITION

尽管进入了输入界面，但是我目前还是没有使用成功，回报的错误是

```
✕ [API Error: User location is not supported for the API use. (Status: FAILED_PRECONDITION)]
```

如上面所说，我的VPS是新加坡的，应该没有在Gemini的封锁区域里。查了一下[官方文档](https://ai.google.dev/gemini-api/docs/troubleshooting?hl=zh-cn)，说是需要在Google AI Studio中启用结算才可以。

