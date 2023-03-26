---
title: Stable Diffusion
tags: software
---

2022年下半年开始，画图AI层出不穷，如DALL-E、NovelAI等，然后到了22年年底，ChatGPT横空出世，再带火了一次聊天AI。目前ChatGPT只有在线版本，不能离线部署，虽然没吃上热饭，至少喝一口汤。

这里是本地部署[Stable Diffusion](https://github.com/CompVis/stable-diffusion)来画图，没有艺术细胞的人也能画出漂亮小姐姐啦。另外也有大佬给Stable Diffusion开发了UI（[stable-diffusion-webui](https://github.com/AUTOMATIC1111/stable-diffusion-webui)），对新手更友好了。

我的破电脑是2070，有8G显存，应该够用的。



## 安装

首先需要电脑环境中有python 3.10.6（我的python是3.11.2，试了一下，装pytorch时会报错，因此老老实实用conda创建一个python 3.10.6的环境了），然后自行安装pytorch。pytorch下载要2.3G，还好我系统盘是个2T的M2 SSD。

```cmd
conda create -n stw python=3.10.6
conda activate stw
python -m pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117 --no-cache-dir
conda install git
conda install numpy matplotlib
pip install ftfy regex tqdm
```

这样就安装好依赖了。


然后下载stable-diffusion-webui并解压
```
https://hub.gitmirror.com/https://github.com/AUTOMATIC1111/stable-diffusion-webui/archive/refs/tags/v1.0.0-pre.zip
```

切换到目录并运行webui-user.bat，另外由于我使用了conda环境，修改webui-user.bat中的python路径为环境中python，找路径方法如下（python脚本）。还是遇到了不少坑的，比如国内的某堵高墙，git clone超时，我找了一些镜像仓库来clone后成功了。还有比如numpy的版本需要1.24以下等等。

```python
import sys
print(sys.path)
```

然后是我具体的webui-user.bat设置
```cmd
@echo off

set PYTHON=C:\Users\pzw-pc\.conda\envs\stw\python.exe
set GIT=C:\Users\pzw-pc\.conda\envs\stw\Library\bin\git.exe
set VENV_DIR=C:\Users\pzw-pc\.conda\envs\stw
set COMMANDLINE_ARGS= --xformers

call webui.bat
```



## 模型文件

下载模型文件，并放到models/Stable-diffusion下

```
https://huggingface.co/CompVis/stable-diffusion-v-1-4-original/resolve/main/sd-v1-4.ckpt
```

下面这个是[CivitAI](https://civitai.com/)图像社区的共享模型中心

```
https://civitai.com/models/6424/chilloutmix
```


## LORA

模型风格微调文件，生成的图像会接近对应的LORA文件。放到models/Lora下。


## 修复面部

修复面部会用到[CodeFormer](https://github.com/sczhou/CodeFormer)项目，具体在webui部署好后，可以自动下载。

