---
title: KataGo和Sabaki安装（Windows版本）
tags: default
key: katago
---

## 背景
[KataGo](https://github.com/lightvector/KataGo)是一款采取和[AlphaGo Zero](https://deepmind.com/blog/article/alphago-zero-starting-scratch)相同策略的围棋AI程序，号称平民狗。而[Sabaki](https://github.com/SabakiHQ/Sabaki)是一个围棋AI引擎的GUI，两者搭配，就相当于一个围棋AI的GUI软件。

## 准备工作
首先请先确定自己能用哪个版本的KataGo，是OpenCL还是CUDA还是Eigen，三个版本的对比在[这里](https://github.com/lightvector/KataGo#opencl-vs-cuda-vs-eigen)。下载自己需要的版本，然后再下载[model文件](https://github.com/lightvector/KataGo/releases/tag/v1.4.5)。KataGo提供了2020年6月的训练结果，建议上说GPU越强，就用blocks数字越大的。最后去下载sabaki程序。

分别是以下三个。
```
katago-v1.7.0-gpu-opencl-windows-x64.zip
g170-b40c256x2-s5095420928-d1229425124.bin.gz
sabaki-v0.51.1-win-x64-portable.exe
```

即下载下来以上三个文件。


## KataGo配置
解压katago-v1.7.0-gpu-opencl-windows-x64.zip。适用cmd，到达kataGo程序文件夹下，使用以下命令进行配置，我将g170-b40c256x2-s5095420928-d1229425124.bin.gz文件放在kataGo下的model文件夹中了。

```
katago.exe genconfig -model model\g170-b40c256x2-s5095420928-d1229425124.bin.gz -output gtp_g170_b40.cfg
```

然后需要回答几个问题，分别是用什么规则，这里选择chinese，即中国规则。然后是时间规则，n就好了。之后是是否允许KataGo在对手读秒时进行计算，默认是否，那么我们也选n就好了。然后是选择使用运算的设备，一般是选分数最大的，因为我只有一块GPU，所以直接选0了。之后是设置内存，KataGo最大会用3GB，那直接给好了。后面的问题是和GPU相关的，都默认直接回车确定就好了。


![1](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/katago-1.jpg)


## Sabaki配置
语言可以在File-Preferences下修改为中文。然后选择引擎-显示引擎侧边栏。点击左边播放符号-管理引擎-新增，按下图方式配置好。

![2](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/katago-2.jpg)

```
gtp -model model\g170-b40c256x2-s5095420928-d1229425124.bin.gz -config gtp_g170_b40.cfg
```

点击新增完成保存。

在引擎处选择新增的kataGo引擎，右键选择设为分析器，就可以使用了。
![3](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/katago-3.jpg)


## 其他
还有一个更傻瓜式的版本，[KataTrain](https://github.com/sanderland/katrain)，已经内置了KataGo，直接下载下来就能用！
![katatrain](https://raw.githubusercontent.com/sanderland/katrain/master/screenshots/analysis.png)



[-_-]:xx