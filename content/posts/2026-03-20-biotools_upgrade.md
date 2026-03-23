---
title: BioTools超进化
tags: software
---


[BioTools](https://github.com/pzweuj/biotools)经过持续的迭代，现在拥有了更多的工具。

一开始，BioTools的初衷是所有工具都是在本地进行运算的，所有信息都不会离开你自己的电脑。但是，实际情况是，我们增加了下面这些工具，会向外部服务器发送信息。当然，在这些工具页面里，我们都著名了请求将被发送到外部服务器。

BioTools仍然提供着[在线使用](https://use.biotools.space/)的方式。

### 向外请求

首先是[Mutalyzer](https://use.biotools.space/tools/mutalyzer)，这是一个将HGVS变异描述验证、标准化、格式转换和序列映射的工具。请求会被发送到Mutalyzer官方，然后根据返回结果，自行渲染。

然后是[剪接位点有害性预测工具](https://use.biotools.space/tools/spliceai)，这个工具可以使用SpliceAI和Pangolin模型预测变异对剪接位点的影响，支持GRCh37和GRCh38。请求将被发送到博德研究所的SpliceAI Lookup。

还有是[Transvar](https://use.biotools.space/tools/transvar)，这个工具可以反向注释pHGVS和cHGVS等，支持GRCh37和GRCh38。请求将被发送我自己在抱抱脸部署的space。

### 跳转

当然，现在BioTools还包含了下面三个外部工具的跳转。

1. [ManeLoca](https://use.biotools.space/tools/maneloca)，这会跳转到我编写的MANE转录本坐标查询器
2. [DeepHPO](https://use.biotools.space/tools/deephpo)，这会跳转到我编写的基于LLM的病历信息HPO术语提取
3. [Warfarin](https://use.biotools.space/tools/warfarin)，这会跳转到我编写的华法林剂量计算器

### 实验性更新

BioTools现在加入了Page Agent，可以以自然语言的方式，进行完整的页面操控。
