---
title: FusionDraw 基因融合示意图绘制工具
tags: coding
---

经常收到帮忙绘制融合图的需求，之前我都是用Arriba画的。因此，SchemaBio也有了[draw_fusion](https://github.com/SchemaBio/draw_fusion)这个项目。

问题是，我希望这个事情更轻量化。draw_fusion现在带有了一套笨重的R环境，同时需要GTF来获取信息，而不在GTF上的点，就画不出来。

我的想法是，我们不过是要画一幅图，我自己觉得合理即可，管它实际有没有这个基因。因此诞生了[FusionDraw](https://github.com/pzweuj/FusionDraw)。

### 你自己可以用

![fusiondraw_1](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/fusiondraw_1.png)

无需多说，直接线上使用[FusionDraw](https://fusiondraw.biotools.space)。注意，这个项目所有的绘图信息，需要用户自己人工确认，因为后台没有内置GTF（只内置了cytoband），也就是绘图信息都是依赖自己去核对的。

当然，这产生的缺陷是**无法按照外显子的实际长度对渲染的矩形进行scale**，要是你的需求是精细化的话，那我还是建议你回归到[draw_fusion](https://github.com/SchemaBio/draw_fusion)。

> 当然，FusionDraw的结果是一个SVG图，你可以让Agent去改的



### 你家的Agent也可以用

让用户自己拼出完整的 PlotSpec，仍然要求他们理解不少内部字段。因此项目还提供了 [fusiondraw-svg](https://github.com/pzweuj/FusionDraw/blob/main/skills/fusiondraw-svg/SKILL.md) SKILL，告诉 Agent 应该如何收集信息、何时追问，以及什么情况下绝对不能猜。

一个自然语言请求可以从“画一张 EML4-ALK 融合图”开始。Agent 接下来应该确认保留的外显子段；如果用户要求真实基因组视图，再收集顶层 assembly、染色体、断点坐标和链方向；如果 resolution 不完整，就先停下来补齐 `region`、`codingRegion` 和相应的 exon/intron 编号。

> 绘图信息都是依赖自己去核对

这一点，显然我们懒人也是不用的，也全权交给你的大模型即可。