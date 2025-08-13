---
title: 两分钟复刻世和的MINERVA评分计算器
tags: coding
---



世和的报告中搞了一个MINERVA评分，宣称是根据关键伴随基因的变异情况，知道早中期EGFR突变阳性非小细胞肺癌患者术后个性化辅助用药方案的选择。

看了一下出处是一篇合作文章，[PMID:34750392](https://pubmed.ncbi.nlm.nih.gov/34750392/)。前端见风就是雨，想着学人加这个。先不说Panel有没有覆盖人家模型所有的位点，也先不管世和有没有申请专利，首先看看计算方案能不能做出来。

本来打算是找文章看看的（后面看了，文章上也给了公式），但人家世和够大方，上线了一个在线的[计算器网页](https://minerva.geneseeq.com/)，我拉了源码看了下，直接就是静态的js，那么丢给AI即可。



## 两分钟重构

把html源码扒下来后，丢给[windsurf](https://windsurf.com/editor)（模型：claude 4），提示他使用python来重构这个计算器，两分钟后就拿到结果代码了。

代码请点击：[minerva_geneseeq.py](https://github.com/pzweuj/practice/blob/master/python/Minerva_Geneseeq/minerva_geneseeq.py)。



## 后记

感觉现在的AI编程软件的用户粘性确实不强。到目前为止，已经用过Cursor、Trae、Kiro、Augment、ClaudeCode、WindSurf、Cline等等，都是哪里有优惠就投奔到哪里。感觉上重要的还是内置的模型，无疑Claude是最好的。

欢迎使用[我的AFF链接](https://windsurf.com/refer?referral_code=7xn7c4yozdjv4fug)注册WindSurf。



