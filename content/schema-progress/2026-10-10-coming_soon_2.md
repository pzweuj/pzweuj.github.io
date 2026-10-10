---
title: 其实道爷我已经成了
tags: default
---

## 当前进度

项目开发做完了，整条链路已经走通。现在剩下的是前端，还在收细节。

YiJian、Octopus、Sepiida、schema-germline 这四个仓库一直是 Public 的。正式发布我按 12 月左右来排。

下面是现在的界面。前端还在调，不是最终状态。右上角的积分和任务费用不用看，开源版本没有这块。

另外，ACMG 规则 2027 年要改。AI 评估会在这个阶段做完。

## 总览

打开一个已完成的任务，先落到总览。样本跑完之后，检出有多少、哪些结果加载了、质控怎么样、性别对不对得上，都在这一页。

![pic1](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/yijian_20261010_1.png)

上面是四个汇总：已加载检出、已加载回报、分析成员、参考基因组。左边按结果类型列数量和数据状态，包括 SNP/InDel、CNV 区段、CNV 外显子、STR、MEI、线粒体、UPD、ROH。右边是质控，从过滤前后的 reads、比对率、捕获效率，到深度、覆盖率、重复率和 Q30。最下面核对登记性别和数据里的性别。

## SNP / InDel

变异判读按结果类型分栏。SNP/InDel 这一栏是一张表，可以筛选，也可以导出。

![pic2](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/yijian_20261010_2.png)

每一行是一个变异。列上是基因、位置、参考和变异、杂合性、VAF、深度、gnomAD 的总体和东亚频率、ClinVar，再往右是 ACMG 评定、转录本和 cDNA 变化。

## 展开

点开一行，右边拉出这条变异的详情。

![pic3](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/yijian_20261010_3.png)

这里放的是证据，不是一个直接拍板的结论。基因和疾病的关联、遗传模式、疾病标识、dbSNP、HGNC 都在这一侧。证据不够、形不成 ACMG 分类的，就标成证据不足。
