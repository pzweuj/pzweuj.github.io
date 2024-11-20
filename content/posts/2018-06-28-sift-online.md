---
title: 使用在线的sift
tags: software
---

就是这个网址：
[sift](http://sift.bii.a-star.edu.sg/)。
SIFT根据序列同源性和氨基酸的物理性质预测氨基酸替换是否影响蛋白质功能。
应用于天然存在的非同义多态性和实验室诱导的错义突变。

![sift](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/sift-online.PNG)

这里提供了好几个在线的工具。这次我需要的是预测一个插入突变，所以选择下面的Classify coding indels。
然后选上参考基因组，看一下例子。
需要逗号分隔。
比如说下面这个：

X,99663401,99663401,1,GA

意思是，chr,start,end,orientation(1正链，-1负链)，indel(如果是ins就填插入的碱基，如果是del就填/)。

下面可以选择额外的信息，我是这样选的：

![sift](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/sift-select.PNG)

得到结果，可以看到预测是有害的。

![sift](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/sift-results.PNG)

然后点击Protein Sequence Change这一列可以看到氨基酸的改变情况哦。

![sift](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/sift-aachange.PNG)


[T_T]:不知道前景在干嘛。