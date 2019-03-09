---
title: 用R处理ab1数据
tags: coding
---
不多说，直接上代码！

```r
#加载sangerseqR包
library(sangerseqR)

#读入数据
seq = readsangerseq('input.ab1')

#读取碱基数据，0.33指的是将达到主峰0.33的次峰定义为杂合子峰
bc = makeBaseCalls(seq, ratio = 0.33)

#读主峰
primarySeq(seq)

#读次峰
secondarySeq(seq)

#输出可视化图像
chromatogram(bc, showcalls = 'both')
```

结果如下：
![sangerseqR-results](https://github.com/pzweuj/pzweuj.github.io/raw/master/downloads/images/sangerseqR.png)


[^_^]:我最近学的东西太少，快不够用了
