---
title: 用bedtools统计深度和覆盖度
tags: coding
key: bedtools-genomecov
---

[bedtools](https://bedtools.readthedocs.io/en/latest/) 也是一个常驻我环境变量的软件了。

突然有一个需求，是统计任意深度的百分比，一时半会，没想到有什么现成的工具可以用，所以就用bedtools来统计了。
操作很简单，需要的输入文件是一个排过序的bam。

看一下统计原理:
![bedtools-genomecov](https://bedtools.readthedocs.io/en/latest/_images/genomecov-glyph.png)

可以看出，当加入参数-d -split的时候，可以得到比较好进行下一步统计的结果。

```bash
bedtools genomecov -ibam xxx.sorted.bam -g reference.fa -d -split > results.txt
```
例如：
得到下面这样的结果

chr1	100	101	300
chr1	101	102	311
chr1	102	103	400
chr1	103	104	320

这时候就可以用python脚本来统计了
```python
# 读入文件
inputdata = open("results.txt", "r")
n = 0
m = 0
for line in inputdata:
	lineAfterSplit = line.split("\t")
	chrom = lineAfterSplit[0]
	start = lineAfterSplit[1]
	end = lineAfterSplit[2]
	depth = int(lineAfterSplit[3])
	
	# 覆盖度统计
	if depth != 0:
		n += 1
	
    # >300x的覆盖度
	if depth >= 300:
		m += 1
	
	# 平均深度
	sum_depth += depth
		
# 最后除以区域总长，这个自己想方法算
percent = float(n / (103 - 100 + 1))
# 平均深度
avg_depth = float(sum_depth / (103 - 100 + 1))

```

这时，只要depth不等于0的，就统计，最后用n值除以统计区域的总长度，就是覆盖度了。
同理多定义几个统计的变量，就可以得到想要的深度的百分比了。比如说上面统计深度大于300的百分比。



[-_-]:LoveJing