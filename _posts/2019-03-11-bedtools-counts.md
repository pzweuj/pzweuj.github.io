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
bedtools genomecov -ibam xxx.sorted.bam -g reference.fa -bga > results.txt
```
例如：
得到下面这样的结果

chr1	100	101	300

chr1	101	102	311

chr1	102	103	400

chr1	103	104	320



这时候就可以用python脚本来统计了
```python
def getChromCoverage(inputfile, chrom):
	inputFile = open(inputfile, "r")
	length = 0
	covLength = 0
	counts = 0
	for line in inputFile:
		lineAS = line.split("\t")
		chromosome = lineAS[0]
		start = lineAS[1]
		end = lineAS[2]
		count = lineAS[3].split("\n")[0]
		if chromosome == chrom:
			length_tmp = int(end) - int(start)
			length += length_tmp
			counts += (int(count) * length_tmp)
			if count != "0":
				covLength_tmp = length_tmp
				covLength += covLength_tmp

	coverage = "%.2f" % ((float(covLength) / float(length)) * 100) + "%"
	depth_cov = "%.2f" % (float(counts) / float(covLength))
	depth_all = "%.2f" % (float(counts) / float(length))
	inputFile.close()

	return [coverage, depth_cov, depth_all]

```

这时，只要depth不等于0的，就统计，最后用n值除以统计区域的总长度，就是覆盖度了。
同理多定义几个统计的变量，就可以得到想要的深度的百分比了。比如说上面统计深度大于300的百分比。



[-_-]:LoveJing