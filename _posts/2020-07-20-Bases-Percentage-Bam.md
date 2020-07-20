---
title: 统计Bam文件单个位点碱基比例
tags: coding
---

现有一个bam文件，需要统计一个位点的碱基比例。

在查询相关资料时，[有人](https://www.jianshu.com/p/0697dbfe316c)是用[bam-readcount](https://github.com/genome/bam-readcount)这个软件做的。使用现成的软件当然好，但是我们的需求是要帮他人升级这个功能，所以还是使用现有的软件而不是引入新软件比较好，因此这里使用samtools，基本上大家都会有。

使用samtools的mpileup功能。比方说我要chr1:1111111这个位点，那么我就只输出这个位置的pileup格式文件就好了。

```bash
samtools mpileup -f ref.fa -r chr1:1111111-1111111 in.bam > out.pileup
```

然后结果格式是这样的
```
chr1	1111111    A       11      ,,,..T,....     BHIGDGIJ?FF
```

在这里，第一列是contig名称，第二列是位置，第三列是参考碱基，第四列是reads数，第五列是结果，第六列是质量值。

我们可以通过解析第五列获得要的结果。在samtools的文档中提到，“.”表示与参考序列正链匹配，"，"表示与参考序列负链匹配，因此用python这样写了

```python
pileup = open("out.pileup", "r")

for line in pileup:
	linesplit = line.split("\t")
	chrom = linesplit[0]
	position = linesplit[1]
	ref = linesplit[2]
	A = "A: " + str(linesplit[4].count("A") + linesplit[4].count("a"))
	T = "T: " + str(linesplit[4].count("T") + linesplit[4].count("t"))
	C = "C: " + str(linesplit[4].count("C") + linesplit[4].count("c"))
	G = "G: " + str(linesplit[4].count("G") + linesplit[4].count("g"))

	if ref == "A":
		A = "A: " + str(linesplit[4].count(".") + linesplit[4].count(","))
	elif ref == "T":
		T = "T: " + str(linesplit[4].count(".") + linesplit[4].count(","))
	elif ref == "C":
		C = "C: " + str(linesplit[4].count(".") + linesplit[4].count(","))
	elif ref == "G":
		G = "G: " + str(linesplit[4].count(".") + linesplit[4].count(","))
	else:
		continue
    
	output = "\t".join([chrom, position, A, T, C, G])
	print(output)
```

上面的输出结果

```
chr1	1111111       A: 10    T: 1    C: 0    G: 0
```