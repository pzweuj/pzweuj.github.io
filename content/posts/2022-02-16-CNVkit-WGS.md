---
title: CNVkit分析WGS
tags: software
---

[CNVkit](https://cnvkit.readthedocs.io/en/stable/)一般用来分析肿瘤样本的拷贝数变异（使用配对样本或者正常样本建立参考基线的）。实际上，CNVkit也提供了全基因组胚系CNV分析的方法。

一般来说，WGS遗传样本不会做参考样本（也有会用同批次其他WGS样本作为参考的），同时分析多个样本时，运行命令如下
```bash
cnvkit.py batch \
	sample1.bam sample2.bam sample3.bam \
	-m wgs -f reference.fa \
	--annotate refFlat.txt \
	-t target.bed --target-avg-size 1000 \
	-p 16 -d output_dir \
	--segment-method hmm -n
```

其中，annotate参数需要输入一个注释文件，可以是refFlat格式。refFlat文件可来源于UCSC。如hg19的refFlat可在[这里](https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/)找到。target文件非必须，但是最好还是加入target来提高WGS分析的速度。target-avg-size是划分bin的大小，划得越小时，假阳性会越多；越大假阴性越多。WGS胚系分析[建议](https://github.com/etal/cnvkit/issues/484)使用hmm作为segment方法，当有对照样本时，在-n后指定，无对照样本则-n后留空。

作图请参考[官方文档](https://cnvkit.readthedocs.io/en/stable/plots.html)。

散点图
```bash
cnvkit.py scatter \
	sample1.cnr -s sample1.cns \
	-i sample1 --segment-color red -g BRAF --by-bin \
	-o sample1.BRAF.pdf
```

热图
```bash
cnvkit.py heatmap \
	*.cns -o samples.heatmap.pdf
```

转换为vcf格式
```bash
cnvkit.py export vcf sample1.cns > sample1.cnvkit.vcf
```