---
title: 用cufflinks统计RPKM值
tags: software
---

一般的，做完rna-seq的比对部分，就需要找出每个基因或者转录本的counts/RPKM/FPKM值。
之前介绍过featureCounts，统计counts非常快。
这次要说的是[cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)，用来统计RPKM和FPKM值。

首先下载cufflinks
```bash
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar -zxvf cufflinks-2.2.1.Linux_x86_64.tar.gz
```
然后把路径加入到环境中。

由于cufflinks需要的bam文件必须是排序过的，所以在采取hisat2进行比对的流程后，必须用samtools进行排序。

用cufflinks进行统计：
```bash
cufflinks -p 8 \
    -g reference.gtf \
    -o output1_dir \
    input1.sort.bam
```
然后就能得到这个样本fpkm统计文件以及组装后的gtf文件。

可以对这些gtf进行合并
```bash
cuffmerge -o merge_output_dir \
    -p 8 \
    -g reference.gtf \
    -s reference.fa \
    GTF_list.txt
```
会生成一个merge.gtf文件。就是合并好的转录本。


事实上，我觉得使用cufflinks主要是为了fpkm和rpkm的统计，之后的差异分析，我更倾向于用DESeq2来做。

当然，目前我更喜欢的流程还是hisat2+featureCounts+DESeq2的流程。

最后，博客的更新又会变慢了。因为要好好复习了！

[-_-]:井我老婆