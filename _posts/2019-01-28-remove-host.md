---
title: 宏基因组，除去宿主序列
tags: default
---

~~偷偷摸摸写成未来的日期，假装天天更新。~~

需要的软件：

[bowtie2](https://github.com/BenLangmead/bowtie2)

[samtools](https://github.com/samtools/samtools)

[bedtools](https://github.com/arq5x/bedtools2)


首先，是将序列比到宿主基因组上
```bash
# 建索引
bowtie2-build host.fna host

# 比对
bowtie2 -x host \
	-1 sample_1.fastq.gz \
	-2 sample_2.fastq.gz \
	-S sample.sam --threads 16

# 转成bam
samtools view -bS sample.sam > sample.bam
```

然后，除掉比对上了的
```bash
samtools view -b -f 12 -F 256 sample.bam > sample.unmapped.bam

# -f 表示提取； 12 表示未比对上的reads和未比对上的pair
# -F 表示不要提取； 256 主要比对上了的
```
那个数字，可以在这里查询： [flag](http://broadinstitute.github.io/picard/explain-flags.html)


最后，把bam转回去fastq。
```bash
# samtools根据名字排序
samtools sort -n sample.unmapped.bam -O BAM -o sample.unmapped.sort.bam

# bedtools 转格式
bedtools bamtofastq -i sample.unmapped.sort.bam \
	-fq sample_remove_host_1.fastq \
	-fq2 sample_remove_host_2.fastq
```


P.S.

bowtie2中其实有一个参数，可以直接把未比对上的输出出来。
```bash
bowtie2 -x host \
	-1 sample_1.fastq.gz \
	-2 sample_2.fastq.gz \
	-S sample.sam --threads 16 \
	--un-conc sample.unalign.sam
```



[-_-]:继续努力