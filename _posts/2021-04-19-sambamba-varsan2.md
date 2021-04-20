---
title:  sambamba与varscan2的使用
tags: software
---


现在的肿瘤方向分析流程我用的是

比对：BWA

排序：samtools

去重：gatk Markduplicates

校正：gatk BaseRecalibrator + gatk ApplyBQSR

变异检测：gatk Mutect2



尝试一下另外一条路线

比对：BWA

排序：[sambamba](https://lomereiter.github.io/sambamba/)

去重：sambamba

校正：不做

变异检测：[varscan2](http://varscan.sourceforge.net/)



## sambamba

用sambamba的原因主要是因为比samtools快。

直接下载编译好的版本，解压就能用

```bash
wget https://github.com/biod/sambamba/releases/download/v0.8.0/sambamba-0.8.0-linux-amd64-static.gz
gunzip sambamba-0.8.0-linux-amd64-static.gz
```

常用的view与sort功能
```bash
sambamba view input.sam -S -h -f bam -o output.bam
sambamba sort output.bam -t 8 -o output.sort.bam
```

与bwa的pipe
```bash
bwa mem -t 8 -M -R "@RGxxxxx" ref.fa read1.fq.gz read2.fq.gz \
	| sambamba view /dev/stdin -S -h -f bam -o output.bam
```

sambamba去重和samtools去重，有空测试一下差异
```bash
sambamba markdup -t 8 -p input.bam output.bam
samtools markdup -@ 8 -O BAM input.bam output.bam
```



## VarScan2

VarScan2也是[直接下载jar包](https://sourceforge.net/projects/varscan/files/)就能用了。

VarScan2 somatic 貌似只支持配对样本。

```bash
samtools mpileup -B -f ref.fa -q 15 -d 10000 tumor.bam > tumor.pileup
samtools mpileup -B -f ref.fa -q 15 -d 10000 normal.bam > normal.pileup
java -jar VarScan2.jar somatic \
	normal.pileup tumor.pileup output.prefix \
	--min-coverage-normal 10 --min-coverage-tumor 20 \
	--min-var-freq 0.02 --strand-filter 1
```

pipe

```bash
samtools mpileup -B -f ref.fa \
	-q 15 -d 10000 normal.bam tumor.bam \
	| java -jar VarScan2.jar somatic \
	-mpileup output.prefix \
	--min-coverage-normal 10 --min-coverage-tumor 20 \
	--min-var-freq 0.02 --strand-filter 1
```

