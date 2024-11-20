---
title: GATK 检测 Germline CNV
tags: software
---


适用于GATK 4.2 版本以上，流程参考于[这篇](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants)文章。



## 软件安装

GATK建议使用conda部署环境，其中gatkcondaenv.yml文件可以在下载的[gatk.zip](https://github.com/broadinstitute/gatk/releases/tag/4.2.2.0)中找到。
```bash
conda env create -f gatkcondaenv.yml
conda activate gatk
```

或者使用GATK官方docker
```bash
docker pull broadinstitute/gatk:4.2.2.0
```

测试了一下，普通环境下主要是缺少gcnvkernel这个模块。



## 下载测试数据

测试数据照旧使用xhmm的数据。这些数据是使用hs37d5参考基因组进行比对，经过了排序、去重以及GATK BQSR的标准bam文件准备流程。

```bash
wget https://statgen.bitbucket.io/xhmm/EXAMPLE_BAMS.zip
unzip EXAMPLE_BAMS.zip
```



## 运行

GATK文章中分为批次模式和单样本模式。其中单样本模式需要先使用批次模式建立model。这里使用批次模式进行。

需要输入参考基因组、bed文件、输出文件夹、bam文件存放文件夹（仅包含一批次的bam及bai）、ploidy文件。

其中ploidy文件格式参考

| CONTIG_NAME | PLOIDY_PRIOR_0 | PLOIDY_PRIOR_1 | PLOIDY_PRIOR_2 | PLOIDY_PRIOR_3 |
| ----------- | -------------- | -------------- | -------------- | -------------- |
| 1           | 0.01           | 0.01           | 0.97           | 0.01           |
| 2           | 0.01           | 0.01           | 0.97           | 0.01           |
| X           | 0.01           | 0.49           | 0.49           | 0.01           |
| Y           | 0.50           | 0.50           | 0.00           | 0.00           |



运行参数

```bash
reference=$1
bed=$2
output=$3
bamDir=$4
ploidy=$5
```



### 准备interval list

使用PreprocessIntervals准备intervallist
```bash
gatk PreprocessIntervals \
    -R $reference \
    -L $bed \
    --bin-length 0 \
    --padding 0 \
    -imr OVERLAPPING_ONLY \
    -O $output/target.prep.interval_list
```



### 注释GC，可选

使用AnnotateIntervals获得GC比例
```bash
gatk AnnotateIntervals \
    -L $bed \
    -R $reference \
    -imr OVERLAPPING_ONLY \
    -O $output/target.annotated.tsv
```



### 获得counts结果

使用循环，CollectReadCounts获得每个bam的counts
```bash
for i in $bamDir/*.bam;
do
    bamName=$(basename "$i" .bam);
    echo $bamName;

    # collect reads per bin
    gatk CollectReadCounts \
        -L $output/target.prep.interval_list \
        -R $reference \
        -imr OVERLAPPING_ONLY \
        -I $i \
        --format TSV \
        -O $output/counts/$bamName.tsv;    
done
```



### 获得query

获取所有counts结果，形成query供后续使用
```bash
ls $output/counts/*.tsv | while read line;
do
    echo -n "-I $line " >> $output/tmp;
done
sampleQuery=$(cat $output/tmp)
```



### 过滤interval，可选

可选操作，使用FilterIntervals与GC比例结果过滤interval，但是这步我运行报错了，因此没有选。
```bash
gatk FilterIntervals \
    -L $output/target.prep.interval_list \
    --annotated-intervals $output/target.annotated.tsv \
    $sampleQuery \
    -imr OVERLAPPING_ONLY \
    -O $output/cohort.gc.filtered.interval_list
```



### 获取ploidy

使用DetermineGermlineContigPloidy获取ploidy，如果上面运行成功，这一步的interval_list可用上一步的结果。
```bash
gatk DetermineGermlineContigPloidy \
    -L $output/target.prep.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    $sampleQuery \
    --contig-ploidy-priors $ploidy \
    --output $output/ploidy \
    --output-prefix ploidy \
    --verbosity DEBUG
```



### 使用批次模式获得结果

获得批次的CNV分析结果
```bash
gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L $output/target.prep.interval_list \
    $sampleQuery \
    --contig-ploidy-calls $output/ploidy/ploidy-calls \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output $output/cohort \
    --output-prefix cohort \
    --verbosity DEBUG
```



### 结果导出

使用PostprocessGermlineCNVCalls导出结果。GATK这个流程的问题在于结果居然是凭index而不是样本名称来导出的，需要确认样本的index才能导出对应样本。这里我直接重命名了。
```bash
echo -e "#SM\tSAMPLEID" > $output/sample.list
for ((i=0; i<$((`ls $bamDir | wc -l` / 2)); i++));
do
    echo "sample_${i}";
    gatk PostprocessGermlineCNVCalls \
        --model-shard-path $output/cohort/cohort-model \
        --calls-shard-path $output/cohort/cohort-calls \
        --sample-index $i \
        --contig-ploidy-calls $output/ploidy/ploidy-calls \
        --output-genotyped-intervals $output/output/sample_${i}_cohort.vcf.gz \
        --output-genotyped-segments $output/output/sample_${i}_segment.cohort.vcf.gz \
        --output-denoised-copy-ratios $output/output/sample_${i}_ratio.txt;

    # rename / name from bam SM tag
    sampleName=`zcat $output/output/sample_${i}_cohort.vcf.gz | grep "#CHROM" | cut -f 10`;
    echo -e $sampleName"\tSample_"$i >> $output/sample.list;
    mv $output/output/sample_${i}_cohort.vcf.gz $output/output/${sampleName}_cohort.vcf.gz;
    mv $output/output/sample_${i}_cohort.vcf.gz.tbi $output/output/${sampleName}_cohort.vcf.gz.tbi;
    mv $output/output/sample_${i}_segment.cohort.vcf.gz $output/output/${sampleName}_segment.cohort.vcf.gz;
    mv $output/output/sample_${i}_segment.cohort.vcf.gz.tbi $output/output/${sampleName}_segment.cohort.vcf.gz.tbi;
    mv $output/output/sample_${i}_ratio.txt $output/output/${sampleName}_ratio.txt;
done
```



最终结果可查看${sampleName}_segment.cohort.vcf.gz。总的来说，GATK流程比exomeDepth以及XHMM都要慢，而且CPU占用高，另外获得结果更多（考虑exomeDepth及XHMM结果一致，GATK的假阳性结果可能多得多）。
