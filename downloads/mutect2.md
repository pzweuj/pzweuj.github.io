用GATK4 Mutect2来找体细胞变异

软件：
gatk4
igv

准备数据：
配对的已经比对好的
NC.bam
TC.bam

1,使用Mutect2所有的参数来找体细胞变异
这里需要使用的是配对的NC.bam和TC.bam, panel of normals(PoN)和常用的生殖细胞变异位点。

```
gatk --java-options "-Xmx2g" Mutect2 \
-R hg38/Homo_sapiens_assembly38.fasta \
-I tumor.bam \
-I normal.bam \
-tumor HCC1143_tumor \
-normal HCC1143_normal \
-pon resources/chr17_pon.vcf.gz \
--germline-resource resources/chr17_af-only-gnomad_grch38.vcf.gz \
--af-of-alleles-not-in-resource 0.0000025 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L chr17plus.interval_list \
-O 1_somatic_m2.vcf.gz \
-bamout 2_tumor_normal_m2.bam 
```
最终会得到一个未过滤的体细胞变异文件1_somatic_m2.vcf.gz，一个重组的bam文件2_tumor_normal_m2.bam。
当然还会有两个索引文件。

2,创建PoN
首先，需要对所有的普通样本运行Mutect2的只有肿瘤模式。
例如，下面的这个-tumor 后面，其实是普通样本的名字。-I的其实是普通样本。
```
gatk Mutect2 \
-R ~/Documents/ref/hg38/Homo_sapiens_assembly38.fasta \
-I HG00190.bam \
-tumor HG00190 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L chr17plus.interval_list \
-O 3_HG00190.vcf.gz
```
接下来，把他们合并起来
```
gatk CreateSomaticPanelOfNormals \
-vcfs 3_HG00190.vcf.gz \
-vcfs 4_NA19771.vcf.gz \
-vcfs 5_HG02759.vcf.gz \
-O 6_threesamplepon.vcf.gz
```

这样就得到了PoN 6_threesamplepon.vcf.gz

3,有样本污染的情况。
对肿瘤样本使用GetPileupSummaries
主要需要一个参考的vcf。
```
gatk GetPileupSummaries \
-I tumor.bam \
-V resources/chr17_small_exac_common_3_grch38.vcf.gz \
-O 7_tumor_getpileupsummaries.table
```
然后计算一下污染
```
gatk CalculateContamination \
-I 7_tumor_getpileupsummaries.table \
-O 8_tumor_calculatecontamination.table
```

如果污染太多。。。
可以使用picard的CrosscheckFingerprints

4，过滤
```
gatk FilterMutectCalls \
-V somatic_m2.vcf.gz \
--contamination-table tumor_calculatecontamination.table \
-O 9_somatic_oncefiltered.vcf.gz
```

5,可选，过滤（可能）假的reads
```
gatk CollectSequencingArtifactMetrics \
-I tumor.bam \
-O 10_tumor_artifact \
–-FILE_EXTENSION ".txt" \
-R ~/Documents/ref/hg38/Homo_sapiens_assembly38.fasta
```

也可以用picard
```
java -jar picard.jar \
CollectSequencingArtifactMetrics \
I=tumor.bam \
O=10_tumor_artifact \
FILE_EXTENSION=.txt \
R=~/Documents/ref/hg38/Homo_sapiens_assembly38.fasta
```

然后使用FilterByOrientationBias过滤
```
gatk FilterByOrientationBias \
-A G/T \
-A C/T \
-V 9_somatic_oncefiltered.vcf.gz \
-P tumor_artifact.pre_adapter_detail_metrics.txt \
-O 11_somatic_twicefiltered.vcf.gz
```

