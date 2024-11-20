---
title: GATK推荐的找体细胞突变流程
tags: software
---
>GATK4推荐流程。
>其实是Mutect2的使用教程。
>我还没用过。。
>放假前的最好一更啦。

#1 首先把原始数据处理成可以用 的bam
参考推荐的数据准备流程。

#2 如果是单肿瘤组织测序。
像这样。直接就从bam得到vcf了。需要参考基因组文件，还有要关注的区域如chr17plus.interval_list。-tumor后面是名字。
```
gatk Mutect2 \
-R ~/Documents/ref/hg38/Homo_sapiens_assembly38.fasta \
-I HG00190.bam \
-tumor HG00190 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L chr17plus.interval_list \
-O 3_HG00190.vcf.gz
```

#3 如果是配对测序，就是普通组织和肿瘤组织。
首先要建立PoN，就是Panel of Normals。
对每个普通组织的bam做一次下面的。
```
gatk Mutect2 \
-R ~/Documents/ref/hg38/Homo_sapiens_assembly38.fasta \
-I HG00190.bam \
-tumor HG00190 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L chr17plus.interval_list \
-O 3_HG00190.vcf.gz
```
然后
```
gatk CreateSomaticPanelOfNormals \
-vcfs 3_HG00190.vcf.gz \
-vcfs 4_NA19771.vcf.gz \
-vcfs 5_HG02759.vcf.gz \
-O 6_threesamplepon.vcf.gz
```
合起来，就得到PoN 6_threesamplepon.vcf.gz了。

#4 然后是配对组织的vcf生成。
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

出来的都是没过滤过的vcf。

#5 过滤掉污染
```
gatk GetPileupSummaries \
-I tumor.bam \
-V resources/chr17_small_exac_common_3_grch38.vcf.gz \
-O 7_tumor_getpileupsummaries.table

gatk CalculateContamination \
-I 7_tumor_getpileupsummaries.table \
-O 8_tumor_calculatecontamination.table

gatk FilterMutectCalls \
-V somatic_m2.vcf.gz \
--contamination-table tumor_calculatecontamination.table \
-O 9_somatic_oncefiltered.vcf.gz
```

#6 可选操作
```
gatk CollectSequencingArtifactMetrics \
-I tumor.bam \
-O 10_tumor_artifact \
–-FILE_EXTENSION ".txt" \
-R ~/Documents/ref/hg38/Homo_sapiens_assembly38.fasta

gatk FilterByOrientationBias \
-A G/T \
-A C/T \
-V 9_somatic_oncefiltered.vcf.gz \
-P tumor_artifact.pre_adapter_detail_metrics.txt \
-O 11_somatic_twicefiltered.vcf.gz
```


[T_T]:放假啦，过年啦
