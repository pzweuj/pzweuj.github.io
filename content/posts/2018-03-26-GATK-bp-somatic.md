---
title: GATK4推荐的call体细胞突变流程（脚本）
tags: coding
---
这是整理过的脚本！

说一说我的目录结构：
```
|--~
   |--Project             # 存放项目
   |  |--Somatic          # 每个项目单独一个文件夹
   |     |--Bam           # 存放最终生成的bam文件
   |     |--Bin           # 存放使用的脚本
   |     |--Temp          # 存放中间文件，最终可以删掉
   |     |--Vcf           # 存放最终的Vcf文件
   |
   |--Database            # 数据库
   |--Scripts             # 脚本大杂烩，做啥都套模板就行了
   |--Software            # 软件	
```


GATK4 best practise somatic
```shell
# normal和tumor样本
normal=~/Project/Somatic/Bam/normal.bam
tumor=~/Project/Somatic/Bam/tumor.bam

# 需要的数据集
reference=~/Database/GATK/library/b37/human_g1k_v37_decoy.fasta
indel1=~/Database/GATK/library/b37/1000G_phase1.indels.b37.vcf
indel2=~/Database/GATK/library/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
snp=~/Database/GATK/library/b37/1000G_phase1.snps.high_confidence.b37.vcf
af=~/Database/GATK/library/b37/af-only-gnomad.raw.sites.b37.vcf.gz

# 软件
gatk=~/Software/GenomeAnalysisTK/gatk-4.0.2.1/gatk

# 建立panel of normal，注意，对于配对样本，-I后的是正常样本，-tumor后的是正常样本的RGSM，如果有多个正常样本，重复多次此步骤
mkdir ../Temp

$gatk Mutect2 \
	-R $reference \
	-I $normal \
	-tumor normal \                       # 这里必须用RGSM
	--germline-resource $af \
	-O ../Temp/normal_for_pon.vcf.gz

# 用多个-vcfs合并成PON，如果只有一个正常样本，就用一个。
$gatk CreateSomaticPanelOfNormals \
	-vcfs ../Temp/normal_for_pon.vcf.gz \
	-O ../Temp/pon.vcf.gz

# 找变异
$gatk Mutect2 \
	-R $reference \
	-I $tumor \
	-I $normal \
	-tumor tumor \
	-normal normal \
	-pon ../Temp/pon.vcf.gz\
	--germline-resource $af \
	--af-of-alleles-not-in-resource 0.0000025 \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	-O ../Vcf/somatic_m2.vcf.gz \
	-bamout ../Bam/tumor_normal_m2.bam

# 一次过滤
$gatk GetPileupSummaries \
	-I $tumor \
	-V $snp \
	-V $indel1 \
	-V $indel2 \
	-O ../Temp/summaries.table

$gatk CalculateContamination \
	-I ../Temp/summaries.table \
	-O ../Temp/calculatecontamination.table

$gatk FilterMutectCalls \
	-V ../Vcf/somatic_m2.vcf.gz \
	--contamination-table ../Temp/calculatecontamination.table \
	-O ../Vcf/somatic_oncefiltered.vcf.gz

# 二次过滤
$gatk CollectSequencingArtifactMetrics \
	-I $tumor \
	-O ../Temp/tumor_artifact \
	-EXT ".txt" \
	-R $reference

$gatk FilterByOrientationBias \
	-V ../Vcf/somatic_oncefiltered.vcf.gz \
	--artifact-modes 'G/T' \
	-P ../Temp/tumor_artifact.pre_adapter_detail_metrics.txt \
	-O ../Vcf/somatic_twicefiltered.vcf.gz

# 最后删除中间文件
rm -rf ../Temp
```


[T_T]:#疯了