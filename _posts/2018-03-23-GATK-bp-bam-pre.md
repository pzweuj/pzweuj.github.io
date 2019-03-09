---
title: GATK4推荐的数据预处理（脚本）
tags: coding
---
这是整理过的脚本！

说一说我的目录结构：
```
|--~
   |--Project             # 存放项目
   |  |--MakePreBam       # 每个项目单独一个文件夹
   |     |--Bam           # 存放最终生成的bam文件
   |     |--Bin           # 存放使用的脚本
   |     |--Temp          # 存放中间文件，最终可以删掉
   |     |--Rawdata       # 存放原始的fastq文件
   |
   |--Database            # 数据库
   |--Scripts             # 脚本大杂烩，做啥都套模板就行了
   |--Software            # 软件	
```


GATK4 best practise prepare bam
```shell
# 双端测序fq
sample1=~/Project/MakePreBam/Rawdata/Sample_1.fastq.gz
sample2=~/Project/MakePreBam/Rawdata/Sample_2.fastq.gz

# 需要的数据集
reference=~/Database/GATK/library/b37/human_g1k_v37_decoy.fasta
indel1=~/Database/GATK/library/b37/1000G_phase1.indels.b37.vcf
indel2=~/Database/GATK/library/b37/Mills_and_1000G_gold_standard.indels.b37.vcf

# 软件
bwa=~/Software/bwa-0.7.17/bwa
picard=~/Software/picard/picard.jar
gatk=~/Software/GenomeAnalysisTK/gatk-4.0.2.1/gatk

# 建立bwa的索引，如果建立过可以忽略
$bwa index $reference

# 比对到参考基因组
$bwa mem -t 6 $reference $sample1 $sample2 > ~/Project/MakePreBam/Temp/sample.sam

# 排序
$picard SortSam \
	I=~/Project/MakePreBam/Temp/sample.sam \
	O=~/Project/MakePreBam/Temp/sample.sorted.bam \
	SO=coordinate

# 去重复序列
$picard MarkDuplicates \
	I=~/Project/MakePreBam/Temp/sample.sorted.bam \
	O=~/Project/MakePreBam/Temp/sample.marked_dup.bam \
	M=~/Project/MakePreBam/Temp/sample_marked_dup.txt \
	REMOVE_DUPLICATES=false

# 加头信息。RGSM是样本名称，在somatic中需要，不能乱取。
$picard AddOrReplaceReadGroups \
	I=~/Project/MakePreBam/Temp/sample.marked_dup.bam \
	O=~/Project/MakePreBam/Temp/sample.addhead.bam \
	RGID=001 \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	RGSM=sample

# 创建索引
$picard BuildBamIndex \
	I=~/Project/MakePreBam/Temp/sample.addhead.bam


# indel区域重比对
$gatk BaseRecalibrator \
	-I ~/Project/MakePreBam/Temp/sample.addhead.bam \
	-R $reference \
	--known-sites $indel1 \
	--known-sites $indel2 \
	-O ~/Project/MakePreBam/Temp/sample_recal_data.table

# 确认比对区域重新校正
$gatk ApplyBQSR \
	-R $reference \
	-I ~/Project/MakePreBam/Temp/sample.addhead.bam \
	--bqsr-recal-file ~/Project/MakePreBam/Temp/sample_recal_data.table \
	-O ~/Project/MakePreBam/Bam/sample.final.bam

# 删除中间文件，如果要保留可以忽略
rm -rf ~/Project/MakePreBam/Temp
```


[T_T]:跟异性交往，要发乎情，止乎礼