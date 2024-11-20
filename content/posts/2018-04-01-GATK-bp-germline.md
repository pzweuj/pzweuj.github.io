---
title: GATK4推荐的call生殖细胞突变流程（脚本）
tags: coding
---
这是整理过的脚本！

说一说我的目录结构：
```
|--~
   |--Project             # 存放项目
   |  |--Germline         # 每个项目单独一个文件夹
   |     |--Bam           # 存放最终生成的bam文件
   |     |--Bin           # 存放使用的脚本
   |     |--Temp          # 存放中间文件，最终可以删掉
   |     |--Vcf           # 存放最终的Vcf文件
   |
   |--Database            # 数据库
   |--Scripts             # 脚本大杂烩，做啥都套模板就行了
   |--Software            # 软件	
```


GATK4 best practise germline
```shell
# bam文件
sample=~/Project/Germline/Bam/sample.bam

# 需要的数据集
reference=~/Database/GATK/library/b37/human_g1k_v37_decoy.fasta

# 软件
gatk=~/Software/GenomeAnalysisTK/gatk-4.0.2.1/gatk

# 建立panel of normal，注意，对于配对样本，-I后的是正常样本，-tumor后的是正常样本的RGSM，如果有多个正常样本，重复多次此步骤
mkdir ../Temp

$gatk --java-options "-Xmx4g" HaplotypeCaller  \
	-R $reference \
	-I $sample \
	-O ../Temp/output.g.vcf.gz \
	-ERC GVCF

$gatk --java-options "-Xmx4g" GenotypeGVCFs \
	-R $reference \
	-V ../Temp/output.g.vcf.gz \
	-O ../Vcf/output.vcf.gz

rm -rf ../Temp
```


[T_T]:#疯了