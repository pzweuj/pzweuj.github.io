---
title: GATK流程
tags: default
---
>不管啦，先放点代码上来，我真的在慢慢整理。

```
#!/usr/bin/bash
# 20171114
# FastQC查看质量

# 低质量数据过滤1
# -f 保留开始碱基个数，默认1；-l 保留结尾碱基个数，默认1
# -z 输出为gz
# fastx_trimmer -h -f -l -z -i test.fastq -o testQC

# 低质量数据过滤2
# -q 设置最小quality；-p 设置要保留的最少碱基百分比
# -z 输出为gz；-v 输出最终碱基个数
# fastq_quality_filter -h -q 30 -p 90 -z -i testQC -o finalQC -v

# 比对到参考基因组
# 生成sam文件
bwa index ucsc.hg19.fasta
bwa mem -M -t 1 ucsc.hg19.fasta read1.fastq.gz read2.fastq.gz | gzip -3 > readall.sam.gz

# bam文件预处理
# 重新排序
java -Xmx30g -jar picard.jar SortSam I=readall.sam.gz O=readall.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR='pwd'/tmp

# 去除duplicates
java -Xmx6g -jar picard.jar MarkDuplicates I=readall.bam O=readall_dedup.bam M=readall.metric VALIDATION_STRINGENCY=LENIENT TMP_DIR='pwd'/tmp

# 加头处理
java -Xmx6g -jar picard.jar AddOrReplaceReadGroups I=readall_dedup.bam O=readall_final.bam RGID=group1 RGLB=lib1 RGPU=unit1 RGPL=illumina RGSM=readall VALIDATION_STRINGENCY=LENIENT TMP_DIR='pwd'/tmp

# 生成索引
java -Xmx6g -jar picard.jar BuildBamIndex I=readall_final.bam O=readall_final.bai TMP_DIR='pwd'/tmp VALIDATION_STRINGENCY=LENIENT 

# 确认比对区域
java -Xmx6g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 1 -R ucsc.hg19.fasta -I readall_final.bam -o readall.intervals -known Mills_and_1000G_gold_standard.indels.hg19.vcf

# 对比对区域重新比对
java -Xmx6g -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I readall_final.bam -targetIntervalsreadall.intervals -o readall_realigned.bam -known Mills_and_1000G_gold_standard.indels.hg19.vcf

# 检测突变
java -jar -Xmx30g GenomeAnalysisTK.jar -T UnifiedGenotyper -R ucsc.hg19.fasta -L exome.bed -I readall_final.bam -glm BOTH -stand_emit_conf 10 -stand_call_conf 30 -o noclean.vcf

# 位点初步筛选
java -Xmx15g-jar GenomeAnalysisTK.jar -R ucsc.hg19.fasta -T SelectVariants -V noclean.vcf -ef -o clean.vcf

# 位点注释
mkdir anno
table_annovar.pl clean.vcf humandb/hg19/ -buildver hg19 -out anno/clean -remove -protocol refGene,genomicSuperDups,phastConsElements46way,esp6500siv2_all,exac03,1000g2015aug_eas,1000g2015aug_all,avsnp142,clinvar_20170905,scsnv,revel,mcap,cosmic68wgs,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

```


[T_T]:老井子最近都不说晚安。