---
title: 评估WES流程SNP/INDEL检测准确率
tags: default
---

评估一套WES分析流程检测的准确率(Benchmark)，一般会与NA12878的基准结果进行对比。即先把NA12878的标准检出vcf，与自己分析流程获得的vcf，使用同样的bed先限制区域获得结果。

然后以NA12878.target.vcf的结果作为标准，与NA12878的交集即真阳性位点、NA12878中没有的检出即假阳性位点、NA12878有但自己的vcf没有的检出即假阴性位点。

那么，首先需要找到NA12878的标准检出结果vcf。

### Genome in a Bottle (GIAB)项目
GIAB项目是由美帝国家标准与技术研究院（NIST）领导的，专门提供高质量的基因组参考标准。可以在其[官方网站](https://www.nist.gov/programs-projects/genome-bottle)上找到NA12878的基准数据。也有用NA24149的，还有中国人群家系[NA24631、NA24694、NA24695]。

高置信的vcf在[NCBI的FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/)中，这里下载GRCh37的版本。

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
```

### Benckmark方法
GIAB项目同时提供了一个Benchmark的最佳实践。

最好是先把自己的bed和GIAB的bed取个交集，不然结果就差得远咯
```bash
bedtools intersect -a my.bed -b giab.bed > target.bed
```

#### hap.py

即Github中的[ga4gh/benchmarking-tools](https://github.com/ga4gh/benchmarking-tools/)这个项目。主要使用了Illumina的[hap.py](https://github.com/Illumina/hap.py)脚本进行评估。

以下是使用命令
```bash
export HGREF=reference.fasta
hap.py \
	PG.vcf.gz \
	NA12878.vcf.gz \
	-f target.bed \
	--threads 8 \
	-o output.prefix
```

在输出的csv文件中，会分别输出INDEL和SNP的Benchmark结果。


#### VCFeval

当然，也可以使用[RTG-tools](https://github.com/RealTimeGenomics/rtg-tools)里的VCFeval进行评估。但是它不会分开SNP和Indel，需要先用bcftools将SNP和Indel自行分开。

```bash
# 提取SNP和Indel参考示例，baseline和custom的vcf都需要进行这个操作
bcftools view -v snps NA12878.vcf.gz -o NA12878.snp.vcf.gz
bcftools view -v indels NA12878.vcf.gz -o NA12878.indels.vcf.gz

# 需要先对参考基因组创建一个SDF
rtg format -o RTG_GRCh37_SDF human_v37_decoy.fasta

rtg vcfeval \
	--bed-regions target.bed \
	--evaluation-regions target.bed \
	--template RTG_GRCh37_SDF \
	--baseline NA12878.vcf.gz \
	--calls PG.vcf.gz \
	--sample HG001,NA12878 \
	--output outpt.prefix
```


### 性能要求
华大主导的[人类全基因组遗传变异解读的高通量测序数据规范](https://en.genomics.cn/uploadfiles/2018/12/20181214171110779.pdf)中对各项性能的要求如下

|    类别     | 要求 |
| :---------: | ---- |
|  SNP精确度  | ≥99% |
|  SNP灵敏度  | ≥99% |
| InDel精确度 | ≥92% |
| InDel灵敏度 | ≥96% |



### 变异检测流程

不考虑商业软件或GPU软件的前提下，现阶段的胚系SNP/Indel变异检测常规方案是
1，使用BWA进行比对；

2，使用Samtools排序；

3，使用GATK去重及质量校正等；

4，使用GATK HaplotypeCaller进行变异检测。

也确实有性能表现更优的存在，例如Google的[DeepVariant](https://github.com/google/deepvariant)。不用DeepVariant的理由大概是GATK有更多的工具可以支持HaplotypeCaller的下游分析？

