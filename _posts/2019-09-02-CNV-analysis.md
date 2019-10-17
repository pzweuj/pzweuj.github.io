---
title: CNV分析
tags: software
---

试用了几个分析CNV的软件。

### VarScan2
[VarScan](http://varscan.sourceforge.net/copy-number-calling.html)这个做somatic变异检测的软件也加入了对CNV分析的支持。

```bash
samtools mpileup -q 1 -f ref.fa normal.bam tumor.bam | \
	java -jar VarScan.jar copynumber prefix --mpileup 1 --data-ratio 1
```
如果normal和tumor的数据差异较大，记得调整--data-ratio，默认是1，可选范围为0-1。

以上这一步会生成prefix.copynumber文件，接下来进行分析
```bash
java -jar VarScan.jar copyCaller prefix.copynumber --output-file output.cnv.txt
```

输出的结果比较多，官网上建议用DNAcopy包进行过滤
```R
library(DNAcopy)
cn <- read.table("output.cnv.txt", header=F)
CNA.object <- CNA(cbind( cn$adjusted_log_ratio), cn$chrom,cn$chr_start, data.type="logratio", sampleid="tumor")
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
segs2 = segs$output
write.table(segs2[,2:6], file="out.file", row.names=F, col.names=F, quote=F, sep="\t")
```

### cnvkit
网上搜[cnvkit](https://cnvkit.readthedocs.io/en/stable/)好像一般是用来做大批量数据，这里还是用单个配对样本试一下。
```bash
cnvkit.py access hg19.fa -o access.hg19.bed

cnvkit.py batch tumor.bam --normal normal.bam \
	--target target.bed \
	--annotate refFlat.txt \
	--fasta hg19.fa \
	--access access.hg19.bed \
	--output-reference my_reference.cnn \
	--output-dir results \
	--diagram --scatter
```
这里的[refFlat.txt](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz)可以在官网下载。生成的my_reference.cnn实际上是参考样本组成的，方便以后使用。


### GATK4
我觉得[GATK4](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11147)的流程比较复杂，不是特别好用。
```bash
gatk CreateSequenceDictionary -I hg19.fa -O hg19.dict
gatk BedToIntervalList -I target.bed -O target.list -SD hg19.dict

gatk CollectReadCounts \
	-I normal.bam \
	-L target.list \
	--interval-merging-rule OVERLAPPING_ONLY \
	-O normal.counts.hdf5
	
gatk AnnotateIntervals \
	-R hg19.fa \
	-L target.list \
	--interval-merging-rule OVERLAPPING_ONLY \
	-O annotated_intervals.tsv
	
gatk CreateReadCountPanelOfNormals \
	-I normal.hdf5 \
	--annotated-intervals annotated_intervals.tsv \
	-O cnv.pon.hdf5
	
gatk CollectReadCounts \
	-I tumor.bam \
	-L target.list \
	--interval-merging-rule OVERLAPPING_ONLY \
	-O tumor.counts.hdf5
	
gatk CollectAllelicCounts \
	-I normal.bam \
	-R hg19.fa \
	-L target.list \
	-O normal.allelicCounts.tsv
	
gatk CollectAllelicCounts \
	-I tumor.bam \
	-R hg19.fa \
	-L target.list \
	-O tumor.allelicCounts.tsv

gatk DenoiseReadCounts \
	-I tumor.counts.hdf5 \
	--count-panel-of-normals cnv.pon.hdf5 \
	--standardized-copy-ratios tumor.standardizedCR.tsv \
	--denoised-copy-ratios tumor.denoisedCR.tsv
	
gatk ModelSegments \
	--denoised-copy-ratios tumor.denoisedCR.tsv \
	--allelic-counts tumor.allelicCounts.tsv \
	--normal-allelic-counts normal.allelicCounts.tsv \
	--output-prefix tumor \
	-O output_dir
```
跑上面最后一步的时候还可能报hdf5错误，需要重装（修复）hdf5才能继续。

```bash
gatk CallCopyRatioSegments -I tumor.cr.seg \
	-O tumor.called.seg
```

附带三个作图方法：

[PlotDenoisedCopyRatios](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_copynumber_plotting_PlotDenoisedCopyRatios.php)

[PlotModeledSegments](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_copynumber_plotting_PlotModeledSegments.php)

[PostprocessGermlineCNVCalls](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_copynumber_PostprocessGermlineCNVCalls.php)

### conifer
[conifer](http://conifer.sourceforge.net/)是外显子常用的cnv分析工具。这个软件写的很烂，感觉就算不是版本问题，也会出现很多bug。
```bash
python conifer.py rpkm --probes probes.txt --input sample.bam --output sample.rpkm.txt

python conifer.py analyze --probes probes.txt --rpkm_dir /RPKM/ \
	--output analysis.hdf5 \
	--svd 6 \
	--write_svals singular_values.txt \
	--plot_scree screeplot.png \
	--write_sd sd_values.txt
```
conifer使用的一些包年代久远，一些方法早就已经修改，存在大量[bug](https://www.biostars.org/p/262148/)，参考文章，对于依赖来说最后安装的版本是pytables 2.4.0、hdf5 1.8、numpy 1.9.3、numexpr 2.1。另外conifer_functions.py文件里第113行将samples\[s\]改成rpkm_filename。另外，probes.txt可以在官网中[下载](http://sourceforge.net/projects/conifer/files/probes.txt/download)，也可以自行创建，就是bed文件。
conifer需要8个样本以上进行对比，才能有结果。

### freec
[freec](http://boevalab.inf.ethz.ch/FREEC/)是分析全基因组cnv的工具。
freec首先需要编辑配置文件。
配置文件的格式如下
```
[configtype]
a = xxx
b = xxxx
c = xxxxxx
```
按上面这个格式，参照[这里](http://boevalab.inf.ethz.ch/FREEC/tutorial.html#CONFIG)进行配置文件的编辑。
最后使用
```bash
freec -conf config.txt
```
就可以运行了。