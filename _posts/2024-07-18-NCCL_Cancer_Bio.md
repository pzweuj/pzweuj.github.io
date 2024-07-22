---
title: NCCL 肿瘤生信室间质评
tags: default
---

注：此篇文章介绍2024年卫健委肿瘤生信室间质评的分析方案，纯手工敲命令，未建立流程（意味着用相同版本的相同软件和相同命令能获得相同的结果哦）。

使用[Illumina平台的数据](https://pan.baidu.com/s/1i-KaQ8ViIgX8DelpPYXdvQ?pwd=47ny)作为例子。

检测要求：
1，检测出测序 Panel 中包含的所有体细胞突变(somatic mutation)；
2，Illumina 测序文件：需对突变频率VAF≥1%且原始数据测序深度(Read depth)>500X的位
点进行回报。

叠甲：
此篇文章在2024年7月18日发布，2024年7月17日 质评已截止回报。

## 数据处理

先对202411以及2024NC两个样本均使用相同的方案获得bam文件，以202411为例。使用更新的软件理论可获得更准确的效果或者更佳的效率，但是环境里有啥就用啥了。

```bash
mkdir cleandata bam vcf fusion cnv annotation QC
```

### 数据质控

使用fastp（0.14.1）对原始数据进行质控。
```bash
fastp \
	-i 202411_R1.fq.gz \
	-I 202411_R2.fq.gz \
	-o cleandata/202411_R1.fq.gz \
	-O cleandata/202411_R2.fq.gz \
	-j QC/202411.json \
	-h QC/202411.html \
	-w 16
```

### 比对

使用bwa（0.7.17-r1188）/ samtools（1.5）进行比对。用自己的hg19参考基因组就好，应该没人没有bwa吧。
```bash
bwa mem -t 16 -M -Y \
	-R "@RG\tPL:NCCL\tPU:illumina\tSM:202411\tID:mapping" \
	ucsc.hg19.fasta \
	cleandata/202411_R1.fq.gz cleandata/202411_R2.fq.gz | \
	samtools view -bSh | samtools sort -@ 16 - -o bam/202411.sort.bam
samtools index bam/202411.sort.bam
```

### 标记重复

使用gatk（4.5.0.0，下同）MarkDuplicates进行重复标记。
```bash
gatk MarkDuplicates \
	-I bam/202411.sort.bam \
	-O bam/202411.marked.bam \
	-M bam/202411.marked.txt \
	--CREATE_INDEX true
mv bam/202411.marked.bai bam/202411.marked.bam.bai
```

接下来本应进行BQSR的，但是考虑室间质评的数据质量足够好，省去这一步。预处理完成的bam就可以同时进行下面的分析了。

## 质控统计

使用bamdst（1.0.9）进行统计。
```bash
mkdir QC/202411
bamdst -p Illumina.bed bam/202411.marked.bam -o QC/202411
```

## 变异检测

使用gatk Mutect2 进行SNV/InDel检测。这里也不加入人群数据库和PoN了。常规生产中，还要使用-ip参数将bed放大一定大小来覆盖剪接区，这里也不用。

```bash
# gatk 4.5不需要指定tumor样本名称，默认都是tumor样本
gatk Mutect2 \
	-R ucsc.hg19.fasta \
	-I bam/202411.marked.bam \
	-I bam/2024NC.marked.bam \
	-normal 2024NC \
	--max-reads-per-alignment-start 0 \
	-L Illumina.bed \
	-O vcf/202411.vcf.gz \
	--native-pair-hmm-threads 16

gatk FilterMutectCalls \
	-R ucsc.hg19.fasta \
	-V vcf/202411.vcf.gz \
	-O vcf/202411.filtered.vcf.gz

# 我习惯最后做一个左对齐和切割多方向突变位点
gatk LeftAlignAndTrimVariants \
	-R ucsc.hg19.fasta \
	-V vcf/202411.filtered.vcf.gz \
	-O vcf/202411.filtered.left.vcf \
	--split-multi-allelics
```

使用bcftools（1.17）提取出PASS的位点，并进行深度和突变丰度过滤，满足室间质评要求。为了避免过滤过于严苛，不要直接把深度设置为500。不建议直接提取PASS的突变，这样可能会导致某些突变假阴性，根据NCCL的规则，去掉对照背景就好。

```bash
bcftools filter -i 'FORMAT/DP>=200 && FORMAT/AF>=0.01' -o vcf/202411.pass.vcf -O v vcf/202411.filtered.left.vcf
cat vcf/202411.pass.vcf | grep -v "normal_artifact" > vcf/202411.final.vcf
```

经过上述处理，获得29个体细胞突变。

下面可以各自用自己的[注释方案](https://pzweuj.github.io/2024/04/02/VEP.html)进行注释，[按照3'原则进行校正](https://pzweuj.github.io/2021/04/15/NCCL-pancancer.html)等等咯，VEP等工具在设置参数里可以指定3'原则，在结果中的HGVSg即是已校正的基因组坐标了。最终的结果是必然需要人工复核的，不要想着结果直出直接回报。经过核对，上述的29个突变均可回报。


## 融合检测

使用manta（1.6.0）进行融合检测。

```bash
# manta需要bgzip压缩的bed文件，当然，也可以选择不传入
bgzip Illumina.bed -c > Illumina.bed.gz
tabix -p bed Illumina.bed.gz

# 配置manta
configManta.py \
	--bam bam/2024NC.marked.bam \
	--tumorBam bam/202411.marked.bam \
	--exome --runDir fusion \
	--referenceFasta ucsc.hg19.fasta \
	--callRegions Illumina.bed.gz

# 运行
runWorkflow.py -j 16
zcat fusion/results/variants/somaticSV.vcf.gz > fusion/202411.manta.vcf

# 使用snpeff注释
java -jar snpEff.jar hg19 fusion/202411.manta.vcf > fusion/202411.manta.snpeff.vcf
```

## CNV检测

使用CNVkit（0.9.11）进行检测，需要先建立target和antitarget区域。anti region和refFlat都可以在ucsc找到。

```bash
# bed文件预处理
cnvkit.py access ucsc.hg19.fasta \
	-x wgEncodeDukeMapabilityRegionsExcludable.bed \
	-o access.hg19.bed
cnvkit.py target Illumina.bed \
	--annotate refFlat.txt \
	-o Illumina.target.bed \
	--short-names
cnvkit.py antitarget Illumina.target.bed \
	-g access.hg19.bed \
	-o Illumina.antitarget.bed

# 接下来获得coverage文件

cnvkit.py coverage bam/202411.marked.bam Illumina.target.bed -o cnv/202411.targetcoverage.cnn
cnvkit.py coverage bam/202411.marked.bam Illumina.antitarget.bed -o cnv/202411.antitargetcoverage.cnn
cnvkit.py coverage bam/2024NC.marked.bam Illumina.target.bed -o cnv/2024NC.targetcoverage.cnn
cnvkit.py coverage bam/2024NC.marked.bam Illumina.antitarget.bed -o cnv/2024NC.antitargetcoverage.cnn

# 将对照样本形成参考基线
cnvkit.py reference cnv/2024NC.targetcoverage.cnn cnv/2024NC.antitargetcoverage.cnn -f ucsc.hg19.fasta -o cnv/reference.cnn

# 分析CNV
cnvkit.py fix cnv/202411.targetcoverage.cnn cnv/202411.antitargetcoverage.cnn cnv/reference.cnn -o cnv/202411.cnr

# 获得每个划窗的拷贝数值
cnvkit.py call cnv/202411.cnr -o cnv/202411.cns

# 获得segment后的结果
cnvkit.py segment cnv/202411.cnr -o cnv/202411.seg.cnr
cnvkit.py call cnv/202411.seg.cnr -o cnv/202411.seg.cns
```

CNVkit输出的结果非常多，我建议是将划窗结果复制到excel表中，通过log2值计算确切的Copy Number，然后再手动筛选。

根据经验，筛选CN值大于4的外显子即可，从筛选结果中挑选回报基因和确认外显子，貌似未试过回报缺失。
