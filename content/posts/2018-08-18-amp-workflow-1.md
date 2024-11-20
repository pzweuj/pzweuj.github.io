---
title: 扩增子流程以及复现文章（1）
tags: default
---
这一篇其实和之前的DADA2那个差不多，但是为了和后面的保持完整，所以再写一次。

这次使用的数据是一个位于太平洋底部约3公里的水下山脉，它是一个低温（~5-10°C）的热液喷发点。
该扩增子数据集是从山上收集的碎玄武岩中提取的DNA生成的，目的是开始研究深海岩石的微生物群落。
参考文献：[26779122](https://www.ncbi.nlm.nih.gov/pubmed/26779122)
使用的是Illumina MiSeq平台，测的是16S V4序列。
共有20个样品，其中4个是空白对照。

这里除了DADA2，还要使用一个叫[BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)的软件。

下载原始数据
---

```bash
curl -L -o dada2_amplicon_ex_workflow.tar.gz https://ndownloader.figshare.com/files/11342996
```
解压还有放好。
创建一个样本名称的列表。
```bash
ls *_R1.fq | cut -f1 -d "_" > samples.txt
```

去除引物
---

下载的压缩包中有primers.fa文件，是引物序列。接下来就是使用BBTools来比对引物序列然后去除掉。
```bash
for sample in $(cat ../Rawdata/samples.txt); \
do bbduk.sh \
	in=../Rawdata/"$sample"_sub_R1.fq \
	in2=../Rawdata/"$sample"_sub_R2.fq \
	out=../Filtdata/"$sample"_sub_R1_trim.fq.gz \
	out2=../Filtdata/"$sample"_sub_R2_trim.fq.gz \
	literal=GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT \
	k=10 \
	ordered=t \
	mink=2 \
	ktrim=l \
	rcomp=f \
	minlength=220 \
	maxlength=280 \
	tbo=t \
	tpe=t; \
done 2> ../Filtdata/bbduk_primer_trim.txt
```
这里的literal参数后面跟的是正向引物序列和反向引物序列，就从primers.fa文件上的；
k是需要找的kmer大小；
ordered，是否保持reads的顺序不变；
mink，检测读数的大小的最小值；
ktrim，要修建那一边，这里的l就是左边的意思；
rcomp，是否寻找互补的kmer；
tbo，是否过滤overlap序列；
tpe，是否将正向与反向修剪成同样的长度。

dada2
---

```R
library(dada2) # loads DADA2

list.files("Filtdata")

# setting a few variables we're going to use
fnFs <- sort(list.files("Filtdata", pattern="_sub_R1_trim.fq.gz", full.names=TRUE))
fnRs <- sort(list.files("Filtdata", pattern="_sub_R2_trim.fq.gz", full.names=TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path("Cleandata", paste0(sample.names, "_sub_R1_filtered.fq.gz"))
filtRs <- file.path("Cleandata", paste0(sample.names, "_sub_R2_filtered.fq.gz"))

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)
plotQualityProfile(fnFs[17:20])


## 过滤
out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs, truncLen=c(250, 200),
  maxN=0, maxEE=c(2, 2), truncQ=2, rm.phix=TRUE, minLen=175,
  compress=TRUE, multithread=TRUE # 在windows下，multithread设置成FALSE
)
head(out)

class(out) # matrix
dim(out) # 20 2

plotQualityProfile(filtFs)
plotQualityProfile(filtRs)
plotQualityProfile(filtFs[17:20])


## 计算错误率模型
# 分别计算正向和反向
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# 画出错误率统计图
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##消除误差
# 下面其实是一个批量的操作，如果是处理大文件，内存可能不足，更好的做法是一个一个样本的进行
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# 用sample.names来改名
names(derepFs) <- sample.names
names(derepRs) <- sample.names
## dada2核心算法
# 从头OTU方法必须在处理之前对样本进行聚类，因为样本之间没有聚类标签不一致且无法比较，即样本1中的OTU1和样本2中的OTU1可能不相同。
# 而DADA2可以更精确地解析序列变异，可以独立处理然后组合。
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")
# 检查，这里是检测正向的第一个样本
dadaFs[[1]]

## 合并双端
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, trimOverhang=TRUE, minOverlap=170)
# 构造列表
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# 检测长度分布
table(nchar(getSequences(seqtab)))

## 去除嵌合体
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

## 物种分类
# 这里是使用silva数据库进行注释。
# 然后把特征序列都提出来，在对这些序列进行命名，命名的方式为“ASV_x”。
taxa <- assignTaxonomy(seqtab.nochim, "Database/silva_nr_v132_train_set.fa.gz", multithread=T, tryRC=T)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

## 输出DADA2结果
# 这一步输出三个文件，一个是把特征序列都放在一起的fasta文件，一个是注释文件，一个是counts数文件。
# fasta:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "DADA2/ASVs.fa")

# count table:
asv_count <- t(seqtab.nochim)
row.names(asv_count) <- sub(">", "", asv_headers)
write.table(asv_count, "DADA2/ASVs_counts.txt", sep="\t", quote=F)

# tax table:
asv_taxa <- taxa
row.names(asv_taxa) <- sub(">", "", asv_headers)
write.table(asv_taxa, "DADA2/ASVs_taxonomy.txt", sep="\t", quote=F)
```

[-_-]:七夕