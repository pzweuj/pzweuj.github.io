---
title: ExomeDepth检测CNV
tags: software
---



[ExomeDepth](https://github.com/vplagnol/ExomeDepth)是一个基于HMM方法来检测全外显子CNV的软件。软件文档可在[这里](https://cran.r-project.org/web/packages/ExomeDepth/vignettes/ExomeDepth-vignette.pdf)查看。



## 安装exomeDepth

exomeDepth是一个R包，在R中安装

```R
install.packages("ExomeDepth")
```



## 下载测试数据

测试数据使用上次xhmm的数据。
```bash
wget https://statgen.bitbucket.io/xhmm/EXAMPLE_BAMS.zip
unzip EXAMPLE_BAMS.zip
```

这些数据是使用hs37d5参考基因组进行比对，经过了排序、去重以及GATK BQSR的标准bam文件准备流程。



## ExomeDepth流程

以下均是R脚本。



### 导入bam及bed

bed文件需要四列，无标题
```
<chrom>	<start>	<end>	<name>
```

导入数据
```R
library(ExomeDepth)

# 导入bed文件
bedFilePath <- "/path/to/EXOME.bed"
bedFile <- read.table(bedFilePath, head=FALSE)
names(bedFile) <- c("chromosome", "start", "end", "name")
bedFile <- as.data.frame(bedFile)

# 导入bam文件
bamFilePath <- "/path/to/bam/"
bamFile <- list.files(bamFilePath, pattern="*.bam$")
bamFile <- file.path(bamFilePath, bamFile)

# 输出文件夹
outputPath <- "/path/to/output/"
if(!(dir.exists(outputPath))) {
    dir.create(outputPath)
}
```



### 获得counts

对bam获得bed中每个区域的read counts。
```R
my.counts <- getBamCounts(bed.frame=bedFile,
    bam.files=bamFile,
    include.chr=FALSE
)
my.counts.dafr <- as(my.counts, "data.frame")
exomeCount.mat <- as.matrix(my.counts.dafr[, grep(names(my.counts.dafr), pattern="*.bam")])
```



### 注释

可以使用这一步加入注释，但是我后续运行时报错了，因此没有加入
```R
data(exons.hg19)
genes.GRanges.hg19 <- GenomicRanges::GRanges(
    seqnames=exons.hg19$chromosome,
    IRanges::IRanges(start=exons.hg19$start, end=exons.hg19$end),
    names=exons.hg19$name
)
```

同样的，加入常见CNV注释，后续执行也报错了
```R
data(Conrad.hg19)
```



### 循环执行

循环执行即每次挑出一个样本作为case，其余样本作为control，获得每个样本的结果。下面代码中使用AnnotateExtra注释时，报错。但是还可以自行注释，因此影响不大。循环执行的方式适合于在同一批次的数据中找CNV。如果是单样本的方式，则调整下面脚本，固定reference和choice即可。分析获得的结果与上次的xhmm结果相比，多了一个缺失。不过看了下支持reads数只有15条，因此可以排除，实际结果是相同的。
```R
nSamples <- ncol(exomeCount.mat)
for (i in 1:nSamples) {
    my.choice <- select.reference.set(
        test.counts=exomeCount.mat[, i],
        reference.counts=exomeCount.mat[, -i],
        bin.length=(my.counts.dafr$end - my.counts.dafr$start)/1000, n.bins.reduced=10000
    )
    my.reference.selected <- apply(X=exomeCount.mat[, my.choice$reference.choice, drop=FALSE], MAR=1, FUN=sum)
    all.exons <- new("ExomeDepth", test=exomeCount.mat[, i],
        reference=my.reference.selected,
        formula="cbind(test, reference) ~ 1"
    )
    all.exons <- CallCNVs(x=all.exons,
        transition.probability=10^-4,
        chromosome=my.counts.dafr$chromosome,
        start=my.counts.dafr$start,
        end=my.counts.dafr$end,
        name=my.counts.dafr$exon
    )
    # all.exons <- AnnotateExtra(x=all.exons, reference.annotation=Conrad.hg19.common.CNVs,
    #     min.overlap=0.5, column.name="Conrad.hg19"
    # )
    # all.exons <- AnnotateExtra(x=all.exons, reference.annotation=genes.GRanges.hg19,
    #     min.overlap=0.0001, column.name="exons.hg19"
    # )
    output.file <- paste(outputPath, "/", colnames(exomeCount.mat)[i], ".txt", sep="")
    write.table(file=output.file, x=all.exons@CNV.calls, row.names=FALSE, quote=FALSE, sep="\t")
}
```



## docker

建立了一个docker镜像
```bash
docker pull pzweuj/exomedepth:1.1.15
```

使用命令
```bash
docker run --rm -it -v /outsidePath:/insidePath pzweuj/exomedepth:1.1.15 exomeDepthPipe.R -i <bam dir> -o <output dir> -b <bed file>
```

