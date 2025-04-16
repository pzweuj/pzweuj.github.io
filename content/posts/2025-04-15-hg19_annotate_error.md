---
title: hg19的错误注释
tags: default
---

最近收到了关于注释所得的HGVSc以及HGVSp不准确的反馈，排查问题主要源于GRCh37和GRCh38的版本差异。


## 问题

下面的这个突变位点，dbsnp的编号是[rs9332964](https://www.ncbi.nlm.nih.gov/snp/rs9332964)。

它的hg19坐标信息如下：

```hg19
2	31754395	.	C	T
```

它的hg38坐标信息如下：

```hg19
2	31529325	.	C	T
```

它在NCBI中的MANE转录本信息如下：

```txt
NM_000348.4:c.680G>A:p.Arg227Gln:4/5
```


但实际上，在参考基因组**hg19**下，我使用VEP注释所得的结果是：

```txt
NM_000348.4:c.677G>A:p.Arg226Gln:5/6
```

我使用SNPEff注释所得的结果是：
```txt
NM_000348.4:c.677G>A:p.Arg226Gln:4/5
```

--------

可以看到，SNPEff的HGVSc、HGVSp注释错误，VEP在这个基础上，还错了外显子编号。


## 分析

### 已知 1

通过UCSC Browser可以看到，hg19比hg38这个位置少了一个G，剩下的两个G无法构成氨基酸，导致氨基酸的编号少了一位。尽管都是在参考转录本NM_000348.4上，但hg19版本的转录本是liftover出来的，硬生生的把这两个G也抠了出来，认为这里是内含子（存疑）。

![pic1](https://github.com/pzweuj/pzweuj.github.io/raw/master/content/data/images/hg19_annotate_pic1.png)

最终，导致了在同一个转录本中，hg19与hg38注释不同的问题。

### 已知 2

为什么VEP注释所得的区域结果是 5/6（位于5号外显子，共计6个外显子），因为VEP的GRCh37使用的GFF是[GCF_000001405.25_GRCh37.p13_genomic.gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz)，这上面就是错的。而我的[ManeSelectBed](https://github.com/pzweuj/ManeSelectBed)项目是通过读取行数得到Exon数目的，也错了，有空得改改 🙂。


### 已知 3

然后，为什么HGVS是错的，从[SNPEff的信息](https://pcingola.github.io/SnpEff/snpeff/human_genomes/)里可以看出一点端倪（VEP应该同理）。

>Variant annotations using RefSeq may not be precise at the exact loci where the RefSeq transcript doesn't match the genome reference. This is yet another consequence of the previous items, but since the transcript do not match the reference genome, and variant annotations are based on the reference genome, the variant annotaion predictions might be off at those genomic loci.



根据NCBI对MANE转录本的描述，[MANE转录本与GRCh38是完全一致的](https://www.ncbi.nlm.nih.gov/refseq/MANE/)。

>MANE Select transcripts exactly match the sequence of the GRCh38 human reference genome assembly. 


因此，当我们使用GRCh37时，由于GRCh37可能与GRCh38存在差异，即与MANE转录本会存在差异，对于转录本来说，是发生了突变；HGVS的此时是根据区域的参考基因组坐标注释（GFF文件），和实际检出坐标重新计算的，导致错误。

### 已知 4

NCBI页面上的外显子编号好解析，因为它显示就是GRCh38优先。那么为什么SNPEff却能注释出正确的外显子编号呢。我的SNPEff注释是使用hg19版本，而不是GRCh37。我猜测（没时间做测试了）SNPEff的注释是拟合了GENCODE和NCBI，不太确定，但是GENCODE的GRCh37版本，尽管在liftover后也把前面exon1里缺失了一个P（30）的区域进行了切割，但实际的外显子编号标记是\_1\_0和\_1\_1。这也是混乱的地方，在参考基因组水平看这里有一个内含子，但在转录本水平却又在CDS上。

### 已知 5

现在hg19的这个位置，与hg38相比是一个缺失，即与MANE转录本相比也是一个缺失，所以实际是包含在转录本里的。如果单纯采用将\_1\_0和\_1\_1合并起来的方式，会造成的问题是，这个转录本少了1bp，除不了3了。

其他类似于插入的情况，就相当于转录本多了几bp。


## 问题修正

我不准备修正这个问题，这个问题产生自参考基因组落后于参考转录本，只要更新参考基因组就迎刃而解，编写补丁修正这个问题治标不治本。

内部需要做的是整理好这些相关资料，在必要时可以输出；固定好参考基因组，不要反复横跳造成多次检测相同位点回报不同的问题即可。

[已反馈给VEP官方](https://github.com/Ensembl/ensembl-vep/issues/1863)，看看他们有没有解决方案。

