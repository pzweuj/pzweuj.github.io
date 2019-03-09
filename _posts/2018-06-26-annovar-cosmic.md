---
title: 新建用于annovar的cosmic数据库
tags: default
---

众所周知，annovar提供的cosmic数据库到v70就没了，原因是cosmic开始区分学术账户和商业账户，开始限制下载。

目前cosmic的最新版本是v85。
annovar还是提供了制作cosmic数据库的方法：

第一步，当然先去cosmic下载数据：
[Cosmic Download](https://cancer.sanger.ac.uk/cosmic/download)

前提是，你需要有一个拥有下载权限的账户。
annovar需要的是以下这4个文件：

CosmicCodingMuts.vcf.gz
CosmicNonCodingVariants.vcf.gz
CosmicMutantExport.tsv.gz
CosmicNCV.tsv.gz

我下载了b37版本的，下载下来之后，先解压
```
gunzip *.gz
```

第二步，使用annovar提供的脚本：

[prepare_annovar_user.pl](http://www.openbioinformatics.org/annovar/download/prepare_annovar_user.pl)

```
prepare_annovar_user.pl -dbtype cosmic \
	CosmicMutantExport.tsv \
	-vcf CosmicCodingMuts.vcf \
	> hg19_cosmic85_coding.txt

prepare_annovar_user.pl -dbtype cosmic \
	CosmicNCV.tsv \
	-vcf CosmicNonCodingVariants.vcf \
	> hg19_cosmic85_noncoding.txt
```



[T_T]:要走要走啦