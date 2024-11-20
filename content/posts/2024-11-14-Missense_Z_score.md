---
title: 错义突变Z评分
tags: default
---

### 数据库

错义突变Z评分表示对于给定基因，观察到的错义突变（单个氨基酸取代）的预期数量的偏差，其中正值表示错义突变缺失，负值表示错义突变富集。

需要将这个内容注释到结果表中中，首先从UCSC下载该数据库。

```bash
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/gnomAD/pLI/missenseByTranscript.v4.1.bb
```

这个数据库提供了转录本水平的Z-scores，依据的GnomAD版本是v4.1。

也可以下载基因水平的Z-scores，但依据的GnomAD版本是v2.1.1，相对较旧。

```bash
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/gnomAD/pLI/missenseByGene.bb
```

因为我们最终是按转录本或基因将这个Z-scores匹配到表格中的，因此不用考虑是hg19还是hg38。

### 处理

下载下来是UCSC的bigbed格式，为了方便之后使用，先转换为bed格式。可以使用UCSC的[bigBedToBed工具](https://hgdownload.cse.ucsc.edu/admin/exe)进行转换，注意下载自己适用的版本。

```bash
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod a+x bigBedToBed
```

然后进行转换
```bash
bigBedToBed missenseByTranscript.v4.1.bb missenseByTranscript.hg38.v4.1.bed
```

最后使用自己的方案将Z-scores根据基因/转录本进行注释。

### 注释

如果你需要的是转录本水平，还有是使用VEP进行注释，可以使用[我的插件](https://github.com/pzweuj/VEP_Plugins_Self)**MissenseZscoreTranscript**，会注释出一个MissenseZscore字段（依赖Feature字段进行注释）。

```bash
vep -i input.vcf -o output.vcf --plugin MissenseZscoreTranscript,/path/to/missenseByTranscript.hg38.v4.1.bed
```

