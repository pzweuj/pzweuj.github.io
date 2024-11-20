---
title: maftools瀑布图上侧柱形图改为Multi Hits
tags: coding
---

maftools的瀑布图的上、左、右三个柱形图（条形图）都可以拿自己的数据进行修改，对于左右两侧的条形图，只需要基因可以匹配得上即可；对于上侧的柱形图，只需要样本可以匹配的上即可。参考[这篇文章](https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html)。

以TCGA-LAML数据为例，首先读入maftools内置的TCGA数据

```R
library(maftools)

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml.clin <- system.file('extdata', "tcga_laml_annot.tsv", package = "maftools")
laml <- read.maf(maf = laml.maf, clinicalData = laml.clin, verbose = FALSE)
oncoplot(maf = laml, draw_titv = TRUE)
```

这时的瀑布图默认上侧是TMB
![onco1](https://github.com/pzweuj/pzweuj.github.io/raw/master/content/data/images/maftools_oncoplot_mh1.png)


接下来，需要一个记录了multi hit的dataframe，必须包含样本列，具体做法是统计同一个样本中，同一个基因出现的次数，出现大于等于2就认为该基因是Multi Hit，然后再将Multi Hit基因的数量加起来，形成结果表。

```R
multi_hit_dt <- laml@data[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)]
multi_hit_dt_count <- multi_hit_dt[, Multi_Hit := .N, by = .(Hugo_Symbol, Tumor_Sample_Barcode)]
multi_hit_dt_count <- unique(multi_hit_dt_count)
```

接下来，可以将这个data.table转换为dataframe，因为dataframe对样本的顺序无要求，而data.table作为输入时，需要顺序一致。最后画出以Multi Hit为上侧柱形图的瀑布图。

```R
multi_hit_df <- as.data.frame(multi_hit_dt)
oncoplot(
  maf = laml,
  topBarData = multi_hit_df
)
```

![onco2](https://github.com/pzweuj/pzweuj.github.io/raw/master/content/data/images/maftools_oncoplot_mh2.png)

