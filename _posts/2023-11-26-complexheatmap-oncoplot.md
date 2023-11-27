---
title: complexHeatmap 瀑布图
tags: coding
---

使用[maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html)来绘制瀑布图，很多内容没办法自定义，这里用[complexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)来画图。

还是先使用maftools来读入maf文件，这样需要调整的格式比较少，我把UTR和同义突变都纳入统计，实际不用
```R
library(ComplexHeatmap)
library(maftools)
library(dplyr)
library(tidyr)

# 读入数据
clin <- "clinical.txt"
# 需要统计的突变类型
nonSyn_list <- c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'In_Frame_Del',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Nonstop_Mutation',
  'Intron_Mutation',
  'Start_Lost_Mutation',
  'Other_Mutation',
  'Synonymous_Mutation',
  "5\\'UTR'",
  "3\\'UTR'"
)

all <- read.maf("sample.maf",
  clinicalData = clin,
  vc_nonSyn = nonSyn_list
)

# 可以使用下面命令把特殊的样本挑出来
all <- subsetMaf(all, clinQuery = "Cancer_Type != 'Lung'")
```

然后对multi hit进行统计。maftools统计multi hit的方式是同一个样本中的同一个基因有多个突变就会记为Multi Hit。由于我没有只留着体细胞突变，出现的Multi Hit基因特别多是必然的。

```R
df <- all@data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")]
df_count <- df %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  summarise(Variant_Classification = paste(Variant_Classification, collapse = ";"))
df_uniq <- df %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  summarise(Variant_Classification = paste(unique(Variant_Classification), collapse = ";"))
df_mh <- df %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  summarise(Variant_Classification = ifelse(n() > 1, "Multi_Hit", toString(Variant_Classification[1])))
```

整理一下临床信息
```R
clinical_data <- all@clinical.data
clinical_data$Age_tmp <- NULL
clinical_data <- t(clinical_data)
colnames(clinical_data) <- clinical_data[1,]
clinical_data <- clinical_data[-1,]
```

按临床信息进行样本排序，并且获得输入complexheatmap中的注释
```R
df_wide_top <- df_wide_top[, order(clinical_data["Cancer_Type", ])]
clinical_data <- clinical_data[, colnames(df_wide_top)]
ha_clinical_data <- HeatmapAnnotation(
  df = data.frame(
    Cancer_Type = as.factor(clinical_data["Cancer_Type", ]),
    Age = as.factor(clinical_data["Age", ])
  )
)
```

统计每个样本发生Multi Hit的基因数
```R
mh_count <- df_mh %>%
  filter(Variant_Classification == "Multi_Hit") %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(count=n())
```

转换Multi Hit到complexheatmap的注释
```R
mh_count <- t(mh_count)
colnames(mh_count) <- mh_count[1,]
mh_count <- mh_count[-1,]
mh_count <- as.numeric(mh_count)
barplot <- anno_barplot(mh_count, border = FALSE, labels = "Multi_Hit")
ha_multi_hit <- HeatmapAnnotation(Multi_Hit = barplot, col = list(Multi_Hit = "#8E9AAD"))

# 最后需要转换为宽格式
df_wide <- spread(df_mh, Tumor_Sample_Barcode, Variant_Classification)
df_wide <- as.data.frame(df_wide)
row.names(df_wide) <- df_wide$Hugo_Symbol
df_wide$Hugo_Symbol <- NULL
df_wide[is.na(df_wide)] <- ''
```

统计top20的基因
```R
gene_mutation_count <- rowSums(df_wide != "")
top_genes <- names(sort(gene_mutation_count, decreasing = TRUE))[1:20]
df_wide_top <- df_wide[top_genes, ]
```

绘制瀑布图
```R
# 绘图前先给每个类型配置好颜色
col = c(
  "Frame_Shift_Ins" = "#c70039",
  "Frame_Shift_Del" = "#ff5733",
  "In_Frame_Ins" = "#ff8d1a",
  "In_Frame_Del" = "#ffc300",
  "Missense_Mutation" = "#eddd53",
  "Nonsense_Mutation" = "#add45c",
  "Nonstop_Mutation" = "#57c785",
  "Start_Lost_Mutation" = "#00baad",
  "Splice_Site" = "#2a7b9b",
  "Synonymous_Mutation" = "#E6E6FA",
  "Intron_Mutation" = "#E1FFFF",
  "Other_Mutation" = "#F0F8FF",
  "Multi_Hit" = "#8E9AAD"
  )


alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
  In_Frame_Ins = alter_graphic("rect", fill = col["In_Frame_Ins"]),
  Intron_Mutation = alter_graphic("rect", fill = col["Intron_Mutation"]),
  In_Frame_Del = alter_graphic("rect", fill = col["In_Frame_Del"]),
  Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
  Synonymous_Mutation = alter_graphic("rect", fill = col["Synonymous_Mutation"]),
  Splice_Site = alter_graphic("rect", fill = col["Splice_Site"]),
  Nonsense_Mutation = alter_graphic("rect", fill = col["Nonsense_Mutation"]),
  Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),
  Other_Mutation = alter_graphic("rect", fill = col["Other_Mutation"]),
  Nonstop_Mutation = alter_graphic("rect", fill = col["Nonstop_Mutation"]),
  Start_Lost_Mutation = alter_graphic("rect", fill = col["Start_Lost_Mutation"]),
  Multi_Hit = alter_graphic("rect", fill = col["Multi_Hit"])
)
```

然后画图
```R
oncoPrint(df_wide_top, alter_fun = alter_fun, col = col, top_annotation = ha_multi_hit, bottom_annotation = ha_clinical_data)
```



效果与原始数据有关
![oncoplot](https://github.com/pzweuj/pzweuj.github.io/raw/master/downloads/images/complexheatmap_oncoplot.png)

