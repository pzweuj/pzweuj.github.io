---
title: 扩增子流程以及复现文章（2）
tags: default
---
这一篇是画图啦。

读入数据
---
其中sample_info是记录了每个样本的状态之类信息的文件。
```R
library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")

countdata <- read.table("DADA2/ASVs_counts.txt", header=T, row.names=1, check.names=F)
taxdata <- as.matrix(read.table("DADA2/ASVs_taxonomy.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

sample_info <- read.table("Rawdata/sample_info.txt", header=T, row.names=1, check.names=F)
```

前处理
---
```R
## 把样本分一下类，之前说过了是有空白对照的
# 把对照组和实验组的总count数统计出来，方便之后做差异分析
blank_ASV_counts <- rowSums(countdata[,1:4])
sample_ASV_counts <- rowSums(countdata[,5:20])

# 由于实验组是16个样本，对照组是4个样本，所以把实验组的结果除以4来标准化
norm_sample_ASV_counts <- sample_ASV_counts/4

# 找出大概的差异基因
blank_ASVs <- names(blank_ASV_counts[blank_ASV_counts * 10 > norm_sample_ASV_counts])
length(blank_ASVs)

# 去掉空白之后剩下的
colSums(countdata[!rownames(countdata) %in% blank_ASVs, ]) / colSums(countdata) * 100

# 利用刚刚的结果过滤一次
filt_count_tab <- countdata[!rownames(countdata) %in% blank_ASVs, -c(1:4)]
# 过滤之后剩下的样本的信息
filt_sample_info_tab<-sample_info[-c(1:4), ]

# 设置作图颜色
filt_sample_info_tab$color[filt_sample_info_tab$char == "water"] <- "blue"
filt_sample_info_tab$color[filt_sample_info_tab$char == "biofilm"] <- "darkgreen"
filt_sample_info_tab$color[filt_sample_info_tab$char == "altered"] <- "chocolate4"
filt_sample_info_tab$color[filt_sample_info_tab$char == "glassy"] <- "black"
filt_sample_info_tab$color[filt_sample_info_tab$char == "carbonate"] <- "darkkhaki"

filt_sample_info_tab
```

Alpha分析
---
```R
# 稀释曲线
rarecurve(t(filt_count_tab), step=100, col=filt_sample_info_tab$color, lwd=2, ylab="# of ASVs", xlab="# of Sequences")
abline(v=(min(rowSums(t(filt_count_tab)))))

# 富集还有丰度分析
# 创建新的phyloseq对象
count_tab_phy <- otu_table(filt_count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(taxdata)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

# 画出Chao1和Shannon
# 横坐标为样本
plot_richness(ASV_physeq, color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$char)])) +
  theme(legend.title = element_blank())
# 横坐标为状态
plot_richness(ASV_physeq, x="type", color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$char)])) +
  theme(legend.title = element_blank())
```

Beta分析
---
```R
# 使用DESeq2进行差异分析
deseq_counts <- DESeqDataSetFromMatrix(filt_count_tab, colData = filt_sample_info_tab, design = ~type)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# 生成表格
vst_trans_count_tab <- assay(deseq_counts_vst)

# 计算出差异矩阵
euc_dist <- dist(t(vst_trans_count_tab))

# 聚类
euc_clust <- hclust(euc_dist, method="ward.D2")
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- filt_sample_info_tab$color[order.dendrogram(euc_dend)]
labels_colors(euc_dend) <- dend_cols

# 画出树
plot(euc_dend, ylab="VST Euc. dist.")

# 使用差异矩阵
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(filt_sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# 使用phyloseq画PCoA图
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues

# PCoA图
plot_ordination(vst_physeq, vst_pcoa, color="char") + 
  labs(col="type") + geom_point(size=1) + 
  geom_text(aes(label=rownames(filt_sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$char)])) + 
  theme(legend.position="none")
```


参考文献：[Example marker-gene workflow](https://astrobiomike.github.io/amplicon/workflow_ex)

[-_-]:dd