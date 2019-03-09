---
title: GO/KEGG 富集分析
tags: coding
---

>可以说是自己慢慢琢磨着写的第一个R程序了。。主要用的是大佬写的clusterProfiler这个包。


首先就是装这些包。
```R
#source("https://bioconductor.org/biocLite.R")
#biocLite("DOSE") 
#biocLite("topGO")
#biocLite("clusterProfiler")
#biocLite("pathview")
```


载入包。
```R
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
```

导入数据。原始数据是一个基因列表。只有一列就是基因，有标题。
```R
data <- read.table("gene.list",header=TRUE)
data$GeneName <- as.character(data$GeneName)
```

利用org.Hs.eg.db转换基因名。
```R
transID = bitr(data$GeneName,
	fromType="SYMBOL",
	toType=c("ENSEMBL", "ENTREZID"),
	OrgDb="org.Hs.eg.db"
)
```

建立文件夹来装结果。
```R
dir.create("GO")
dir.create("KEGG")
```

GO_CC注释。
```R
CC <- enrichGO(transID$ENTREZID,
	"org.Hs.eg.db",
	keyType="ENTREZID",
	ont="CC",
	pvalueCutoff=0.05,
	pAdjustMethod="BH",
	qvalueCutoff=0.1
)
CC <- setReadable(CC, OrgDb=org.Hs.eg.db)

svg(file="./GO/GO_CC_Dotplot.svg", bg="transparent")
dotplot(CC, showCategory=12, colorBy="pvalue", font.size=8, title="GO_CC") # + theme(axis.text.y = element_text(angle = 45))
dev.off()

svg(file="./GO/GO_CC_Barplot.svg", bg="transparent")
barplot(CC, showCategory=12, title="GO_CC", font.size=8)
dev.off()

svg(file="./GO/GO_CC_Network.svg", bg="transparent")
plotGOgraph(CC)
dev.off()

write.table(as.data.frame(CC@result), file="./GO/GO_CC.xls", sep="\t", row.names=F)
```

GO_MF注释。
```R
MF <- enrichGO(transID$ENTREZID, "org.Hs.eg.db", keyType="ENTREZID", ont="MF", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
MF <- setReadable(MF, OrgDb=org.Hs.eg.db)

svg(file="./GO/GO_MF_Dotplot.svg", bg="transparent")
dotplot(MF, showCategory=12, colorBy="pvalue", font.size=8, title="GO_MF") # + theme(axis.text.y = element_text(angle = 45))
dev.off()

svg(file="./GO/GO_MF_Barplot.svg", bg="transparent")
barplot(MF, showCategory=12, title="GO_MF", font.size=8)
dev.off()

svg(file="./GO/GO_MF_Network.svg", bg="transparent")
plotGOgraph(MF)
dev.off()

write.table(as.data.frame(MF@result), file="./GO/GO_MF.xls", sep="\t", row.names=F)

```

GO_BP注释。
```R
BP <- enrichGO(transID$ENTREZID, "org.Hs.eg.db", keyType="ENTREZID", ont="BP", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
BP <- setReadable(BP, OrgDb=org.Hs.eg.db)

svg(file="./GO/GO_BP_Dotplot.svg", bg="transparent")
dotplot(BP, showCategory=12, colorBy="pvalue", font.size=8, title="GO_BP") # + theme(axis.text.y = element_text(angle = 45))
dev.off()

svg(file="./GO/GO_BP_Barplot.svg", bg="transparent")
barplot(BP, showCategory=12, title="GO_BP", font.size=8)
dev.off()

svg(file="./GO/GO_BP_Network.svg", bg="transparent")
plotGOgraph(BP)
dev.off()

write.table(as.data.frame(BP@result), file="./GO/GO_BP.xls", sep="\t", row.names=F)
```

KEGG注释。
```R
kegg <- enrichKEGG(transID$ENTREZID, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
kegg <- setReadable(kegg, OrgDb=org.Hs.eg.db, keytype="ENTREZID")

svg(file="./KEGG/KEGG_Dotplot.svg", bg="transparent")
dotplot(kegg, showCategory=12, colorBy="pvalue", font.size=8, title="KEGG") # + theme(axis.text.y = element_text(angle = 45))
dev.off()

svg(file="./KEGG/KEGG_Barplot.svg", bg="transparent")
barplot(kegg, showCategory=12, title="KEGG", font.size=8)
dev.off()

write.table(as.data.frame(kegg@result), file="./KEGG/kegg.xls", sep="\t", row.names=F)

dir.create("./KEGG/MAP")
kegg_df = as.data.frame(kegg)

for(i in kegg_df$ID){
  pathview(gene.data=transID$ENTREZID,
           pathway.id=i,
           species="hsa",
           kegg.native=TRUE,
           kegg.dir="./KEGG/MAP"
  )
}

print("TASK DONE")
```

还有很多奇奇怪怪的地方，需要后续的维护和修改。

[T_T]:老婆＃我想你做我老婆的。