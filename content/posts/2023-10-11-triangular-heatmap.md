---
title: ggplot2绘制三角热图
tags: coding
---

目标是使用ggplot2绘制出类似maftools中somaticInteractions()模块产生的三角热图。
maftools产出如下：

![maftools\_th](https://github.com/pzweuj/pzweuj.github.io/raw/master/content/data/images/maftools_th.png)


maftools这副三角热图用于展示共突变基因的相关度。首先构建出一个文本文件，格式类似下面这个表格，这里使用的是fisher精确检验，绘制热图时，主要使用的是P值，但同时会根据OR值是否大于1来判断是正相关还是负相关。

```
x	y	value	OR
PMS2	TP53	0.0036828676433501074	0.2
PMS2	TET2	0.4199168731045784	0.6547619047619048
PMS2	SH2B3	0.5570301431667071	1.7105263157894737
PMS2	SETD2	0.5895552875452755	2.4922279792746114
PMS2	SETBP1	1.0	1.2195121951219512
PMS2	PMS2	NA	1
SETBP1	TP53	1.0	0.8695652173913043
SETBP1	TET2	0.5199504039165167	2.074074074074074
SETBP1	SH2B3	0.7342468809923415	1.7142857142857142
SETBP1	SETD2	NA	1
SETBP1	SETBP1	NA	1
SETD2	TP53	1.0	0.5467980295566502
SETD2	TET2	NA	1
SETD2	SH2B3	NA	1
SETD2	SETD2	NA	1
SH2B3	TP53	0.4397247870116061	0.6043956043956044
SH2B3	TET2	0.7364614903174054	1.8313253012048192
SH2B3	SH2B3	NA	1
TET2	TP53	0.23820858847966137	0.4230769230769231
TET2	TET2	NA	1
TP53	TP53	NA	1
```



注意，因为是三角热图，只画半边，最好在处理数据时，绘图前，就对基因进行排序。下面是ggplot2绘图脚本。

```R
library(ggplot2)

# 读入数据
data <- read.csv("test.txt", sep="\t", header=T)
data_log <- data
data_log$value <- log10(data_log$value) * -1
data_log <- as.data.frame(data_log)
data_log$value <- ifelse(data_log$OR < 1, data_log$value * -1, data_log$value)

# 绘图
ggplot(data = data_log, aes(x = x, y = y)) +
  geom_tile(aes(fill = value)) +
  geom_text(data = subset(data, value < 0.05), aes(label = "*"), color = "black") +
  geom_text(data = subset(data, value < 0.1 & value >= 0.05), aes(label = "·"), color = "black") +
  coord_equal(clip = "off") +
  xlab("") + ylab("") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.2)) +
  scale_fill_gradientn(colors = c("#8B6914", "#F5F5F5", "#00868B"),
                       na.value = "#FFFFFF",
                       limits = c(-3, 3),
                       oob = scales::squish,
                       labels = c(">3(Mutually exclusive)", "2", "1", "0", "1", "2", ">3(Co-occurence)"),
                       name = "-log10(p_value)") +
  scale_x_discrete(position = "top") +
  geom_text(x = Inf, y = 5, label = "* P < 0.05", hjust = 0.2, vjust = 1, size = 3, color = "black", show.legend = FALSE) +
  geom_text(x = Inf, y = 4, label = "· P < 0.1", hjust = 0.2, vjust = 1, size = 3, color = "black", show.legend = FALSE)
```

可能需要微调一下最后加上去的两个P值图层的y值。

![ggplot2\_th](https://github.com/pzweuj/pzweuj.github.io/raw/master/content/data/images/ggplot2_th.png)
