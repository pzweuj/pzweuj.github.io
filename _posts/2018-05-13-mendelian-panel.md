---
title: 弄了个孟德尔遗传病的panel
tags: database
---

>其实就是几个数据库的使用而已啦。


说起孟德尔遗传病（遵循孟德尔遗传定律的基因病），首先想到的数据库肯定是[OMIM](https://www.omim.org)。
![omim](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/OMIM.PNG)
我们需要的是申请下载OMIM的数据。申请地址点[这里](https://omim.org/downloads/)。
主要需要的是genemap2.txt这个文件。

因为做的是panel，所以一定要有据可循，要做到李菊福。
所以，我们可以去用免费的疾病数据库[Clinvar](https://www.ncbi.nlm.nih.gov/clinvar/)。
![clinvar](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/Clinvar.PNG)
可以在Clinvar中筛选出致病的位点，但是这样不能保证都是孟德尔遗传病，所以，要增加一个筛选条件，就是提交者为OMIM。

我发现了一个不错的网站，利用Clinvar的数据进一步归类分析。叫做[Clinvar Miner](https://clinvarminer.genetics.utah.edu)。
![clinminer](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/Clinvar_Miner.PNG)

所以我直接在Clinvar Miner中选择了OMIM作为提交者的位点。
然后再选择致病位点。大概是24000多个。
直接下载完整列表，由于我只是想做一个简单的panel，所以只挑选了其中有rsid的部分，总共有19000多个。

接下来，利用annovar的注释rsid功能，把这19000多个位点的位置信息注释出来，并且挑选出东亚人突变频率小于5%的位点（罕见病）。

然后，再利用新版的Clinvar注释，把其中提示Benign的位点都剔除掉，最后剩下8600多个位点。

再然后，利用genemap2这个文件，把OMIM的ID和疾病名称注释上去。这时已经是一个可用的panel了。
然后可以利用[CHPO](http://www.chinahpo.org)数据库，找到疾病的中文名字。
![chpo](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/CHPO.PNG)
CHPO的数据也是可以下载的，不过申请起来很麻烦。可以尝试用爬虫抓取。

[T_T]:思路清晰工作效率高