---
title: 亲子鉴定相关分析
tags: default
---

### 介绍

短串联重复序列（short tandem repeats，STR）也称微卫星DNA（microsatellite DNA）, 通常是基因组中由1~6个碱基单元组成的一段DNA重复序列，可以用于亲子鉴定。本人无医学、法律相关背景，下文只是学习笔记，无参考意义。

计算方式分为二联体（被检孩子+被检父亲或母亲）和三联体（被检孩子+被检父亲+被检母亲）。除了亲子关系的鉴定，还有[祖孙关系的鉴定](http://www.moj.gov.cn/news/download/file/file/20190815/1565869722167038371.pdf)。（多个被检者是双胞胎多胞胎或者近亲的这种情况是无法获得准确结果的）


### STR

《司法鉴定文书规范》中有一个由19个STR组成的检测panel。除此以外，还有[27个](https://patentimages.storage.googleapis.com/a5/f5/80/a131efdefd0ac4/CN104946632A.pdf)的，[23个](http://www.microread.com/foServen/7-893-335.html)的，[22个](https://www.jfsmonline.com/article.asp?issn=2349-5014;year=2018;volume=4;issue=3;spage=122;epage=128;aulast=Fu)的，还有ThermoFisher的[GlobalFIlter](https://assets.thermofisher.cn/TFS-Assets/GSD/brochures/globalfiler-str-brochure.pdf)、[Identifiler](https://assets.thermofisher.cn/TFS-Assets/LSG/manuals/cms_041201.pdf)、[物证鉴定中心的专利](https://patentimages.storage.googleapis.com/1d/4b/7b/a2bfaee19b154c/CN105441534B.pdf)等等等等。因为位点比较少，所以一般不会用NGS测。

选取的常染色体基因座需要满足这些条件：基因座的定义和特征需要已有文献报道，具有种属特异性、灵敏性，已有可供使用的群体遗传数据，遗传方式符合孟德尔定律，串联重复单位需要是四或者五核苷酸（不过我看到有些用三核苷酸重复的）。还有来自[CODIS系统](https://en.wikipedia.org/wiki/Combined_DNA_Index_System)（共有13个核心STR）还是非CODIS系统这些内容。

19 STR Panel，一般来说还会加入Amelogenin性别识别位点。二联体检测最少使用18个STR，但是由于STR均存在突变现象，检测数量也是越多越好（但是多了又不适用于多重荧光PCR了或者毛细管电泳了）。

D19S433、D5S818、D21S11、D18S51、D6S1043、D3S1358、D13S317 、D7S820、D16S539 、CSF1PO 、PentaD 、vWA 、D8S1179、TPOX 、PentaE、TH01、D12S391、D2S1338、FGA、Amel。

然后还需要人群频率，大概可以在这些数据库找一找，最好当然是找到中国人群的：[STRBase数据库](https://strbase.nist.gov/str_fact.htm,美国)、[STRidER,亚洲](https://strider.online/frequencies)、[Promega,亚洲](https://www.promega.com.cn/products/pm/genetic-identity/population-statistics/allele-frequencies/)。这里有[一篇文献](https://europepmc.org/article/pmc/pmc6779668)的数据，频率来自2367个南方汉族；Github上找到的[这个](https://github.com/T1me/MPTK/blob/main/frequency.yml)；另外也有用[这篇文献](https://pubmed.ncbi.nlm.nih.gov/31905040/)的，点击这里下载文献的[Supplemental](https://www.tandfonline.com/doi/suppl/10.1080/03014460.2019.1705391/suppl_file/iahb_a_1705391_sm0479.zip)。

接下来就是计算累计亲权指数CPI（Combined Parentage Index ），可以看[这篇](https://www.promega.com/-/media/files/resources/conference-proceedings/ishi-15/parentage-and-mixture-statistics-workshop/introductiontoparentagestatistics.pdf?la=en)，还有参考内容中的技术规范。



#### 亲权指数计算

根据技术规范，计算三联体亲权指数，写了一个python方法，脚本看[这里](https://github.com/pzweuj/practice/blob/master/python/PaternityIndex/calculate_CPI.py)。**未验证写法是否准确。**尽管脚本中三联体写法是疑父，已知母亲，孩子，实际上把疑母写到疑父的位置，把已知父亲写到已知母亲的位置，再把分析方式由male改为female也能成功运行（这种情况比较少见）。

还有什么排除指数之类的，未进行编写。累计非父排除率不能小于0.9999。

当累计亲权指数小于0.0001时，支持被检测男子**不是**孩子生物学父亲的假设（女子同理）；累计亲权指数大于10000时，支持被检测男子**是**孩子生物学父亲的假设（女子同理）。当累计亲权指数等于10000时，父权概率等于CPI / (CPI + 1) ≈ 99.99%



### SNP

选取中国人群多态性良好的单核苷酸位点，一次性检测5000到10000+个，可能能获得更准确的结果。但是《法医SNP分型与应用规范》只列举了几十个SNP位点（使用[微测序法](https://pubmed.ncbi.nlm.nih.gov/9199933/)检测）。SNP检测可用于个人识别和亲子鉴定。

在个人识别中，需要计算LR值（Likehood Ratio），LR等同于累积匹配概率CPM的倒数（CPM = PM1×PM2×...×PMn），统计学上LR大于1时，认为支持同一人假设，小于1则反之。当LR在数值上超过全球人口总数，表明证据有足够强度支持同一人假设，但规范未给出PM的计算公式。

一些文章中符合孟德尔定律时采取与STR相同的计算方法（具体可参阅参考内容中的技术规范），在不符合孟德尔定律时，则PI等于0.00001。另外，有基于[PMID6881203](https://pubmed.ncbi.nlm.nih.gov/6881203/)的，也是分为二联体三联体，在符合孟德尔定律时采取与STR相同算法，在不符合时PI等于0。即只要有一个点不符合孟德尔定律最终的CPI也会是0。这时还要统计PI为0的位点数，大于等于3时才认为排除亲子关系；而小于3不等于0的情况则认为无法确定亲子关系（引用该文献的一篇[专利](https://patentimages.storage.googleapis.com/fa/20/a3/263b3c9cc6ce07/CN102978286A.pdf)采用的位点数是60个）。

[华大的一篇专利](https://patentimages.storage.googleapis.com/2a/33/06/79d2af22b5539e/WO2016049878A1.pdf)选了70个位点，然后通过统计不符合孟德尔定律的位点数，当个数为0时认为是亲子关系。

[司法部研究所](https://patents.google.com/patent/CN107012226A/)设计的试剂盒有273个位点。[公安部物鉴中心](https://patents.google.com/patent/CN102115788A/)的有48个位点。




### STR+SNP
可用参考ThermoFisher的[PrecisionID](https://assets.thermofisher.com/TFS-Assets/GSD/Technical-Notes/precision_id_str_snp_combo_2019_technical_note.pdf)。当STR与SNP联合使用时，需要证据证明选取的SNP与STR之间是相互独立的，同时SNP之间也需要是独立的。计算方式未提及。



### NIPPT

Noninvasive Prenatal Paternity Testing (NIPPT)指无创**产前**亲子鉴定，与肿瘤线测血浆中游离的肿瘤DNA类似，这个项目是检测孕妇血浆中游离的胎儿DNA达到无创的效果。一般需要孕周在8周以上才能进行检测（备注：如不影响孕妇生命安全，合法人流需要在24周以下进行）。

与普通的SNP（也有基于STR的方法）亲子鉴定相比，区别主要在实验过程中，来源于孕妇的外周血需提取出白细胞和血浆，其中白细胞用于检测孕妇基因型，血浆用于检测胎儿基因型。由于血浆中游离DNA片段比较短（一般是170bp以下），引物设计一般设计为扩增120bp左右（有些文献中认为需小于100bp）。血浆结果检测出来会包含母亲和胎儿的混合基因型，因此还需要进行筛选。对于cfDNA的概率，[CN108875307A这篇专利](https://patents.google.com/patent/CN108875307A/)说得比较清楚。

![CN108875307A](https://patentimages.storage.googleapis.com/17/cb/e7/9b284360c809fc/CN108875307AD00102.png)

Github上有几个现成的项目，可以参考

[csbl-usp/NIPT](https://github.com/csbl-usp/NIPT)

[VUmcCGP/sanefalcon](https://github.com/VUmcCGP/sanefalcon)

[nshomron/hoobari](https://github.com/nshomron/hoobari)





### 参考

[Population structure of Han population in China revealed by 41 STR loci](https://www.tandfonline.com/doi/pdf/10.1080/03014460.2019.1705391)

[Informatics-based, highly accurate, noninvasive prenatal paternity testing](https://www.nature.com/articles/gim2012155.pdf)

[Noninvasive Prenatal Paternity Testing with a Combination of Well-Established SNP and STR Markers Using Massively Parallel Sequencing](https://www.mdpi.com/2073-4425/12/3/454/htm)

[Noninvasive Prenatal Paternity Testing (NIPAT) through Maternal Plasma DNA Sequencing: A Pilot Study](https://pubmed.ncbi.nlm.nih.gov/27631491/)

[Identification of sequence polymorphisms at 58 STRs and 94 iiSNPs in a Tibetan population using massively parallel sequencing](https://www.nature.com/articles/s41598-020-69137-1.pdf)

[亲权鉴定技术规范 GB-T 37223-2018](http://zjcfs.com/style/download/%E6%8A%80%E6%9C%AF%E8%A7%84%E8%8C%83/%E7%89%A9%E8%AF%81%E6%AF%92%E7%89%A9/GB-T%2037223-2018%E4%BA%B2%E6%9D%83%E9%89%B4%E5%AE%9A%E6%8A%80%E6%9C%AF%E8%A7%84%E8%8C%83.pdf)

[亲权鉴定技术规范 SF/ZJD0105001-2016](https://www.moj.gov.cn/government_service/download/file/file/20190816/1565937635270095145.pdf)

[亲子鉴定文书规范 SF/ZJD0105004-2015](http://www.moj.gov.cn/news/download/file/file/20190815/1565869707758003500.pdf)

[法庭科学DNA亲子鉴定规范 GA/T 965-2011](http://116.52.249.81/submodule/Editor/uploadfile/20170329170253616.pdf)

[法医SNP分型与应用规范 SF/ZJD0105003-2015](https://www.moj.gov.cn/news/download/file/file/20190815/1565869735018092089.pdf)

[个体识别技术规范 SF/Z JD0105012-2018](http://www.moj.gov.cn/government_public/download/file/file/20181221/1545362858363089134.pdf)

[WO2016049877A1](https://patents.google.com/patent/WO2016049877A1/)

