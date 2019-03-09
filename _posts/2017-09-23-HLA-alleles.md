---
title: HLA分型与资料数据库
tags: database
---
> HLA是一个非常难的项目，查了无数资料和网站，发现了这个，资料很丰富。
[HLA-alleles](http://hla.alleles.org/nomenclature/index.html)。

[HLA](https://en.wikipedia.org/wiki/Human_leukocyte_antigen)
系统是目前所知人体最复杂的多态系统。
自1958年发现(Jean Dausset)第一个HLA抗原，到20世纪70年代，
HLA便成为免疫遗传学、免疫生物学和生物化学等学科的一个重要新兴研究领域。
现在，已基本弄清其系统的组成、结构和功能，阐明了其理化性质和生物学作用。
这些研究成果不仅具有重要的理论意义，而且具有巨大的生物医学价值。

在生信的领域，主要就是通过确定HLA基因的分型，来通过分型设定对应的治疗方案。
比如说HLA-B\*58:01，HLA-B\*27。这是比较常见的相关区域。

HLA-alleles这个网页提供了每个HLA区域的[分型表](http://hla.alleles.org/alleles/text_index.html)。

比如说[HLA-B区](https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/B_nuc.txt)
![hla-b-alignments](https://github.com/pzweuj/pzweuj.github.io/raw/master/downloads/images/HLA-B-alignments-nuc.png)
可以看到，这里用了B\*07:02:01:01这个分型作为参照组，下面的都是对照组，
其中，“—”表示的是缺失，“\*”表示的相同。然后和参照组不一样的就会直接的列出AGTC。

然后，HLA怎么确定分型呢。首先，第一步是设计相应的引物。比如，要想确定B区的分型，首先要设计能扩展B区的引物，然后通过Sanger测序，得到下机的数据。在输出测得的序列，通过和HLA-B区的这个分型表进行比对（BLAST），就可以得到分型结果。

另外，也有比较省时间的方式。
第一步是要自己弄出这个分型表，可以去EMBL下载HLA的序列信息。点击[这里](ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/)。

通过muscle等比对方法，就可以比对出差异。然后，将差异的位点以及相应的位置找出来。然后在Sanger测序之后，只要去查看关键的位置（可以利用R的SangerSeqR这个包实现自动化），就可以得到分型。这个方法有一个刚性的要求，就是必须，引物一定要好！

-------------------------------------------
还有一件事，网页左侧有一个Lead sponsors。他家的HLA分型试剂盒是最领先的最准确的。再配合他家的非卖品分型软件，效果非常好。不说了，涉嫌打广告了已经。


-------------------------------------------

说到HLA分型软件，还有一下免费的作品，比如[Seq2HLA](https://bitbucket.org/sebastian_boegel/seq2hla/src)。但是这些对没什么编程知识的人不太友好。

PS. 要搜索HLA相关的分型的临床表现，可以去之前介绍过的[SNPedia](https://www.snpedia.com/index.php/SNPedia)等网站。


[^_^]:真的不骗你，的确是刚好打开刚好看到。