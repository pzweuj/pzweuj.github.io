---
title: ncbi的genotyping tool
tags: software
---

之前做HBV的分型的时候，发现ncbi有这个分型的工具，[genotyping tool](https://www.ncbi.nlm.nih.gov/projects/genotyping/help.html)。原理是使用blast将query序列与设定好的blast库进行blast获得最相似的序列。

目前可以做分型的病毒有HIV-1，HIV-2，HTLV-1，HTLV-2，HBV，HCV，PV。

使用起来也比较简单，例如，随便复制一段HBV序列，
```
TTTTTCTTGTTGACAAGAATCCTCACAATACCGCAGAGTCTAGACTCGTGGTGGACTTCTCTCAATTTTC
TAGGGGGAACTACCGTGTGTCTTGGCCAAAATTCGCAGTCCCCAACCTCCAATCACTCACCAACCTCTTG
TCCTCCAACTTGTCCTGGTTATCGCTGGATGTGTCTGCGGCGTTTTATCATCTTCCTCTTCATCCTGCTG
CTATGCCTCATCTTCTTGTTGGTTCTTCTGGACTATCAAGGTATGTTGCCCGTTTGTCCTCTAATTCCAG
GATCCTCAACAACCAGCACGGGACCATGCCGGACCTGCATGACTACTGCTCAAGGAACCTCTATGTATCC
CTCCTGTTGCTGTACCAAACCTTCGGACGGAAATTGCACCTGTATTCCCATCCCATCATCCTGGGCTTTC
GGAAAATTCCTATGGGAGTGGGCCTCAGCCCGTTTCTCCTGGCTCAGTTTACTAGTGCCATTTGTTCAGT
GGTTCGTAGGGCTTTCCCCCACTGTTTGGCTTTCAGTTATATGGATGATGTGGTATTGGGGGCCAAGTCT
GTACAGCATCTTGAGTCCCTTTTTACCGCTGTTACCAATTTTCTTTTGTCTTTGGGTATACATTTAAACC
CTAACAAAACAAAGAGATGGGGTTACTCTCTAAATTTTATGGGTTATGTCATTGGATGTTATGGGTCCTT
GCCACAAGAACACATCATACAAAAAATCAAAGAATGTTTTAGAAAACTTCCTATTAACAGGCCTATTGAT
TGGAAAGTATGTCAACGAATTGTGGGTCTTTTGGGTTTTGCTGCCCCTTTTACACAATGTGGTTATCCTG
CGTTGATGCCTTTGTATGCATGTATTCAATCTAAGCAGGCTTTCACTTTCTCGCCAACTTACAAGGCCTT
TCTGTGTAAACAATACCTGAACCTTTACCCCGTTGCCCGGCAACGGCCAGGTCTGTGCCAAGTGTTTGCT
GACGCAACCCCCACTGGCTGGGGCTTGGTCATGGGCCATCAGCGCATGCGTGGAACCTTTTCGGCTCCTC
TGCCGATCCATACTGCGGAACTCCTAGCCGCTTGTTTTGCTCGCAGCAGGTCTGGAGCAAACATTATCGG
GACTGATAACTCTGTTGTCCTATCCCGCAAATATACATCGTTTCCATGGCTGCTAGGCTGTGCTGCCAAC
TGGATCCTGCGCGGGACGTCCTTTGTTTACGTCCCGTCGGCGCTGAATCCTGCGGACGACCCTTCTCGGG
```
上面这是一段[NC_003977.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_003977.2?report=fasta)。将这段序列放入[genotyping](https://www.ncbi.nlm.nih.gov/projects/genotyping/formpage.cgi)中，选择HBV分型，点击Subtype，得到结果，表面这个作为HBV参考基因组的序列其实是HBV D型。
![genotyping](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/downloads/images/genotyping.jpg)

由图其实可以知道，对HBV分型的方式是对序列划分窗口然后每一段拿去和参考的已经分型样本进行blast。


[.](https://www.zhihu.com/question/274328777)


[-_-]:jing

