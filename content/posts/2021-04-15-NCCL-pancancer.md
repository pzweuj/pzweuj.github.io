---
title:  怎么规范写突变
tags: default
---



精准与规范同行，每次的NCCL室间质评，又到了卫健委教你写突变结果的时间。



### SNV与INDEL

| Chr  | Start    | End      | Ref  | Alt  | Gene | Type | Transcript  | cHGVS     | pHGVS       | VAF%  | Consequence           | Affected_Exon |
| ---- | -------- | -------- | ---- | ---- | ---- | ---- | ----------- | --------- | ----------- | ----- | --------------------- | ------------- |
| 7    | 55249071 | 55249071 | C    | T    | EGFR | SNV  | NM_005228.5 | c.2369C>T | p.Thr790Met | 23.12 | Missense_substitution | 19/27         |



SNV类型中，start与end均要是参考基因组上的位置，两个值必须一致。Transcript要使用[Clinvar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)最新版的转录本，**必须**包含版本号，即例子中的**.5**。cDNA突变与氨基酸突变均需要采用[HGVS](https://varnomen.hgvs.org/recommendations/general/)的规则，目前的几大注释软件中，annovar需要使用annotate_variation.pl加上-hgvs参数，而snpEff、VEP、Nirvana等均默认HGVS规则。而基因名则需要参考[HGNC](https://www.genenames.org/)的基因名。由于[Clinvar的文献](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753237/)里写了，Clinvar的转录本是参考RefSeqGene以及[LRG](ftp://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_xrefs.txt)的，因此我们可以从LRG、RefSeq以及HGNC中获得参考转录本的编号，形成一个转录本与基因的对照字典（排序分先后，LRG优先）。对于pHGVS应该还是可以缩写如p.T790M的。Affected Exon则是在当前转录本中的外显子排位以及外显子总数。

而对于SNV来说，Consequence有Synonymous_substitution（同义突变）、Missense_substitution（错义突变）、Nonsense_substitution（无义突变）。



| Chr  | Start    | End      | Ref                | Alt  | Gene | Type     | Transcript  | cHGVS            | pHGVS                     | VAF%  | Consequence      | Affected_Exon |
| ---- | -------- | -------- | ------------------ | ---- | ---- | -------- | ----------- | ---------------- | ------------------------- | ----- | ---------------- | ------------- |
| 7    | 55242467 | 55242484 | AATTAAGAGAAGCAACAT | -    | EGFR | Deletion | NM_005228.5 | c.2237_2254del18 | p.Glu746_Ser75 2delinsAla | 23.12 | Inframe_deletion | 19/27         |

对Deletion缺失类型来说，start与end必须完全与ref序列一致（使用"samtools faidx hg19.fa chr7:55242467-55242484"可获得AATTAAGAGAAGCAACAT），而alt为“-”。Deletion的Consequence有Inframe_deletion（框内删除）以及Frameshift_deletion（移码删除）。




| Chr  | Start    | End  | Ref  | Alt  | Gene | Type      | Transcript  | cHGVS             | pHGVS                | VAF%  | Consequence       | Affected_Exon |
| ---- | -------- | ---- | ---- | ---- | ---- | --------- | ----------- | ----------------- | -------------------- | ----- | ----------------- | ------------- |
| 7    | 55249010 | -    | -    | GTT  | EGFR | Insertion | NM_005228.5 | c.2308_2309insGTT | p.Asp770delinsGlyTyr | 23.12 | Inframe_insertion | 19/27         |

对Insertion插入类型来说，start是插入位置坐标，end与ref均是“-”。Insertion的Consequence有Inframe_insertion（框内插入）以及Frameshift_insertion（移码插入）。




| Chr  | Start    | End      | Ref  | Alt  | Gene | Type    | Transcript  | cHGVS             | pHGVS       | VAF%  | Consequence          | Affected_Exon |
| ---- | -------- | -------- | ---- | ---- | ---- | ------- | ----------- | ----------------- | ----------- | ----- | -------------------- | ------------- |
| 9    | 36923458 | 36923459 | GG   | AT   | PAX5 | Complex | NM_016734.3 | c.803_804delinsTA | p.Ala268Asp | 23.12 | **Complex_mutation** | 7/10          |

对于以上这种MNV来说，不能把点拆成多个SNV来写，类型是Complex。【另外，Consequence是Complex_mutation（复杂突变）】。以上表格是NCCL原文原例，实际上，这个例子是也和NCCL下文解析中描述相悖。**原文大意：对于仅影响单个密码子的替换，只能选Nonsense_substitution、Missense_substitution、Synonymous_substitution三种；而在影响多密码子时才选Complex_mutation**。因此上述例子其实应该是一个Missense_substitution。



突变类型除SNV、Insertion、Deletion外的，均写为Complex。Consequence除以上描述的情况以及Splice_Site_mutation（剪接点改变）外，记为Other。



另外，正链基因（以正链为模板转录，出现在重复序列上的突变根据转录本上的 3’rule 发生在 3’端，在基因组层面，也按照转录本 上 3’端书写，即突变发生在重复序列的**右侧**）；负链基因（以负链为模板转录，出现在重复序列上的突变根据转录本上的 3’rule，发生在 3’端，在基因组层面，也按照转录 本上 3’端书写，即突变发生在重复序列的**左侧**）。



### SV

SV分为四种，倒位（Inversion）、易位(Translocation)、大片段插入/缺失（Insertion/Deletion）。其中大片段的插入和缺失是指大于 1kb 的片段。可以使用lumpy等软件进行检测。

| Chr  | Breakpoint1 | Breakpoint2 | Type      | VAF%  | Gene1 | Gene2 | Annotation      |
| ---- | ----------- | ----------- | --------- | ----- | ----- | ----- | --------------- |
| 2    | 29432591    | 42530406    | Inversion | 10.25 | ALK   | EML4  | EML4-ALK Fusion |

同一染色体例子如上表，即Inversion。Insertion 或 Deletion结果也按上表，但Gene1，Gene2相同，Annotation不写即可。




| Chr  | Gene1 | Breakpoint1 | Type          | Chr  | Gene2 | Breakpoint2 | VAF%  | Annotation     |
| ---- | ----- | ----------- | ------------- | ---- | ----- | ----------- | ----- | -------------- |
| 2    | ALK   | 29432591    | Translocation | 3    | TFG   | 100455436   | 10.25 | ALK-TFG Fusion |

不同染色体例子如上表，即Translocation。

其实在我看来还不如都按照表2来写，Chr写同一个就好了，统一风格。



### CNV

（这个在本次NCCL里没有，我脑补的）CNV就比较简单了，只有CNV gain和CNV loss两种。只需要写以下内容即可，建议使用cnvkit检测。

| Chr  | Gene | Region   | Type     | VAF    |
| ---- | ---- | -------- | -------- | ------ |
| 13   | FLT3 | Exon1-24 | CNV gain | 4 copy |


