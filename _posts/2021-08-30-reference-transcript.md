---
title: 参考转录本
tags: default
---

大部分时候，报告结果上呈现的cDNA及氨基酸突变等，都是以某个转录本为参考的。不同转录本的位置会有差异。

因此，对每个基因，固定一个参考转录本非常重要。



## Nirvana方法

[illumina的Nirvana文档](https://illumina.github.io/NirvanaDocumentation/core-functionality/canonical-transcripts)中，提到了选择参考转录本的方式。

1. 只把RefSeq中NM或NR开头的转录本作为候选；
2. 对转录本按以下优先度顺序排序：
	i. 来源于[LRG](https://www.lrg-sequence.org/)；
	II. CDS长度降序；
	III. 转录本长度降序；
	iv. 编号升序；
3. 使用排第一的结果。



## 杂七杂八

不想过于复杂，按以下方法进行。

下载LRG、MANE Select、RefSeq、Clinvar、HGNC等参考。

```bash
# LRG
wget ftp://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_xrefs.txt
# MANE Select
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.95/MANE.GRCh38.v0.95.summary.txt.gz
# RefSeq
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene
# Clinvar
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
# HGNC
wget https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt

```

然后编写代码，提取出每个数据库对于每个基因的参考转录本。优先度以LRG > MANE select > RefSeq > Clinvar > HGNC排序。

最后，对自己小panel内的核心基因进行**人工校对**，这一步不可缺少。举个例子，[MET](http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_662.xml)在LRG中的参考转录本是NM_001127500.1，而MANE select是NM_000245.4。按照上面的优先度原则，应该选择NM_001127500.1，但事实上我们会用NM_000245.4作为参考。
