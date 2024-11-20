---
title: DNA序列转成氨基酸序列
tags: coding
---

实际上是抄来的，然后加了点如果不满足3的倍数时候报Warning。


```python
def translate(seq):
	table = { 
		"ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
		"ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
		"AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
		"AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
		"CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
		"CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
		"CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
		"CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
		"GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
		"GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
		"GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
		"GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
		"TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
		"TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
		"TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
		"TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W",
	}
	protein = ""
	if len(seq) % 3 == 0:
		for i in range(0, len(seq), 3):
			codon = seq[i: i + 3]
			protein += table[codon]
	else:
		print ("Warning, sequence can not be translated completely")
		for i in range(0, len(seq), 3):
			if i + 3 < len(seq):
				codon = seq[i: i + 3]
				protein += table[codon]
	return protein
```



[-_-]: 123,