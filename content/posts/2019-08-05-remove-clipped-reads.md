---
title: 去除bam文件中的clipped reads
tags: coding
---

bam文件里面会存在soft clipped 和 hard clipped等结果，有时我们只需要完美匹配结果，可以对bam文件进行处理。（当然，也可以在比对时就设定相关参数，比如bowtie2的end-to-end）

去除soft clipped 和 hard clipped 的方法来自[biostar](https://www.biostars.org/p/137461/)。

```bash
samtools view sample.bam | awk '$6 ~ /H|S/{print $1}' | sort -k1,1 | uniq > sample.names.txt
samtools view sample.bam | sort -k1,1 > sample.tmp.sam
samtools view -H sample.bam > sample.new.sam
join -t '	' -v 1 -1 1 -2 1 sample.tmp.sam sample.names.txt >> sample.new.sam
samtools view -bS sample.new.sam > sample.noclip.bam
```

注意，join -t 命令那里是一个tab。


[-_-]:jing

