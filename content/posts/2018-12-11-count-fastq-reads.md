---
title: 用shell命令统计fastq的reads数
tags: coding
---

显然，fastq的明显特征有每个reads都是@开头，可以用这一点写个python脚本。
另外，也可以根据4行一个reads这一点来统计。

对于fastq文件：
```bash
echo $(cat xxx.fastq|wc -l)/4|bc
```

对于fastq.gz文件:
```bash
echo $(zcat xxx.fastq.gz|wc -l)/4|bc
```


[-_-]:前景