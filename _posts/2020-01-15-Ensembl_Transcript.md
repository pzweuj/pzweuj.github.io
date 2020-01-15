---
title:  根据转录本编号，提取ensembl的cdna序列
tags: coding
---

本来想通过解析网页来做这个，但是人家ensembl都提供了下载了，还这么做就觉得太蠢了。。。

hg19下载这个文件
```bash
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz
```
里面有所有的ensembl的cdna转录本。

下载程序

```bash
wget https://github.com/pzweuj/practice/blob/master/python/Ensembl/EnsemblTS.py
```

修改好程序中hg19的路径。

然后参照说明使用即可，例如，查询ENTS00001转录本
```bash
python EnsemblTS.py -i ENTS00001
```

如果需要批量查询，则需要将素有查询的ID放置于txt文件内，每个转录本号占一行
```bash
python EnsemblTS.py -I transcriptID.txt -o output.txt
```

就可以生成需查询的转录本了。







[^_^]: 继续努力，继续挖井，继续吸花