---
title: 适用于VEP的突变位点侧翼注释
tags: coding
---

## 背景

现有需求是注释出突变位点的上游序列和下游序列，本来的想法是根据位点的Pos，使用samtools faidx从参考基因组中提取出对应的区域再输出。因为samtools的速度很快，效率并不低，但这样毕竟是要多进行一步。

查了一下VEP的文档，发现并没有对应的参数，也没有现成的插件（要不就是我没找到）。所以，自行编写了[一个插件](https://github.com/pzweuj/VEP_Plugins_Self/blob/main/plugins/FlankingSequence.pm)。

对于一个给定的输入，可以输出下面的格式，即发生了一个G>T突变，上游序列是GCCCATCTGTC；下游序列是TCTCTCTGATC。

```
GCCCATCTGTC[G/T]TCTCTCTGATC
```

## 下载插件

下载插件并放在自己的VEP插件目录里

```bash
wget https://github.com/pzweuj/VEP_Plugins_Self/raw/refs/heads/main/plugins/FlankingSequence.pm
```

## 使用插件

使用下面的命令进行注释

```bash
./vep -i input.vcf -o output.vcf -fa hg38.fasta --plugin FlankingSequence,10
```

结果将会被注释到```FlankingSequence```字段中。可通过**10**这个参数对侧翼序列的长度进行调整。

