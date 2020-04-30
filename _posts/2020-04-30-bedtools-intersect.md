---
title: 才发现bedtools intersect -v是这样的
tags: default
---

使用bedtools intersect -v 找B.bed未覆盖A.bed的位置。输出结果里面，理论上A的总长度应该等于B的总长度加上新生成的C的总长度，由于结果总是不对，去官方文档看了一下，发现原来-v参数不能做到我想要的操作。

![intersect](https://bedtools.readthedocs.io/en/latest/_images/intersect-glyph.png)

```bash
bedtools intersect -a A.bed -b B.bed -v > C.bed
```
得到的结果，实质上是整个区域都没有overlap时，才会输出的。

那么应该怎么才能获得所有没有覆盖的位置呢，可以使用bedtools subtract。
![subtract](https://bedtools.readthedocs.io/en/latest/_images/subtract-glyph.png)

```bash
bedtools subtract -a A.bed -b B.bed > C.bed
```

这才是我想要的。