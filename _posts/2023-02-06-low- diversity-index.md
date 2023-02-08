---
title: 低样本量下的测序barcode选择
tags: default
---

对于华大平台来说，当样本量较少时，标签（测序barcode）的低多样性会导致标签的碱基cycle碱基不平衡，测序出错，最终数据无法拆分。碱基平衡对于测序非常重要，最佳情况下，A、T、C、G碱基在文库中各占25%。对于低多样性文库，一般加入Phix文库来使文库平衡。

[illumina的文档](https://support.illumina.com/bulletins/2017/02/how-much-phix-spike-in-is-recommended-when-sequencing-low-divers0.html)说由于PhiX不含标签信息，因此无法用于平衡标签的信号，所以在测序前，对于barcode的选择也是非常重要的。

当然，由于illumina测序是检测红色(A/C)、绿色(G/T)两种[光信号的强弱](https://perkinelmer-appliedgenomics.com/2020/09/07/tech-tips-low-level-multiplexing/)（Hiseq、Miseq），因此，A/C和G/T可以认为是同一类型，只需要保持每个cycle的AC比例和GT比例在50%左右即平衡。在Miniseq、NextSeq、NovaSeq中，C'是红光、T是绿光、G不发光、A是红+绿同时发光（黄）。

当前有若干个barcode，那么从中挑选出8个barcode进行上机。我想到大概的做法是随机选择8个barcode，然后每两个之间计算汉明距离，最后求和。重复随机挑1000次（足够多次），最后输出汉明距离和最大的组合。可以补充设定每两个barcode之间的汉明距离至少为3，同时舍弃如果在同一位置中，ATCG任意碱基比例小于0.125的barcode组合。

```python
# 汉明距离
def hamming_distance(phase1, phase2):
    if len(phase1) != len(phase2):
        raise ValueError("长度不同，无法计算")
    zipped = zip(phase1, phase2)
    sum = 0
    for z in zipped:
        if z[0] != z[1]:
            sum += 1
    return sum
```

从一个barcode字典中选择若干个barcode，然后计算两两之间的汉明距离并求和

```python
# 字典格式：{"T-1": "ATCAGTGC", "T-2": "CATGCATC", "T-3": "CGATCGAT"}
def randomSelectAndCal(d, amount):
    randomSelectList = random.sample(d.keys(), amount)
    select2 = list(itertools.combinations(randomSelectList, 2))
    hmdSum = 0
    for s in select2:
        hmd = hamming_distance(d[s[0]], d[s[1]])
        hmdSum += hmd
    return [randomSelectList, hmdSum]
```


也有现成的包或软件可以用：

[Barcosel](http://ekhidna2.biocenter.helsinki.fi/barcosel/)，是一个在线的barcode选择器。

[DNABarcodes](https://bioconductor.org/packages/release/bioc/html/DNABarcodes.html)，能用于输出最佳组合的R包。





