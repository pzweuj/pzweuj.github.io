---
title: 基于人群频率的祖源分析
tags: coding
---

## GnomAD方案

GnomAD提供预训练模型进行祖源分析，以下是原文。

[Using the gnomAD genetic ancestry principal components analysis loadings and random forest classifier on your dataset](https://gnomad.broadinstitute.org/news/2021-09-using-the-gnomad-ancestry-principal-components-analysis-loadings-and-random-forest-classifier-on-your-dataset/).

这篇文章会依据原文提供的模型进行分析，如果需要自行构建模型，建议参考原文。



## 其他方案摸索

我不打算用GnomAD的方案，而是自建一个人群分类器。安捷伦的CREV4中提供了一个包含[1300多个SNP](https://github.com/pzweuj/pzweuj.github.io/raw/refs/heads/master/content/data/document/AIM_SNPs.txt)的列表，这些SNP来自不同的研究，被认为可以作为[Ancestry Informative Markers](https://www.genome.gov/genetics-glossary/Ancestry-informative-Markers)。



### 注释

我使用SNP-nexus进行在线注释，获得不同的人群频率，我注释了HapMap。

HapMap中包含的人群有

| 人群分类 | 描述                                   |
| -------- | -------------------------------------- |
| ASW      | 美国西南部的非洲血统                   |
| CEU      | 具有北欧和西欧血统的犹他州居民         |
| CHB      | 中国北京汉族                           |
| CHD      | 美国科罗拉多州丹佛市的华人             |
| GIH      | 美国德克萨斯州休斯顿的古吉拉特印第安人 |
| HCB      | 美国科罗拉多州丹佛市的华人(没注释到内容)             |
| JPT      | 日本东京人                             |
| LWK      | Luhya，韦布耶，肯尼亚                  |
| MEX      | 美国加利福尼亚州洛杉矶的墨西哥血统     |
| MKK      | 肯尼亚基尼亚瓦的马赛人                 |
| TSI      | 意大利的托斯卡纳人                     |
| YRI      | 尼日利亚伊巴丹的约鲁巴人               |

可以看出这个分类人群的太少，缺失的地域多，不过用来可以用来看方案是否可行。

因为SNP-nexus的人群频率库不够新，也可以使用自己的注释方案进行rsid注释，比如annovar、VEP等，这样可以注释到更新的GnomAD v4.1。


### 基因型字典

根据哈迪温伯格平衡，将结果整理为一个python的字典，

```
wt_frq = (1 - p) ** 2
het_frq = 2 * p * (1 - p)
hom_frq = p ** 2
```

格式如

```python
'rs11111': {
    "ASW": {"Wt": 0.0004, "Het": 0.0096, "Hom": 0.0900},
    "CHB": {"Wt": 0.0006, "Het": 0.0104, "Hom": 0.0890}
},
'rs222222': {}
```



### 样本字典

将样本的检测结果整理为一个字典，格式如

```python
'rs11111': 'Wt',
'rs222222': 'Het'
```



### 分类代码


伪代码，需根据实际结果调参。

```python
# 计算样本所属的人群概率
def classify_population(sample_genotypes, pop_freq, epsilon=1e-10):
    # 初始化每个人群的对数概率
    population_log_probs = {pop: 0 for pop in pop_freq[list(pop_freq.keys())[0]].keys()}

    # 遍历样本的每个位点
    for rsid, genotype in sample_genotypes.items():
        if rsid in pop_freq:
            # 获取该位点的人群频率数据
            rsid_freq = pop_freq[rsid]
            for pop, probs in rsid_freq.items():
                # 获取该基因型的概率
                prob = probs.get(genotype, 0)  # 如果基因型不存在，概率为 0
                if prob == 0:
                    prob = epsilon  # 平滑处理，避免概率为 0
                population_log_probs[pop] += log(prob)
        else:
            # 如果位点不存在于 pop_freq 中，记录调试信息
            print(f"警告: 位点 {rsid} 不在 pop_freq 中")

    # 将对数概率转换为概率（避免数值下溢）
    max_log_prob = max(population_log_probs.values())  # 找到最大的对数概率
    scaled_log_probs = {pop: log_prob - max_log_prob for pop, log_prob in population_log_probs.items()}
    population_probs = {pop: np.exp(log_prob) for pop, log_prob in scaled_log_probs.items()}

    # 归一化概率
    total_prob = sum(population_probs.values())
    if total_prob > 0:  # 避免除零错误
        population_probs = {pop: prob / total_prob for pop, prob in population_probs.items()}
    else:
        # 如果所有概率都为 0，则均匀分配概率
        num_pops = len(population_probs)
        population_probs = {pop: 1.0 / num_pops for pop in population_probs}

    # 找到最可能的人群
    most_likely_population = max(population_probs, key=population_probs.get)

    return {
        "probabilities": population_probs,
        "most_likely_population": most_likely_population
    }
```
