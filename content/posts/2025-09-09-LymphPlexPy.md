---
title: 使用LymphPlexPy对DLBCL进行分型
tags: software
---

## LymphGen

弥漫大B细胞淋巴瘤（DLBCL）是成人淋巴瘤中最常见的一种类型，[LymphGen方法](https://www.sciencedirect.com/science/article/pii/S1535610820301550)将DLBCL进行细分成了MCD、N1、EZB、BN2、ST2、TP53、NOS（other）等亚型。


## LymphPlex

而[LymphPlex方法](https://www.nature.com/articles/s41392-023-01358-y)进行了优化，使用36个基因通过聚类进行亚型分型。

LymphPlex在Github中[开源了代码](https://github.com/difuSJTU/LymphPlex)，整体基于R语言编写。


## LymphPlexPy

为了进一步精简流程，优化环境变量，我对LymphPlex进行重构，代码以MIT许可开源。

[LymphPlexPy](https://github.com/pzweuj/LymphPlexPy)即是LymphPlex的Python重构版本。

### 输入文件

使用方法很简单，只需要将自己的数据构筑出一个json文件即可。基因突变/融合了就是1，没有突变/融合就是0，简单明了（当然，建议过滤掉UTR、内含子、同义突变等）。

```json
{
  "BCL2": 0,
  "BCL6": 1,
  "ARID1A": 0,
  "B2M": 0,
  "BTG1": 0,
  "BTG2": 0,
  "CCND3": 0,
  "CD70": 0,
  "CD79B": 0,
  "CIITA": 0,
  "CREBBP": 0,
  "DDX3X": 0,
  "DTX1": 0,
  "DUSP2": 0,
  "EP300": 0,
  "EZH2": 0,
  "FAS": 0,
  "GNA13": 0,
  "IRF4": 0,
  "IRF8": 0,
  "KMT2D": 0,
  "MPEG1": 0,
  "MYD88": 0,
  "NOTCH1": 0,
  "NOTCH2": 1,
  "PIM1": 0,
  "PRDM1": 0,
  "SGK1": 0,
  "SOCS1": 0,
  "STAT3": 0,
  "STAT6": 0,
  "TBL1XR1": 0,
  "TET2": 0,
  "TNFAIP3": 0,
  "TNFRSF14": 0,
  "ZFP36L1": 0
}
```

### 运行命令

运行命令如下：

```bash
python lymphplex_cli.py -i sample.json -o results.json -s PATIENT_001
```


### 输出结果

输出结果同样是一个json文件，可以方便对接自己的下游流程。

```json
 "results": [
    {
      "sample_id": "PATIENT_001",
      "predicted_subtype": "BN2",
      "classification_method": "seeds",
      "confidence_level": "high",
      "tp53_mutation": "unknown",
      "tp53_override": false,
      "myc_fusion": "unknown",
      "double_hit_lymphoma": false,
      "gene_mutations": {}
]
```

## 更多信息

更多信息和软件用法请查看[仓库页](https://github.com/pzweuj/LymphPlexPy)。

