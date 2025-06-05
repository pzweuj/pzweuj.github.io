---
title: 使用OncoSense快速提取指南靶药信息
tags: coding
---

二月底时整理了[靶向药的名称](https://pzweuj.github.io/posts/Cancer_Drug_Name)，当时就有一个想法：我现在有了靶药的名称，是不是就能按图索骥在指南里找到对应的用药信息了，利用LLM来达成量子速读。

## 介绍

[OncoSense](https://github.com/pzweuj/OncoSense)这个程序采用的是策略是：

1. 将靶药名称库转换为向量数据库；
2. 将输入的指南/文献向量化，与靶药向量库比较相似度，获得可能被提及的靶药，并将对应的块作为上下文；
3. 使用LLM进行信息提取，获得这个指南/文献中对应的靶药信息。

OncoSense目前存在以下问题：

1. 依赖于药物数据库，无法识别新药，需要及时更新数据库；
2. 依赖于向量化的相似度阈值，当阈值较高或较低，匹配就可能错误，当前默认设定了一个测试相对稳定的阈值；
3. 当更换了向量化模型，切记需要重新建立靶药向量库。


OncoSense使用3个模型，来自通义千问全家桶：

1. text-embedding-v3：用于向量化（现在最新版本是v4，请注意如果更新向量模型，必须重置数据库）；
2. qwen2.5-vl-7b-instruct：用于OCR；
3. qwen-max：用于提取并整理信息。

你也可以修改config为自己喜欢的模型。


## 安装

拉取github仓库，装上依赖库

```bash
git clone https://github.com/pzweuj/OncoSense.git
pip install -r requirements.txt
```

在config.py里配置自己的阿里云百炼API KEY

```python
API_CONFIG = {
    "dashscope_api_key": "your_api_key_here",
    # ...
}
```


## 使用

在默认阈值和参数下，直接使用

```bash
python main.py input.pdf -o output.json
```


## 结果示例

OncoSense输出一个json格式的结果，你可以根据自己的需求二次提取并整理信息。

```json
{
  "extracted_drugs": [
    {
      "drug_name": "吉非替尼",
      "english_name": "Gefitinib",
      "targets": ["EGFR"],
      "sensitivity_info": [
        {
          "mutation": "EGFR L858R",
          "sensitivity": "敏感",
          "evidence_level": "1A"
        }
      ],
      "page_numbers": [3, 7],
      "confidence_score": 0.85
    }
  ]
}
```

实际测试中，分析NCCN非小细胞指南 (2025.3)，共计285页，耗时约9分钟，测试结果文件[点击此处](https://github.com/pzweuj/OncoSense/blob/main/test/%5BNCCN%5D%5BNon-Small%20Cell%20Lung%20Cancer%5D%5B3.2025%5D.json)。


## 叠甲

项目仅为研究使用，结果均由AI生成，需要人工仔细甄别。

