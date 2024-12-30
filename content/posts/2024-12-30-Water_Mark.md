---
title: 文本文件水印
tags: coding
---


## 背景

突然的一个思考，如果售卖的数据库可以直接提供文本格式给客户下载，怎么避免客户直接把这个数据库传输给其他人？


## 方案讨论
### 方案一

法律威慑。在合同上写明如果二次分发，有权XXX，赔偿XXX。同时在文本文件抬头和结尾补充醒目的字样进一步威慑。


### 方案二

零宽度字符水印。对于每个客户，构建唯一的编号，然后将这个编号编码为一串零宽度的字符串，嵌入到这个文本文件中。


```python
ZERO_WIDTH_SPACE = '\u200B'  # 表示 0
ZERO_WIDTH_NON_JOINER = '\u200C'  # 表示 1
ZERO_WIDTH_JOINER = '\u200D'  # 分隔符
INVISIBLE_SEPARATOR = '\u2060'  # 不可见分隔符

# 将客户编号转换为不可见字符
def embed_watermark(text, client_id):
    # 将客户编号转换为二进制字符串
    binary_id = ''.join(format(ord(char), '08b') for char in client_id)
    # 将二进制字符串转换为不可见字符
    watermark = ''.join(
        ZERO_WIDTH_NON_JOINER if bit == '1' else ZERO_WIDTH_SPACE
        for bit in binary_id
    )
    # 添加不可见分隔符
    watermark = INVISIBLE_SEPARATOR.join(watermark[i:i+8] for i in range(0, len(watermark), 8))
    # 将水印嵌入文本
    return text + watermark

# 从文本中提取客户编号
def extract_watermark(text):
    # 提取所有不可见字符
    watermark = ''.join([char for char in text if char in (ZERO_WIDTH_SPACE, ZERO_WIDTH_NON_JOINER, INVISIBLE_SEPARATOR)])
    # 移除不可见分隔符
    watermark = watermark.replace(INVISIBLE_SEPARATOR, '')
    # 将不可见字符转换为二进制字符串
    binary_id = ''.join('1' if char == ZERO_WIDTH_NON_JOINER else '0' for char in watermark)
    # 将二进制字符串转换为客户编号
    client_id = ''.join(chr(int(binary_id[i:i+8], 2)) for i in range(0, len(binary_id), 8))
    return client_id


# 示例
original_text = "这是一个示例文本。"
client_id = "110101010101010101010101"

# 嵌入水印
watermarked_text = embed_watermark(original_text, client_id)
print("嵌入水印后的文本:", watermarked_text)

# 提取水印
extracted_id = extract_watermark(watermarked_text)
print("提取的客户编号:", extracted_id)
```

这是GPT给的方案，但是实际中，现在的notepad等文本编辑器都是可以显示出这些特殊字符的，如果采用这个方案，最好是在更隐蔽的地方进行嵌入。

### 方案三

虚构信息。根据自己的数据库类型进行设计。譬如说，我的数据库是基因突变信息数据库，我可以设计一批虚构的突变。只要设计大概24个位点，就可以识别10000000个客户了。

1，在热门的基因中设计，可以确保就算客户对基因进行筛选后，这些特殊设计也会保留下来；

2，和人群数据库，COSMIC数据库等比对，确保这个虚构突变不会出现在任何已知数据库中（即在真实情况中，受检者是不会检出这个突变的，避免影响）；

3，采取这批次位点的排列组合，在内部用于对客户进行识别。



### 方案四

同义词。我个人觉得这个方案最隐匿。还是以基因突变信息数据库为例。选择一批热点突变，譬如EGFR L858R、ALK G12C等，对他们的描述中的一些词语采用不影响意义的同义词。提供给每个客户的数据库的这些位点的描述都是唯一的，通过识别这批位点的特定同义词字符的排列组合，达成对客户的识别。

这个方案隐匿度高，而且不用与人群数据库等进行比较。缺陷是和上述一样，需要根据数据的特点进行设计。

伪代码

```python
# 同义词替换水印方案

# 定义热点突变及其同义词选项
hotspot_mutations = {
    "EGFR L858R": {
        "mutation": ["mutation", "variant", "alteration"],
        "associated with": ["associated with", "linked to", "related to"],
        "non-small cell lung cancer": ["non-small cell lung cancer", "NSCLC"]
    },
    "ALK G12C": {
        "mutation": ["mutation", "variant", "alteration"],
        "associated with": ["associated with", "linked to", "related to"],
        "lung cancer": ["lung cancer", "NSCLC"]
    }
}

# 生成客户唯一标识
def generate_watermark(client_id):
    # 将客户编号转换为同义词选择
    choices = []
    for bit in client_id:
        choices.append(int(bit) % 3)  # 假设每个位置有 3 个同义词选项
    return choices

# 嵌入水印
def embed_watermark(database, client_id):
    choices = generate_watermark(client_id)
    watermarked_db = database.copy()
    for i, (mutation, descriptions) in enumerate(hotspot_mutations.items()):
        for j, (key, synonyms) in enumerate(descriptions.items()):
            # 替换同义词
            watermarked_db = watermarked_db.replace(
                f"{mutation}: {key}", 
                f"{mutation}: {synonyms[choices[i * len(descriptions) + j]]}"
            )
    return watermarked_db

# 提取水印
def extract_watermark(database):
    choices = []
    for mutation, descriptions in hotspot_mutations.items():
        for key, synonyms in descriptions.items():
            # 查找使用的同义词
            for index, synonym in enumerate(synonyms):
                if f"{mutation}: {synonym}" in database:
                    choices.append(str(index))
                    break
    return ''.join(choices)

# 示例
database = """
EGFR L858R: mutation is associated with non-small cell lung cancer.
ALK G12C: mutation is associated with lung cancer.
"""

client_id = "110101010101010101010101"

# 嵌入水印
watermarked_db = embed_watermark(database, client_id)
print("嵌入水印后的数据库:\n", watermarked_db)

# 提取水印
extracted_id = extract_watermark(watermarked_db)
print("提取的客户编号:", extracted_id)
```

## 后记

实际操作中所有方案可以全用上，增强容错。
