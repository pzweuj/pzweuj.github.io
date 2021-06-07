---
title: WordWriter，自动化Word报告生成的python脚本
tags: coding
---

### 介绍

这个脚本主要靠[python-docx](https://python-docx.readthedocs.io/en/latest/)实现。同时使用[pandas](https://www.pypandas.cn/)来处理输入的表格。基于python3。

需要安装
```bash
pip install python-docx
pip install pandas
```

脚本地址：

[WordWriter](https://github.com/pzweuj/WordWriter)



### 脚本原理

使用python-docx，通过导入模板docx文件，寻找模板文件中预先保留的tag，对tag进行替换，最后重新保存docx文件达到自动输出报告的效果。支持对模板中的tag进行处理，替换后的字符串格式会完全跟随模板中对应tag的格式。仅支持docx格式的word文档，不支持doc格式word文档。



### 目前支持

整个表格插入、单元格中的字符串替换、单元格中的图片插入、文本中的字符串替换、文本框中的字符串替换、图片插入、页眉页脚字符串替换。



### 使用方法

在python中导入
```python
from WordWriter3 import WordWriter

# 建立标签与输入内容（替换内容）对照字典
testDict["#[tag]#"] = "replace string"

# 输出
WordWriter(templateDocx, outputDocx, testDict)
```



**为了程序稳定运行，插入的tag需要规范命名。**



#### 整个表格插入

表格tag命名规则：**#[TABLE-xxxx]#**

表格的输入为以tab分割的txt文件路径，即

```python
testDict[#[TABLE-xxxx]#] = "/path/to/table/table.txt"
```



需要插入表格，首先需要在模板中定义一个与插入表格列数一致的表格，并且将tag放置在需要填充的第一行第一列（不包含标题）。如下表。**填充时，文字格式会按照tag所在行每个单元格的格式进行填充（即每一列的格式都会参考tag所在行的单元格格式）。**

| 标题1       | 标题2 | 标题3 |
| ----------- | ----- | ----- |
| #[TABLE-1]# |       |       |
|             |       |       |
|             |       |       |



当然，如果第一列是固定的字段，也可以将tag放置在第一行第二列中（类推）。

| 标题1 | 标题2       | 标题3 |
| ----- | ----------- | ----- |
| 一    | #[TABLE-2]# |       |
| 二    |             |       |
| 三    |             |       |

注意，应在模板中预先建立表格，并新建若干行，如模板表格行数小于输入表格行数，程序会自行创建新行；如模板表格行数大于输入表格行数，程序会自动删除空行。



#### 单元格中的字符串替换

单元格字符串tag命名规则：**#[TBS-xxxx]#**

输入为字符串，即

```python
testDict["#[TBS-xxxx]#"] = "table cell string replace"
```

多用于以下情况：

| 标题1 | 标题2           |
| ----- | --------------- |
| 结果  | #[TBS-results]# |



#### 单元格中的图片插入

单元格图片tag命名规则：**#[TBIMG-xxxx]#**

包含大小格式：**#[TBIMG-xxxx-(20,20)]#**

输入为图片文件路径

```python
testDict["#[TBIMG-1-(20,20)]#"] = "/path/to/picture/picture.png"
```

所有图片均建议加入大小格式，不要使用原始尺寸进行插入。括号中分别是(width,height)。

| #[TBIMG-1-(20,20)]# | #[TBIMG-2-(20,20)]# | #[TBIMG-3-(20,20)]# |
| ------------------- | ------------------- | ------------------- |



#### 文本中的字符串替换

文本中字符串tag命名要求不要与其他内容的tag前缀冲突即可。规则：**#[xxxxx]#**

输入为字符串

```python
testDict["#[testString]#"] = "string replace"
```

替换的内容的格式会完全跟随模板中tag的格式，即在模板中对**tag**进行加粗，替换内容也会加粗；将tag颜色更改为<font color=red>红色</font>，替换内容也会是红色。单元格字符串替换以及文本框字符串替换同理。



#### 文本框中的字符串替换

文本框中的字符串tag命名规则：**#[TX-xxxx]#**

输入为字符串

```python
testDict["#[TX-textbox]#"] = "textbox string replace"
```



#### 图片插入

图片插入tag命名规则：**#[IMAGE-xxxx]#**

包含大小格式：**#[IMAGE-xxxx-(20,20)]#**

输入为图片文件路径

```python
testDict["#[IMAGE-2-(20,20)]#"] = "/path/to/picture/picture.png"
```



#### 页眉页脚字符串替换

页眉tag命名规则：**#[HEADER-xxxx]#**

页脚tag命名规则：**#[FOOTER-xxxx]#**

输入为字符串

```python
testDict["#[HEADER-date]#"] = "2021年6月7日"
testDict["#[FOOTER-project]#"] = "test footer"
```



### 测试

[测试文件](https://github.com/pzweuj/WordWriter/tree/master/test)

测试脚本
```python
from WordWriter3 import WordWriter

# 测试脚本
testDict = {}
testDict["#[HEADER-1]#"] = "模板测试"
testDict["#[HEADER-2]#"] = "2021年6月7日"
testDict["#[NAME]#"] = "测试模板"
testDict["#[fullParagraph]#"] = "这是一段测试段落，通过WordWriter输入。"
testDict["#[TBS-1]#"] = "未突变"
testDict["#[TX-1]#"] = "文本框测试成功"
testDict["#[TX-2]#"] = "文本框测试很成功"
testDict["#[FOOTER]#"] = "页脚测试"

# 此处输入的是文件路径
testDict["#[TABLE-1]#"] = "test/testTable.txt"
testDict["#[IMAGE-1-(30,30)]#"] = "test/testPicture.png"
testDict["#[IMAGE-2]#"] = "test/testPicture.png"
testDict["#[TBIMG-3-(20,20)]#"] = "test/testPicture.png"

# 使用主函数进行报告填充
WordWriter("test/test.docx", "test/testOut.docx", testDict)
```


### 注意事项

有时tag并不能被很好的识别，是因为word将一段不连续的输入的tag理解成为是多个“run”。程序需要一个tag作为一个run时才能识别出来。如果转换docx格式为xml可以看到tag的内容被分在了不同的标签中。

所以，如果遇到这种情况，建议是将不能识别的tag完整的复制到文本文档中，再完整的（一次性的）粘贴回模板中替换掉不识别的tag。在修改格式如颜色字体等的时候也需要是将一个tag完全选中来修改，避免被认为是多次的输入内容。


