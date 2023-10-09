---
title: WordWriter，python自动化Word报告
tags: coding
---

### 介绍

这个脚本主要靠[python-docx](https://python-docx.readthedocs.io/en/latest/)实现。同时使用[pandas](https://www.pypandas.cn/)来处理输入的表格。基于python3。使用复杂度没有[docxtpl](https://docxtpl.readthedocs.io/en/latest/)那么高。解决了之前版本标签跨run后无法识别的痛点。

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

整个表格插入、单元格中的字符串替换、单元格中的图片插入、文本中的字符串替换、文本框中的字符串替换、图片插入、页眉页脚字符串替换、表格列合并等。


### 使用方法

通过识别**#[xxx]#**这种格式的tag，进行内容替换。

#### 基础使用方法
在python中导入
```python
from WordWriter import WordWriter

# 预先建立的word模板
templateDocx = "template.docx"
outputDocx = "output.docx"

# 建立标签与输入内容（替换内容）对照字典
testDict = {}
testDict["#[tag]#"] = "replace string"

# 输出
WordWriter(templateDocx, outputDocx, testDict)
```



#### 整个表格插入

表格tag命名规则：**#[TABLE-xxxx]#**

表格的输入为**以tab分割的文本文件**路径，即

```python
testDict["#[TABLE-xxxx]#"] = "/path/to/table/table.txt"
```

由于txt文件中换行符会识别为下一行，当需要在同一单元格进行换行时，将文本文件中对应的换行符替换为**\\x0a**。

需要插入表格，首先需要在模板中定义一个列数足够的表格，并且将tag放置在需要填充的第一行第一列（不包含标题）。如下表。**填充时，文字格式会按照tag所在行单元格的格式进行填充。**即第一列会按照#[TABLE-1]#的格式，第二列按照格式1，第三列按照格式2。

| 标题1       | 标题2 | 标题3 |
| ----------- | ----- | ----- |
| #[TABLE-1]# | 格式1 | <font color=red>格式2</font> |
|             |       |       |
|             |       |       |



当然，如果第一列是固定的字段，也可以将tag放置在第一行第二列中（类推）。

|      | 标题2       | 标题3 |
| ---- | ----------- | ----- |
| 一   | #[TABLE-2]# | 格式1 |
| 二   |             |       |
| 三   |             |       |

如模板表格行数小于输入表格行数，程序会自行创建新行；如模板表格行数大于输入表格行数，程序会自动删除空行。

当将对应的value设置为\#DELETETHISTABLE#时，会将整个表格对象删除

```python
testDict["#[TABLE-xxxx]#"] = "#DELETETHISTABLE#"
```



#### 图片插入

单元格图片tag命名规则：**#[IMAGE-xxxx]#**

可以包含大小格式：**#[IMAGE-xxxx-(20,20)]#**

输入为图片文件路径

```python
testDict["#[TBIMG-1-(20,20)]#"] = "/path/to/picture/picture.png"
```

所有图片均建议加入大小格式，不要使用原始尺寸进行插入。括号中分别是(width,height)。

| 图1                 | 图2                 | 图3                 |
| ------------------- | ------------------- | ------------------- |
| #[IMAGE-1-(20,20)]# | #[IMAGE-2-(20,20)]# | #[IMAGE-3-(20,20)]# |

当插入的图片不存在时，图片插入模式会更改为文字插入模式，插入图片的路径。因此，可以通过下面的方式，达到让图片插入区为空的效果

```python
testDict["#[TBIMG-1-(20,20)]#"] = ""
```



#### 文本框中的字符串替换

文本框中的字符串tag命名规则：**#[TX-xxxx]#**

输入为字符串，但是需要注意，目前仅支持替换掉文本框中所有的内容。

```python
testDict["#[TX-textbox]#"] = "textbox string replace"
```

**图形**中的文本也等同于文本框中文本。




#### 段落/页眉/页脚/单元格中的字符串

段落/页眉/页脚/单元格中的字符串tag命名规则：**#[xxxx]#**，只要不要和IMAGE、TABLE、TX等保留前缀冲突即可。

输入为字符串，即

```python
testDict["#[xxxx]#"] = "table cell string replace"
```

替换的内容的格式会完全跟随模板中tag的格式，即在模板中对tag进行加粗，如**#[xxxx]#**，替换内容也会加粗；将tag颜色更改为红色，如<font color=red>#[xxxx]#</font>，替换内容也会是红色。单元格字符串替换以及文本框字符串替换同理。


如单元格中的替换：

| 标题1 | 标题2           |
| ----- | --------------- |
| 结果  | #[results]# |

当需要设置空段落时，一般建议将替换成空内容即可，但由于这会保留段落的占位，因此特殊情况下可以将对应的value设置为#DELETETHISPARAGRAPH#来删除整个段落对象

```python
# 建议使用
testDict["#[xxxx]#"] = ""
# 特殊情况
testDict["#[xxxx]#"] = "#DELETETHISPARAGRAPH#"
```

### 更多的格式调整
更多的格式调整，就只能通过python-docx自行打补丁了。建议先使用替换的方式，先输出一个中间文件，再根据自身需求调整格式。

```python
from WordWriter import WordWriter

# 先输为中间文件
WordWriter(input, tmp, resultDict)

# 自定义调整格式的函数
def patch(input, output):
    pass
patch(tmp, output)
```

#### 表格行合并
很多时候，需要合并表格同一列中内容相同的行，需要将word文件重新导入为document对象，设定需要合并表格及对应的列索引。

```python
from docx import Document
from WordWriter import MergeTableRow

input = "tmp.docx"
d = Document(input)
table = d.tables[0]
col_index = 1
MergeTableRow(table, col_index)
```

### 测试

[测试文件](https://github.com/pzweuj/WordWriter/tree/master/test)

测试脚本，将WordWriter.py也放在同一文件夹时
```python
from WordWriter import WordWriter

# 测试脚本
resultsDict = {}
resultsDict["#[testheader1]#"] = "测试页眉1"
resultsDict["#[testheader2]#"] = "页眉测试2"
resultsDict["#[testString]#"] = "，文本替换成功"
resultsDict["#[testfooter]#"] = "测试页脚"
resultsDict["#[TX-testString2]#"] = "，文本框文本替换成功"
resultsDict["#[testTableString1]#"] = "单元格文本替换成功"
resultsDict["#[testTableString2]#"] = "单元格文本替换成功"
resultsDict["#[IMAGE-test1-(30,30)]#"] = "testPicture.png"
resultsDict["#[IMAGE-test2]#"] = "testPicture2.png"
resultsDict["#[IMAGE-test3-(10,10)]#"] = "testPicture.png"
resultsDict["#[TABLE-test1]#"] = "testTable.txt"

# 使用主函数进行报告填充
WordWriter("test.docx", "output.docx", resultsDict)
```
