---
title: WordWriter，一个用来填docx模板的模块
tags: coding
---


###  说明
这个模块的主要目的是用于填充docx模板，原理是在模板中保留对应的tag，然后通过模块去把tag替换掉。理论上，替换内容的格式会完全跟随tag的格式，因此格式的调整只需要在模板中进行。

模块基于[pandas](https://pandas.pydata.org/)和[python-docx](https://python-docx.readthedocs.io/en/latest/index.html)。

需要安装
```bash
pip install pandas
pip install python-docx
```

模块下载：
[python2](https://raw.githubusercontent.com/pzweuj/WordWriter/master/WordWriter.py)

[python3](https://raw.githubusercontent.com/pzweuj/WordWriter/master/WordWriter3.py)


### 模板说明
只支持docx格式的word文档，不支持doc格式word文档。不支持文本框及图形中的文本替换。

#### 表格的插入
需要插入表格，首先需要在模板中定义一个列数一致的表格，并且将tag放置在需要填充的第一行第一列（不包含标题），当然，如果第一列是固定的字段，也可以将tag放置在第一行第二列中（类推）。填充时，格式会按照tag所在行的格式进行填充。表格的tag仅支持#[TABLE-xxxx]#的格式。

例1：
![ww1](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/wordw_1.jpg)

上面的表格中，tag定义为#[TABLE-1]#，后续填充则会从该行开始，以后每一行的格式都会参照tag所在行的格式。表格并不需要设置很多行，只需要设置格式行（标题行除外）即可，程序会自动根据插入的行数来补充表格行数。另外，如果需要每行背景色不一样，目前的建议是在模板中创建一个超行数的表格，并且在模板中调整奇偶行的背景色。程序会在填充后自动删除空白行。

例2：
![ww2](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/wordw_2.jpg)

如上图，只填充后面四列。

#### 单独单元格文字的插入
表格中文字的替换与表格的插入不同，固定的tag仅支持#[TBS-xxxx]#的样式。
例3：
![ww3](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/wordw_3.jpg)


#### 图片的插入
支持插入图片。需要注意的是，目前仅支持图片单独作为一个段落，不支持在文字中插入图片。
图片的tag仅支持#[IMAGE-xxx]#的样式。并且，当需要定义图片插入的长度和高度时，可以定义一个这样的tag：#[IMAGE-xxx-(30,40)]#，即表示为，插入图片的长度为30，高度为40。

#### 页眉页脚的插入
页眉的tag仅支持#[HEADER-xxx]#的样式，页脚的tag仅支持#[FOOTER-xxxx]#的样式，支持不同节中定义不同的页眉和页脚。

#### 段落文本的插入
段落文本支持自定义tag，但是最好还是设置为#[xxxx]#的样式，同时，不要使用以上保留的tag样式。注意段落文本和表格中的文本是不一样的。对段落文本tag设置什么样的格式，替换的内容即会是什么样的格式。
![ww4](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/wordw_4.jpg)

如上图，在段落中设置了一个tag，最终替换的样式也会跟随该tag的样式。


### 模块使用说明
设置好模板，再来说明如何使用模块。
首先需要导入模块，以下为python代码。
```python
import WordWriter
```
新建一个字典，将tag作为字典的key，将需要输入的内容作为字典的value，注意，表格的输入内容是一个**文件**，该文件是以tab分割的文本文件，每一行即为表格的行，每个tab分割的内容即为一个格子的内容。该文本文件的列数需要和填充范围的列数一致，不然会报错。图片的输入内容也是一个**图片文件**。其他的输入内容均是字符串。

demo如下：
```python
testDict = {}
testDict["#[HEADER-1]#"] = "模板测试"
testDict["#[HEADER-2]#"] = "2019年7月18日"
testDict["#[NAME]#"] = "测试模板"
testDict["#[fullParagraph]#"] = "这是一段测试段落，通过WordWriter输入。"
testDict["#[TBS-1]#"] = "表格内容测试"
testDict["#[FOOTER]#"] = "页脚测试"

# 此处输入的是文件路径
testDict["#[TABLE-1]#"] = "testTable.txt"
testDict["#[IMAGE-1-(30,30)]#"] = "testPicture.png"
testDict["#[IMAGE-2]#"] = "testPicture.png"
```

在这之前需要做的，就是将需要输入的结果处理成中间文件，或者处理成字符串。然后放入字典中，最后使用WordWriter进行结果输出。
```python
WordWriter.WordWriter("template.docx", "output.docx", testDict)
```
### 注意
有时tag并不能被很好的识别。可能是因为word将一段不连续的输入的tag理解成为是多个“run”。程序需要一个tag作为一个run时才能识别出来。所以，如果遇到这种情况，建议是将不能识别的tag完整的复制到文本文档中，再完整的（一次性的）粘贴回模板中替换掉不识别的tag。在修改格式如颜色字体等的时候也需要是将一个tag完全选中来修改，避免被认为是多个run。


[^_^]: 准备好了