---
title: 使用Qwen 2.5 VL 转换pdf为markdown
tags: coding
---

[阿里云百炼平台](https://bailian.console.aliyun.com/)现在注册免费送1M tokens，先来薅一下这个羊毛。

我拿自己的pdf文件截图，分别测试了一下72B、7B和3B三个模型。其中3B模型错误的将“谷”识别了为“合”，而72B和7B表现都优秀。

如果进行本地自部署的话，其实我觉得3B都够用。这里我为了平衡费用与准确度，决定使用7B，将pdf转换为一个markdown文件。

方案是先使用[PyMuPDF](https://pymupdf.readthedocs.io/)将pdf的每页分割为图片，然后使用7B对图片进行解析，转换为markdown格式的文字描述。具体效果尚可，可以对接进我的业务。

还有下面这些方案，可以考虑使用：

- 巨硬的[markitdown](https://github.com/microsoft/markitdown)
- [MinerU](https://github.com/opendatalab/MinerU)

## 代码

我写了[一个脚本](https://github.com/pzweuj/practice/tree/master/python/pdf2markdown)来执行pdf转markdown的过程，我转换了一个220页的pdf文件，大概消耗了50000个tokens。Qwen 2.5 VL 7B在百炼平台的费用大概是输入 2元/1M tokens，输出 5元/1M tokens，这样算下来如果没有免费额度的话，转换这个220页的文件大概要花费2毛钱。当然，模型是开源的，硬件性能充足可以自部署。


### 依赖

项目需要安装这些python库

```
openai
PyMuPDF
```

### 使用

使用下面的命令执行

```bash
python pdf2markdown_core.py -i input.pdf -o output.md -t 一级标题 -k 百炼API_KEY
```

### GUI版本

我使用[streamlit](https://streamlit.io/)写了一个GUI。使用命令如下

```bash
python -m streamlit run pdf2markdown_gui.py
```

