---
title: 使用Page-Agent爬取OncoKB（夯爆了）
tags: default
---

从Qwen开源起，阿里一直走在国产AI的第一梯队。[page-agent](https://github.com/alibaba/page-agent)是纯JS实现的GUI agent，使用自然语言操作Web应用，无须后端、客户端或浏览器插件。

page-agent以MIT开源，意味着商业产品可用。



### OncoKB爬取

请注意page-agent的核心逻辑是**操作Web应用**，而不是“爬虫”。我尝试让page-agent帮我整理OncoKB的信息。

这一次和之前的OpenClaw、Manus等逻辑不同，不是让它自己打开网页去整理数据，而是一句一句叫它去做。

首先，我打开了 `https://www.oncokb.org/cancer-genes`，然后让page-agent进入ABL1基因的页面，并给了下面这个指令：

```prompt
接下来你需要帮我收集Annotated Alterations里所有的信息，形成一个文件
```

它顺利地把信息整理了下来。

![page-agent](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/page_agent_1.png)



#### 不足

当然，还是存在不足。我用的是浏览器插件方式安装，它整理的结果只放在了对话框里，并没有按照我的要求整理成文件。

不过，page-agent原生是一个JS脚本，完全可以二次开发或建立插件，让它把特定内容保存下来并整理成特定格式的文件。



### 展望

如上所述，page-agent的核心逻辑是**操作Web应用**。我目前在开发的分析平台完全可以接入page-agent，建立知识库和特定的SKILL，让它自行判断位点，并通过**点击**的方式去决策是否回报。这意味着，我们的分析平台拥有了一个全自动的**解读工程师**。

