---
title: 使用Firecrawl爬取网页
tags: software
---

[Firecrawl](https://www.firecrawl.dev/)是一个可以将网站爬取为一个LLM适用的格式（markdown）的工具，特别适用于静态网站。

因此，我的思路是先用Firecrawl将网页爬取为markdown格式，再用通义千问将我需要的内容提取出来，并整理结果成一个json或txt，最后整理成内部数据库的格式。

大概就是将playwright等获得网页数据的过程，转为了用Firecrawl获取，并且转换成markdown；然后将用Beautifulsoup提取网页数据的过程，转换为了用LLM来提取（当然也可以在爬取过程中就输入prompt让Firecrawl直接将结果输出为自己想要的格式，但这需要不断调试，因此我还是习惯先保存源码到本地）。

Firecrawl目前貌似[尚不支持爬取需要认证的网页](https://github.com/mendableai/firecrawl/issues/1093)，这篇issue说需要保存cookies再进行处理。

## Firecrawl

### 官方API

我们可以申请一个Firecrawl账户，这样就可以获得API KEY并且能调用Firecrawl的API。免费用户可以爬取500个网页以及拥有50万tokens。

```python
from firecrawl import FirecrawlApp
app = FirecrawlApp(api_key="fc-YOUR_API_KEY")
scrape_result = app.scrape_url('firecrawl.dev', formats=['markdown', 'html'])
data = scrape_result.markdown
```

上述代码可以将单个页面转换为markdown格式，另外Firecrawl也提供了一个爬取子页面的方案。

```python
from firecrawl import FirecrawlApp, ScrapeOptions
app = FirecrawlApp(api_key="fc-YOUR_API_KEY")
# Crawl a website:
crawl_result = app.crawl_url(
  'https://firecrawl.dev', 
  limit=10, 
  scrape_options=ScrapeOptions(formats=['markdown', 'html']),
)
print(crawl_result)
```

可以使用prompt来在爬取过程中就只提取感兴趣的信息，但在这里不赘述。

### 自部署

[Firecrawl以AGPLv3开源](https://github.com/mendableai/firecrawl)，当500个网页不能满足需求时，可以进行本地部署（也可以成为付费用户）。关于自部署，可以搜索其他教程。


### 伪代码

以下代码，未经测试。

```python
from firecrawl import FirecrawlApp

def fire_base_use(url, api_key):
    app = FirecrawlApp(api_key=api_key)
    scrape_result = app.scrape_url(url, formats=['markdown'])
    with open(os.path.join('output.md"), "w", encoding='utf-8') as o:
        o.write(data.markdown)
```

## 结果整理

使用Qwen3整理并提取需求的信息。系统提示词就自行编写咯。当然，也可以使用AI为你编写的。


```python
from openai import OpenAI

def get_item_sub_url(sub_name, input_file, output_file, api_key):
    model = "qwen3-30b-a3b"
    api_base_url = "https://dashscope.aliyuncs.com/compatible-mode/v1"

    client = OpenAI(
        api_key = api_key,
        base_url = api_base_url,
    )

    system_prompt = f"""
    从以下提供的Markdown文本中，请你执行以下操作：
    1.  **定位目标链接：** 找到所有URL中包含 `/xxxxx/detail/` 关键字的链接。这些链接主要出现在{sub_name}部分的具体子项中。
    2.  **提取标题和链接：**
        *   对于每个定位到的链接，提取其完整的URL。
        *   提取该链接的**显示文本**（即Markdown链接中方括号 `[]` 里的内容）。
        *   从这个显示文本中，**只截取第一个`\\`（代表换行符）之前的部分作为最终的标题**。
    3.  **输出格式：** 将提取到的标题和链接以以下格式呈现，每个一行，不包含任何额外的描述或文本：
        `[标题]\t[URL]`
    **请确保输出中只包含精确的标题（不含风险描述）和对应的链接，并且严格遵循上述输出格式。**
    """
    
    # 读取输入文件内容
    with open(input_file, 'r', encoding='utf-8') as f:
        markdown_content = f.read()
    
    completion = client.chat.completions.create(
        model = model,
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": markdown_content}
        ],
        stream = True
    )

    # 收集流式响应
    collected_content = ""
    for chunk in completion:
        if hasattr(chunk, 'choices') and len(chunk.choices) > 0:
            if hasattr(chunk.choices[0], 'delta') and hasattr(chunk.choices[0].delta, 'content'):
                if chunk.choices[0].delta.content is not None:
                    collected_content += chunk.choices[0].delta.content

    # 将结果写入输出文件
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(collected_content)
    return collected_content
```
