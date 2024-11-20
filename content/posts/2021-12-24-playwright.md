---
title: 尝试用playwright写爬虫
tags: software
---

[playwright](https://github.com/microsoft/playwright-python)是微软开源的一个自动化测试Chromium、Firefox和WebKit的python工具，很明显，这种工具往往都会被用来做爬虫。



首先需要安装playwright

```bash
pip install --upgrade pip
pip install playwright
playwright install
```

以上步骤会按照python的playwright库，然后按照Chromium、Firefox和WebKit等的测试工具。

[官方文档](https://playwright.dev/python/docs/intro)中介绍了很多使用方式，这里还有一篇介绍playwright、selenium等工具差异的[文章](https://www.testim.io/blog/puppeteer-selenium-playwright-cypress-how-to-choose/)。



以爬取[DECIPHER gene](https://www.deciphergenomics.org/genes)这个页面为例，使用下面命令打开新窗口，默认是使用chromium

```bash
playwright codegen -o test.py
```



在浏览器中输入网址：https://www.deciphergenomics.org/genes，等待加载完成后，在页面数下拉框中选择100，可以看到此时共有页数56页。点击第2页。这时结束，查看自动化生成的代码。

```python
from playwright.sync_api import Playwright, sync_playwright


def run(playwright: Playwright) -> None:
    browser = playwright.chromium.launch(headless=False)
    context = browser.new_context()

    # Open new page
    page = context.new_page()

    # Go to https://www.deciphergenomics.org/genes
    page.goto("https://www.deciphergenomics.org/genes")

    # Select 100
    page.select_option("[aria-label=\"Number of rows per page\"]", "100")

    # Click [aria-label="Page 2 of 56"]
    page.click("[aria-label=\"Page 2 of 56\"]")

    # Close page
    page.close()

    # ---------------------
    context.close()
    browser.close()


with sync_playwright() as playwright:
    run(playwright)
```



以上代码由playwright自动录制生成，为了爬取全部56页，需要进行一些修改，见下面代码的中文注释

```python
from playwright.sync_api import Playwright, sync_playwright
import time

def run(playwright: Playwright) -> None:
    browser = playwright.chromium.launch(headless=False)
    context = browser.new_context()

    # Open new page
    page = context.new_page()

    # Go to https://www.deciphergenomics.org/genes
    page.goto("https://www.deciphergenomics.org/genes")
    # 首次打开网页等待加载
    time.sleep(40)

    # Select 100
    page.select_option("[aria-label=\"Number of rows per page\"]", "100")
    # 获得网页并保存
    html = page.content()
    with open("html/page1.html", "w", encoding="utf-8") as f:
        f.write(html)
    
    # Click [aria-label="Page 2 of 56"]
    # 调整为56页的循环
    n = 2
    while n <= 56:
        page.click("[aria-label=\"Page {} of 56\"]".format(str(n)))
        time.sleep(10)
        with open("html/page{}.html".format(str(n)), "w", encoding="utf-8") as f:
            f.write(html)
        n += 1

    # Close page
    page.close()

    # ---------------------
    context.close()
    browser.close()


with sync_playwright() as playwright:
    run(playwright)
```



使用录制功能就能非常快的编写好爬取代码。如果还需要指定UA等信息，可以参考playwright的文档，也能简单查看

```bash
playwright codegen --help
```

