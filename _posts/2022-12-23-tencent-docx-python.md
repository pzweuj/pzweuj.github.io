---
title: python自动导出腾讯文档excel表格
tags: coding
---

主要思路是通过playwright这种自动测试框架去点击腾讯文档中的excel导出。需要先装好playwright。


## 保存cookies

一些文档是需要登录后才能导出的，这里可使用playwright录制登录的部分，并把cookies进行导出，然后后续使用时通过导入cookies来保持登录状态，预计一段时间会失效，因此需要一段时间就进行一次这个操作，更新cookies。

```bash
playwright codegen -o login.py https://xxx.com --save-stroage loginCookies
```

使用以上命令，进入腾讯文档进行登录，最好勾选保持登录等，最后登录状态就会保存到loginCookies文件中。

## 打开文档并导出

以下代码是打开对应excel表格并进行导出的例子，这里还是加上了[stealth.min.js脚本](https://pzweuj.github.io/2022/09/18/crawler.html)来模拟实人访问，也可以不加。这里也导入了cookies文件。如果导出的表格比较大，容易超时，因此这里将timeout设置为200秒。尚未测试无头模式是否可用。

```python
from playwright.sync_api import Playwright, sync_playwright
import time

# 获取原始表格
def run(playwright: Playwright) -> None:
    browser = playwright.chromium.launch(headless=False)
    context = browser.new_context(accept_downloads=True, storage_state="loginCookies")

    # Open new page
    page = context.new_page()
    # 模拟真人，也可以不加
    page.add_init_script(path="stealth.min.js")
    page.set_default_timeout(200000)

    # 访问网页，需要保存的表格地址
    page.goto("https://docs.qq.com/sheet/XXXXXXXXXXX")

    # 点击导出为本地excel表格
    time.sleep(10)
    page.click("[aria-label=\"file\"] >> :nth-match(div, 4)")
    time.sleep(10)
    page.click("text=导出")
    time.sleep(10)

    with page.expect_download() as download_info:
        page.click("text=本地Excel表格 (.xlsx)")
    
    # 等待表格下载完成
    time.sleep(150)
    download = download_info.value
    # 修改保存的名字，也可以不加
    download.save_as("xxxxx.xlsx")

    # Close page
    page.close()
    context.close()
    browser.close()
```



