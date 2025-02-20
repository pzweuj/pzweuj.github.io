---
title: 批量下载NCCN指南
tags: coding
---

使用playwright批量下载NCCN英文指南。~~后续我想对指南建立向量知识库，然后使用DeepSeek来批量整理其中的靶向用药位点信息。~~注意，该行为违反NCCN的最终用户许可，千万不要做。下面也只是一个测试代码，我也不知道有没有用🤪。

使用下面的代码前，首先需要注册NCCN的账户。

照例，为了反反爬，会用到(stealth.min.js)[https://github.com/requireCool/stealth.min.js]

## 爬取代码

下面是爬取的python代码，修改自己的账户密码。为了避免失败，分两阶段进行，第一阶段只查询pdf的网址并保存，然后在第二阶段再进行下载。如果已获得第一阶段文件，完全可以直接进行第二阶段。


```python
from playwright.sync_api import Playwright, sync_playwright
import os
import time

# 创建保存PDF的目录
if not os.path.exists("pdfs"):
    os.makedirs("pdfs")

# 创建日志文件
log_file = "results.txt"

def collect_pdf_links(playwright: Playwright) -> None:
    """第一阶段：收集PDF链接并保存到results.txt"""
    # 启动浏览器
    browser = playwright.chromium.launch(headless=False)
    context = browser.new_context()
    page = context.new_page()
    page.add_init_script(path="stealth.min.js")
    
    try:
        # 跳转到登录页面
        print("跳转到登录页面...")
        page.goto("https://www.nccn.org/login")
        time.sleep(2)
        
        # 输入账户和密码
        print("输入账户和密码...")
        page.locator("#Username").fill("[账户]")
        time.sleep(5)  # 等待登录完成
        page.get_by_label("Password").fill("[密码]")
        page.get_by_role("button", name="Log in").click()
        time.sleep(5)  # 等待登录完成
        
        # 打开Guidelines页面
        print("打开Guidelines页面...")
        page.goto("https://www.nccn.org/guidelines/category_1")
        time.sleep(5)
        
        # 获取Treatment by Cancer Type下的癌症类型
        print("获取癌症类型信息...")
        cancer_links = page.query_selector_all('a[href*="guidelines-detail"]')
        cancer_info = []
        for link in cancer_links:
            name = link.inner_text().strip()
            href = link.get_attribute("href")
            cancer_info.append({"name": name, "url": href})
        
        print(f"共找到 {len(cancer_info)} 种癌症类型：")
        for info in cancer_info:
            print(f"{info['name']} - {info['url']}")
        
        # 清空日志文件并写入表头
        with open(log_file, "w", encoding="utf-8") as f:
            f.write("癌种 | 版本 | PDF链接\n")
            f.write("-" * 50 + "\n")
        
        # 遍历每个癌症类型，打开详情页并记录PDF链接和版本号
        for cancer in cancer_info:
            print(f"正在处理癌症类型: {cancer['name']} ({cancer['url']})")
            
            # 打开癌种详情页
            page.goto(("https://www.nccn.org" + cancer["url"]).replace("nccn-guidelines/", ""))
            time.sleep(5)
            
            # 查找NCCN Guidelines链接及其旁边的版本号
            guidelines_link = page.get_by_role("link", name="NCCN Guidelines", exact=True)
            if guidelines_link.is_visible():
                pdf_url = guidelines_link.get_attribute("href")
                print(f"找到PDF链接: {pdf_url}")
                
                # 获取版本号
                try:
                    version_element = guidelines_link.locator("xpath=following-sibling::span")
                    version_text = version_element.inner_text().strip()
                    version = version_text.replace("Version", "").strip()  # 提取版本号
                    print(f"找到版本信息: {version}")
                except Exception as e:
                    print(f"未找到版本信息: {e}")
                    version = "Unknown"
                
                # 将结果写入日志文件
                with open(log_file, "a", encoding="utf-8") as log:
                    log.write(f"{cancer['name']} | {version} | {pdf_url}\n")
            else:
                print(f"未找到NCCN Guidelines链接: {cancer['name']}")
                # 将结果写入日志文件
                with open(log_file, "a", encoding="utf-8") as log:
                    log.write(f"{cancer['name']} | Unknown | 未找到PDF链接\n")
        
        print("所有PDF链接已收集完成！")
    
    finally:
        # 关闭浏览器
        page.close()
        context.close()
        browser.close()

def download_pdfs(playwright: Playwright) -> None:
    """第二阶段：读取results.txt并下载PDF"""
    # 启动浏览器
    browser = playwright.chromium.launch(headless=False)
    context = browser.new_context()
    page = context.new_page()
    page.add_init_script(path="stealth.min.js")
    
    try:
        # 跳转到登录页面
        print("跳转到登录页面...")
        page.goto("https://www.nccn.org/login")
        time.sleep(2)
        
        # 输入账户和密码
        print("输入账户和密码...")
        page.locator("#Username").fill("[账户]")
        time.sleep(5)  # 等待登录完成
        page.get_by_label("Password").fill("[密码]")
        page.get_by_role("button", name="Log in").click()
        time.sleep(5)  # 等待登录完成
        
        # 读取results.txt文件
        with open(log_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
        
        # 跳过表头
        lines = lines[2:]
        
        for line in lines:
            cancer_name, version, pdf_url = line.strip().split(" | ")
            pdf_url = pdf_url.strip()
            if pdf_url == "未找到PDF链接":
                print(f"跳过 {cancer_name}，未找到PDF链接")
                continue
            
            print(f"正在下载 {cancer_name} 的PDF文件...")
            
            # 构造文件名
            file_name = f"[NCCN][{cancer_name}][{version}].pdf"
            file_path = os.path.join("pdfs", file_name)
            
            # 使用Playwright下载PDF
            with page.expect_download() as download_info:
                page.goto("https://www.nccn.org" + pdf_url)  # 访问PDF链接
            download = download_info.value
            download.save_as(file_path)
            print(f"PDF下载成功: {file_path}")
    
    finally:
        # 关闭浏览器
        page.close()
        context.close()
        browser.close()

with sync_playwright() as playwright:
    # 第一阶段：收集PDF链接 || 如已完成，注释掉第一阶段
    collect_pdf_links(playwright)
    
    # 第二阶段：下载PDF
    download_pdfs(playwright)
```

