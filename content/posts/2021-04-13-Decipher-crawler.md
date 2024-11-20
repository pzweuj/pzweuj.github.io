---
title:  爬取DECIPHER Genomics数据库
tags: coding
---

[DECIPHER](https://www.deciphergenomics.org)数据库收集了关于拷贝数变异的已知综合征及病例信息。本次爬取主要是收集里面每个基因的pLI、LOEUF、sHet、%HI等值，其他内容并不在目标中。

通过网页的Genes页面，可以看到总共有5424个基因信息。那么通过调用selenium，先将表格显示改为100，就只有55页，这时爬下55个网页即可。


![decipher](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/decipher.png)


使用selenium调用edge浏览器，关于调用edge，可参考[此文章](https://pzweuj.github.io/2021/03/18/CHPO-crawler.html)。本来在显示100行后，打算定位到Next按钮来遍历每一页的，但是尝试了几次都有bug，在爬到第六页时定位时定位到“...”按钮去了，看了半天也不知道为啥，因此直接点击对应的目录数字按钮算了。

```python
# coding=utf-8

from msedge.selenium_tools import Edge
from msedge.selenium_tools import EdgeOptions
import time
import random

options = EdgeOptions()
options.use_chromium = True
options.add_argument("headless")
options.add_argument("--disable-blink-features")
options.add_argument("--disable-blink-features=AutomationControlled")
options.add_argument("start-maximized")
options.add_experimental_option("excludeSwitches", ["enable-automation"])
options.add_experimental_option("useAutomationExtension", False)
msedge = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedgedriver.exe"

driver = Edge(options=options, executable_path=msedge)
script = "Object.defineProperty(navigator, 'webdriver', {get: () => undefined})"
driver.execute_script(script)

url = "https://www.deciphergenomics.org/genes"

driver.get(url)
print("网址:", url)
# 等待加载
time.sleep(40)

# 定位下拉选择框并选择100
driver.find_element_by_xpath('//*[@id="content"]/div/div/div[2]/div/div/div[2]/div/div[1]/div/label/select/option[@value="100"]').click()
time.sleep(10)

# 保存第一页
n = 1
f = open("html/decipher_page" + str(n) + ".html", "wb")
f.write(driver.page_source.encode("gbk", "ignore"))
print("已完成第" + str(n) + "页")
f.close()

# 总共有55页
n += 1
while n <= 55:
    t = random.randint(5, 10)
    print("等待", t, "秒")
    time.sleep(t)
    
    # 点击每一页
    driver.find_element_by_xpath("//a[text()='{page}']".format(page=str(n))).click()
    time.sleep(5)
    f = open("html/decipher_page" + str(n) + ".html", "wb")
    f.write(driver.page_source.encode("gbk", "ignore"))
    print("已完成第" + str(n) + "页")
    f.close()
    n += 1
    
print("任务完成")
```


爬取完成后，进行解析。

```python
from bs4 import BeautifulSoup
import os

results = open("DECIPHER.txt", "w", encoding="utf-8")

n = 1
while n <= 55:
    htmlFile = "html/decipher_page" + str(n) + ".html"
    soup = BeautifulSoup(open(htmlFile), "lxml")

    table = soup.find_all("table")[0]
    tbody = table.select("tbody")[0]
    tr = tbody.find_all("tr")

    for t in tr:
        td = t.find_all("td")
        
        name_col = td[0].find_all("div")
        gene = name_col[0].string
        description = name_col[1].string
        
        loca_col = td[1].div.div.find_all("span")
        chrom = "chr" + loca_col[0].string
        start = loca_col[2].string
        end = loca_col[4].string

        pLI = td[2].span.string
        LOEUF = td[3].span.string
        sHet = td[4].span.string
        HI_precent = td[5].span.string
        
        OMIM_col = td[6]
        if len(OMIM_col.find_all("div")) == 0:
            OMIM = "-"
        else:
            if not "OMIM" in OMIM_col.div.a.string:
                OMIM = "-"
            else:
                try:
                    OMIM = OMIM_col.div.a["href"]
                except:
                    OMIM = "-"

        DDG2P_col = td[7]
        if not DDG2P_col.string == "-":
            DDG2P = DDG2P_col.a.div.string.strip()
        else:
            DDG2P = DDG2P_col.string

        ClinGen_col = td[8]
        if not ClinGen_col.string == "-":
            ClinGen_col = ClinGen_col.a.div.find_all("div")
            ClinGen_list = []
            for c in ClinGen_col:
                ClinGen_list.append(c.text.strip().replace("\n", "").replace(" ", "").replace("Tr", ";Tr"))
            ClinGen = ";".join(ClinGen_list)
        else:
            ClinGen = ClinGen_col.string

        Open = td[9].string
        
        Link_col = td[10].find_all("li")
        link_list = []
        for l in Link_col:
            link_list.append(l.a["href"])
        links = ",".join(link_list)
        
        results.write("\t".join([gene, description, chrom, start, end, pLI, LOEUF, sHet, HI_precent, OMIM, DDG2P, ClinGen, Open, links]) + "\n")
    n += 1

results.close()
```

