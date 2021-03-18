---
title: 爬取CHPO数据库
tags: coding
---

CHPO即[china HPO](http://www.chinahpo.org/)，是在中文人类表型标准用语联盟倡导下建立的一个公共网站，希望提供一个共享的平台有助于研究人员和医学专家共同翻译编辑Human Phenotype Ontology，以形成一个中文版的HPO。

通过查看网页源代码，发现想要的中文翻译等信息不在源码中，因此用常规的requests可能不行，因此使用[selenium](https://pypi.org/project/selenium/)通过模拟浏览器进行。由于我用的是新版edge浏览器，因此需要先到微软找找文档。微软自己有一篇介绍[selenium如何调用edge的文章](https://docs.microsoft.com/zh-cn/microsoft-edge/webdriver-chromium/?tabs=c-sharp)，需要安装指定版本的python包。

```python
pip install msedge-selenium-tools selenium==3.141
```

同时下载与自己目前使用的edge版本相同的[webdriver](https://developer.microsoft.com/zh-cn/microsoft-edge/tools/webdriver/)。


到chpo网站看了看，发现HPO编号从HP:0000002到HP:3000079都有，如果按照5s爬一个网页算，也得大半年才能爬完（单线程）。后来发现，并不是每个HP编号都有相关信息的，但是在chpo网站又不容易获得HP编号，因此就到[HPO网站](https://hpo.jax.org/app/)去找了。

HPO网站里能很容易获得这个genes_to_phenotype.txt文件，通过这个文件可以获得HP编号。

```bash
cut genes_to_phenotype.txt -f3 | uniq > hpolist.txt
```

最终其实只有23万左右的HP编号有用。下面是爬取代码，后期可以使用Thread和queue等库来增加线程数，获得更快速度。但是爬网站还是不能太过分，因此还是只用单线程了。这里爬取过程中顺便解析了，也可以爬完再解析。

```python
from selenium import webdriver
from selenium.webdriver import Edge
from bs4 import BeautifulSoup
import os
import time
import random


def chinahpo(hpo, dic):
    time.sleep(random.randint(5, 10))
    driver = Edge(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedgedriver.exe")
    hpid = hpo.split(":")[1]
    url = "http://www.chinahpo.org/#/searchList?trigger=1&tabType=1&searchContent=HP%3A{hpid}".format(hpid=hpid)

    try:
        driver.get(url)
        strtemp = url
        print("网址：", strtemp)
    except Exception:
        print("get page error", hpo)

    time.sleep(2)
    with open("html/hp_" + hpid + ".html", "a+", encoding="utf-8") as f:
        f.write(str(driver.page_source))

    driver.close()

    dic[hpo] = {
        "url": "",
        "en_name": "",
        "cn_name": "",
        "en_def": "",
        "cn_def": ""
    }


    file = open("html/hp_" + hpid + ".html", "rb")
    html = file.read()
    soup = BeautifulSoup(html, "html.parser")
    m = soup.select("main")
    c = m[0].find_all("div", {"class": "row_list"})

    try:
        p = c[0].select("p")
        dic[hpo]["url"] = url
        dic[hpo]["en_name"] = p[0].string.split("：")[1]
        dic[hpo]["cn_name"] = p[1].string.split("：")[1]
        dic[hpo]["en_def"] = p[2].string.split("：")[1]
        dic[hpo]["cn_def"] = p[3].string.split("：")[1]
        file.close()
    except:
        print("未找到描述:", hpo)
        file.close()
        os.remove("html/hp_" + hpid + ".html")

    return dic
    
hpo_dict = {}
hpoFile = open("hpolist.txt", "r")
for line in hpoFile:
	hpo = line.replace("\n", "")
	hpo_dict = chinahpo(hpo, hpo_dict)
hpoFile.close()

results = open("chpo.txt", "w", encoding="utf-8")
results.write("HPID\tEn\tCn\tEn_def\tCn_def\tURL\n")

for hd in hpo_dict.keys():
    print("正在写入", hd)
    output = [hd, hpo_dict[hd]["en_name"], hpo_dict[hd]["cn_name"], hpo_dict[hd]["en_def"], hpo_dict[hd]["cn_def"], hpo_dict[hd]["url"]]
    results.write("\t".join(output) + "\n")

results.close()
```

