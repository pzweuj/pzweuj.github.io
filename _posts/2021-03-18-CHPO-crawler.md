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

最终其实只有23万左右的HP编号有用。下面是爬取代码，后期可以使用Thread和queue等库来增加线程数，获得更快速度。**但是爬网站还是不能太过分**，设置一个长一点的间隔，只用单线程。爬完再解析。

```python
from msedge.selenium_tools import Edge
from msedge.selenium_tools import EdgeOptions
from bs4 import BeautifulSoup
import os
import time
import random

# 爬取
def chinahpo(hpo):
    time.sleep(random.randint(5, 30))
    options = EdgeOptions()
    options.use_chromium = True
    # options.add_argument("headless")
    # options.add_argument("disable-gpu")
    options.add_argument("--disable-blink-features")
    options.add_argument("--disable-blink-features=AutomationControlled")
    options.add_argument("start-maximized")
    options.add_experimental_option("excludeSwitches", ["enable-automation"])
    options.add_experimental_option("useAutomationExtension", False)
    msedge = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedgedriver.exe"

    driver = Edge(options=options, executable_path=msedge)
    script = "Object.defineProperty(navigator, 'webdriver', {get: () => undefined})"
    driver.execute_script(script)
    driver.execute_cdp_cmd("Network.setUserAgentOverride", {"userAgent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/83.0.4103.53 Safari/537.36"})
    print(driver.execute_script("return navigator.userAgent;"))

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
    fin = open("finish.txt", "a")
    fin.write(hpo + "\n")
    fin.close()

# 解析
def analysis(hpo):
    hpid = hpo.split(":")[1]
    file = open("html/hp_" + hpid + ".html", "rb")
    html = file.read()
    soup = BeautifulSoup(html, "html.parser")
    m = soup.select("main")
    c = m[0].find_all("div", {"class": "row_list"})
    url = "http://www.chinahpo.org/#/searchList?trigger=1&tabType=1&searchContent=HP%3A{hpid}".format(hpid=hpid)

    output = [hpo]

    try:
        p = c[0].select("p")
        output.append(p[0].string.split("：")[1])
        output.append(p[1].string.split("：")[1])
        output.append(p[2].string.split("：")[1])
        output.append(p[3].string.split("：")[1])
        output.append(url)
        file.close()
    except:
        output.append("未找到信息")
        output.append("-")
        output.append("-")
        output.append("-")
        output.append(url)

        # os.remove("html/hp_" + hpid + ".html")
        filtered = open("filtered.txt", "a")
        filtered.write(hpo + "\n")
        filtered.close()

    return "\t".join(output)


# 爬取
hpoFile = open("hpolist.txt", "r")
for line in hpoFile:
    hpo = line.replace("\n", "")
    chinahpo(hpo)
hpoFile.close()


# 解析
hpoFileFinish = open("finish.txt", "r")
r = open("chinahpo.txt", "w", encoding="utf-8")
r.write("HPOID\tEN\tCN\tEN_des\tCN_des\tURL\n")
for line in hpoFileFinish:
    hpo = line.replace("\n", "")
    try:
        hpoString = analysis(hpo)
        r.write(hpoString + "\n")
    except:
        pass
r.close()
hpoFileFinish.close()
```

