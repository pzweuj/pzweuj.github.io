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
cut genes_to_phenotype.txt -f3 | sort | uniq > hpolist.txt
```

最终其实只有8000左右的HP编号有用。下面是爬取代码，后期可以使用Thread和queue等库来增加线程数，获得更快速度。**但是爬网站还是不能太过分**，设置一个长一点的间隔，只用单线程。爬完再解析。

```python
from msedge.selenium_tools import Edge
from msedge.selenium_tools import EdgeOptions
from bs4 import BeautifulSoup
import os
import time
import random
```


设置一个UA池
```python
def randomUA():
    MY_USER_AGENT = [
        "Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; AcooBrowser; .NET CLR 1.1.4322; .NET CLR 2.0.50727)",
        "Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 6.0; Acoo Browser; SLCC1; .NET CLR 2.0.50727; Media Center PC 5.0; .NET CLR 3.0.04506)",
        "Mozilla/4.0 (compatible; MSIE 7.0; AOL 9.5; AOLBuild 4337.35; Windows NT 5.1; .NET CLR 1.1.4322; .NET CLR 2.0.50727)",
        "Mozilla/5.0 (Windows; U; MSIE 9.0; Windows NT 9.0; en-US)",
        "Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; Win64; x64; Trident/5.0; .NET CLR 3.5.30729; .NET CLR 3.0.30729; .NET CLR 2.0.50727; Media Center PC 6.0)",
        "Mozilla/5.0 (compatible; MSIE 8.0; Windows NT 6.0; Trident/4.0; WOW64; Trident/4.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; .NET CLR 1.0.3705; .NET CLR 1.1.4322)",
        "Mozilla/4.0 (compatible; MSIE 7.0b; Windows NT 5.2; .NET CLR 1.1.4322; .NET CLR 2.0.50727; InfoPath.2; .NET CLR 3.0.04506.30)",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; zh-CN) AppleWebKit/523.15 (KHTML, like Gecko, Safari/419.3) Arora/0.3 (Change: 287 c9dfb30)",
        "Mozilla/5.0 (X11; U; Linux; en-US) AppleWebKit/527+ (KHTML, like Gecko, Safari/419.3) Arora/0.6",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.2pre) Gecko/20070215 K-Ninja/2.1.1",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; zh-CN; rv:1.9) Gecko/20080705 Firefox/3.0 Kapiko/3.0",
        "Mozilla/5.0 (X11; Linux i686; U;) Gecko/20070322 Kazehakase/0.4.5",
        "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.8) Gecko Fedora/1.9.0.8-1.fc10 Kazehakase/0.5.6",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/535.11 (KHTML, like Gecko) Chrome/17.0.963.56 Safari/535.11",
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_7_3) AppleWebKit/535.20 (KHTML, like Gecko) Chrome/19.0.1036.7 Safari/535.20",
        "Opera/9.80 (Macintosh; Intel Mac OS X 10.6.8; U; fr) Presto/2.9.168 Version/11.52",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/536.11 (KHTML, like Gecko) Chrome/20.0.1132.11 TaoBrowser/2.0 Safari/536.11",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.71 Safari/537.1 LBBROWSER",
        "Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E; LBBROWSER)",
        "Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; QQDownload 732; .NET4.0C; .NET4.0E; LBBROWSER)",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/535.11 (KHTML, like Gecko) Chrome/17.0.963.84 Safari/535.11 LBBROWSER",
        "Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E)",
        "Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E; QQBrowser/7.0.3698.400)",
        "Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; QQDownload 732; .NET4.0C; .NET4.0E)",
        "Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 5.1; Trident/4.0; SV1; QQDownload 732; .NET4.0C; .NET4.0E; 360SE)",
        "Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; QQDownload 732; .NET4.0C; .NET4.0E)",
        "Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E)",
        "Mozilla/5.0 (Windows NT 5.1) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.89 Safari/537.1",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.89 Safari/537.1",
        "Mozilla/5.0 (iPad; U; CPU OS 4_2_1 like Mac OS X; zh-cn) AppleWebKit/533.17.9 (KHTML, like Gecko) Version/5.0.2 Mobile/8C148 Safari/6533.18.5",
        "Mozilla/5.0 (Windows NT 6.1; Win64; x64; rv:2.0b13pre) Gecko/20110307 Firefox/4.0b13pre",
        "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:16.0) Gecko/20100101 Firefox/16.0",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11",
        "Mozilla/5.0 (X11; U; Linux x86_64; zh-CN; rv:1.9.2.10) Gecko/20100922 Ubuntu/10.10 (maverick) Firefox/3.6.10",
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36",
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/83.0.4103.53 Safari/537.36"
    ]
    ua = random.choice(MY_USER_AGENT)
    return ua
```

再设置一个IP池。如果没有比较多的IP，建议放弃。
```python
def randomIP():
    ips = open("ip.txt", "r")
    ip = []
    for line in ips:
        l = line.replace("\n", "")
        ip.append(l)
    ip_random = random.choice(ip)
    return ip_random
```

然后随机导入。基本形成随机IP和UA的访问。

```python
# 爬取
def chinahpo(hpo):
    s = random.randint(5, 10)
    print("等待 " + str(s) + "秒")
    time.sleep(s)
    ip = randomIP()
    print("使用IP " + ip)
    options = EdgeOptions()
    options.use_chromium = True
    # options.add_argument("headless")
    # options.add_argument("disable-gpu")
    options.add_argument("--proxy-server={ip}".format(ip=ip))
    options.add_argument("--disable-blink-features")
    options.add_argument("--disable-blink-features=AutomationControlled")
    options.add_argument("start-maximized")
    options.add_experimental_option("excludeSwitches", ["enable-automation"])
    options.add_experimental_option("useAutomationExtension", False)
    msedge = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedgedriver.exe"

    driver = Edge(options=options, executable_path=msedge)
    script = "Object.defineProperty(navigator, 'webdriver', {get: () => undefined})"
    driver.execute_script(script)
    UA = randomUA()
    driver.execute_cdp_cmd("Network.setUserAgentOverride", {"userAgent": UA})
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

