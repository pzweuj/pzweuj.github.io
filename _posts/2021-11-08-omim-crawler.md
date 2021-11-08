---
title:  爬取OMIM数据库
tags: coding
---



[OMIM](https://omim.org)数据库收录了大量的基因与遗传病相关信息。

根据观察可通过https://omim.org/entry/\{mimID\}这个规律的链接进行爬取。因此先要找到mimID的列表，可以在[这个页面](https://omim.org/downloads)下载mim2gene.txt，解析后获得。

[OMIM的Robots规则](https://omim.org/robots.txt)允许google和bing爬取entry下的页面，允许的延迟是2s。因此爬取时，设置UA头为bingbot以及使用2s以上的延迟。

使用request进行爬取。因为未知要将结果处理为什么格式，所以不采用边爬取边解析的策略，而是爬取保存后，再进行解析。



### 爬取

以下是爬取python代码。其中，如果需要使用ip代理，使用的ip需要自行替换。可以建立一个ip池，再随机从ip池中抽取，下面的代码**未**涉及这步。另外增大延迟到4~10秒。已完成的结果记录到一个文件中，中断时可跳过已爬取结果，同时记录爬取每个结果的ip，方便从ip池中筛选出优质ip。

```python
import requests
import json
import random
import time

def readMim2Gene(filename):
    fileList = []
    with open(filename, "r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                strs = line.split("\t")
                mimNumber = strs[0]
                mimEntryType = strs[1]
                geneSymbol = "."
                if strs[3] != "":
                    geneSymbol = strs[3]
                if not "moved" in mimEntryType:
                    tempList = [mimNumber, mimEntryType, geneSymbol]
                    fileList.append(tempList)
    return fileList

def spider(mimNumber, ip_use):
    UA = [
        "Mozilla/5.0 (compatible; bingbot/2.0; +http://www.bing.com/bingbot.htm)",
        "Mozilla/5.0 AppleWebKit/537.36 (KHTML, like Gecko; compatible; bingbot/2.0; +http://www.bing.com/bingbot.htm) Chrome/96.0.4664.33 Safari/537.36 Edg/95.0.1020.40"
    ]
    url = "https://omim.org/entry/{}".format(mimNumber)
    proxy = {"https": "https://" + ip_use, "http": "http://" + ip_use}
    ua_use = random.choice(UA)
    header = {"User-Agent": ua_use}
    s = requests.session()
    if ip_use != "":
        s.proxies = proxy
    s.headers = header
    s.keep_alive = False
    
    try:
        html = s.get(url)    
        with open("entry/" + mimNumber + ".html", "w", encoding="utf-8") as f:
            f.write(html.text)
        ip_use_log = open("ip_use.txt", "a")
        ip_use_log.write(mimNumber + "\t" + ip_use + "\n")
        ip_use_log.close()
    except:
        print(mimNumber, "无法连接，无法获取")
        
def running(mim2gene, ip_use=""):
    idList = []
    for line in i:
        ids = line.replace("\n", "").split("\t")[0]
        idList.append(ids)
    mimList = readMim2Gene(mim2gene)
    i.close()
    for m in mimList:
        mimNumber = m[0]
        if not mimNumber in idList:
            print("正在爬取", mimNumber, ip_use)
            spider(mimNumber, ip_use)
            sleepTime = random.randint(4, 10)
            print("等待", sleepTime, "秒")
            time.sleep(sleepTime)
            
running("mim2gene.txt", "192.168.0.1:7777")
```



### 解析

使用下面脚本进行解析，获取信息

```python
import os
from bs4 import BeautifulSoup

def getTitle(soup, mimNum):
    title = soup.title.string.strip().split(mimNum + " - ")[1]
    return title
    
def getTable(soup):
    table = soup.table
    location = phenotype = mimNumber = inheritance = mappingKey = "."
    outputString = ""
    if table:
        trs = table.find_all("tr")
        for tr in trs:
            tds = tr.find_all("td")
            tds_len = len(tds)
            if tds_len == 5:
                gene = geneNum = "."
                location = tds[0].get_text().strip() if tds[0].get_text().strip() != "" else "."
                phenotype = tds[1].get_text().strip() if tds[1].get_text().strip() != "" else "."
                mimNumber = tds[2].get_text().strip() if tds[2].get_text().strip() != "" else "."
                inheritance = tds[3].get_text().strip() if tds[3].get_text().strip() != "" else "."
                mappingKey = tds[4].get_text().strip() if tds[4].get_text().strip() != "" else "."
            elif tds_len == 4:
                location = gene = geneNum = "."
                phenotype = tds[0].get_text().strip() if tds[0].get_text().strip() != "" else "."
                mimNumber = tds[1].get_text().strip() if tds[1].get_text().strip() != "" else "."
                inheritance = tds[2].get_text().strip() if tds[2].get_text().strip() != "" else "."
                mappingKey = tds[3].get_text().strip() if tds[3].get_text().strip() != "" else "."
            elif tds_len == 7:
                location = tds[0].get_text().strip() if tds[0].get_text().strip() != "" else "."
                phenotype = tds[1].get_text().strip() if tds[1].get_text().strip() != "" else "."
                mimNumber = tds[2].get_text().strip() if tds[2].get_text().strip() != "" else "."
                inheritance = tds[3].get_text().strip() if tds[3].get_text().strip() != "" else "."
                mappingKey = tds[4].get_text().strip() if tds[4].get_text().strip() != "" else "."
                gene = tds[5].get_text().strip() if tds[5].get_text().strip() != "" else "."
                geneNum = tds[6].get_text().strip() if tds[6].get_text().strip() != "" else "."
            else:
                continue
            outputString = outputString + "[" + "|".join([location, phenotype, mimNumber, inheritance, mappingKey, gene, geneNum]) +"];"
    return outputString.rstrip(";") if outputString != "" else "."
    
def getClinicalFold(soup):
    clinicalSynopsisFoldList = soup.select("#clinicalSynopsisFold")
    childStringDiv = "."
    if len(clinicalSynopsisFoldList) != 0:
        clinicalDivList = clinicalSynopsisFoldList[0].select(".small")[0]
        childDiv = clinicalDivList.find_all("div", recursive=False)
        childStringDiv = ""
        for i in range(len(childDiv) - 1):
            childDivList = childDiv[i].get_text().replace("\n", "").strip().split(" - ")
            headerDiv = childDivList[0].strip()
            stringDiv = ""
            for i in range(len(childDivList)):
                if i != 0:
                    stringDiv = stringDiv + childDivList[i].strip() + "; "
            stringDiv = headerDiv + ": " + stringDiv.rstrip("; ")
            childStringDiv = childStringDiv + stringDiv + " | "
        childStringDiv = childStringDiv.rstrip(" | ")
    return childStringDiv
    
def getDescription(soup):
    descriptionFoldList = soup.select("#descriptionFold")
    descriptionFold = "." if len(descriptionFoldList) == 0 else descriptionFoldList[0].get_text().strip().replace("\n", " \\ ")
    return descriptionFold
    
def main(inputDir, outputFile):
    outputs = open(outputFile, "w", encoding="utf-8")
    outputs.write("#MimNum\tTitle\tTable\tDescription\tClinical\tURL\n")
    for i in os.listdir(inputDir):
        if ".html" in i:
            mimNum = i.split(".html")[0]
            url = "https://omim.org/entry/{}".format(mimNum)
            html = open(inputDir + "/{}.html".format(mimNum), "r", encoding="utf-8")
            soup = BeautifulSoup(html, "lxml")
            html.close()
            title = getTitle(soup, mimNum)
            table = getTable(soup)
            clinical = getClinicalFold(soup)
            des = getDescription(soup)
            output = "\t".join([mimNum, title, table, des, clinical, url])
            outputs.write(output + "\n")
    outputs.close()
    
main("entry", "output.txt")
```



### 参考

[网络爬虫学习笔记](https://blog.csdn.net/fanyucai1/article/details/103987949)

[抓取omim数据库数据](https://www.fee.im/2020/03/crawling-data-from-omim/)