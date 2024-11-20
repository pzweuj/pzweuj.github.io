---
title: 爬取oncokb
tags: coding
---

爬虫务必设置一个大的时间间隔。

确认爬取路线：

1，从网站中获取所有基因的列表；

2，爬取基因页面；

3，解析基因页面并获得位点列表；

4，爬取位点页面；

5，解析并整理所有内容。

### 基因列表
oncokb在[这个页面](https://www.oncokb.org/cancerGenes)中提供了基因列表的下载。



### 爬取基因页面

oncokb的页面URL结构非常整洁，照旧使用playwright框架来爬取，由于每个页面下可能有多个标签，标签的内容需要点击标签才会加载出来，因此这里每个标签都尝试点击一次并保存一次。

```python
from playwright.sync_api import Playwright, sync_playwright, expect
import time
import random

def run(playwright: Playwright) -> None:
    # 基于列表
    finished = []
    finished_num = 0
    with open("finished.txt", "r", encoding="utf-8") as f:
        for line in f:
            finished.append(line.strip())
            finished_num += 1
    
    genelist = []
    with open("cancerGeneList.tsv", "r", encoding="utf-8") as g:
        for line in g:
            if not line.startswith("Hugo Symbol"):
                lineList = line.strip().split("\t")
                if not lineList in finished:
                    genelist.append(lineList[0])
    # 随机打乱列表
    random.shuffle(genelist)

    browser = playwright.chromium.launch(headless=True)
    context = browser.new_context()
    
    # 注入反反爬脚本
    page = context.new_page()
    page.add_init_script(path="stealth.min.js")
    
    # 基因遍历
    for gene in genelist:
        page.wait_for_load_state("domcontentloaded")
        page.goto("https://www.oncokb.org/gene/{gene}".format(gene=gene), timeout=0)
        time.sleep(random.randint(12, 20))
        try:
            page.get_by_role("tab", name="Annotated Alterations").click()
            time.sleep(2)
            html = page.content()
            with open("html/gene/" + gene + ".alterations.html", "w", encoding="utf-8") as f:
                f.write(html)
        except:
            pass
        try:
            page.get_by_role("tab", name="Therapeutic").click()
            time.sleep(2)
            html = page.content()
            with open("html/gene/" + gene + ".therapeutic.html", "w", encoding="utf-8") as f:
                f.write(html)
        except:
            pass
        try:
            page.get_by_role("tab", name="Diagnostic").click()
            time.sleep(2)
            html = page.content()
            with open("html/gene/" + gene + ".diagnostic.html", "w", encoding="utf-8") as f:
                f.write(html)
        except:
            pass
        try:
            page.get_by_role("tab", name="Prognostic").click()
            time.sleep(2)
            html = page.content()
            with open("html/gene/" + gene + ".prognostic.html", "w", encoding="utf-8") as f:
                f.write(html)
        except:
            pass
        try:
            page.get_by_role("tab", name="FDA-Recognized Content").click()
            time.sleep(2)
            html = page.content()
            with open("html/gene/" + gene + ".fda.html", "w", encoding="utf-8") as f:
                f.write(html)
        except:
            pass
        html = page.content()
        with open("html/gene/" + gene + ".empty.html", "w", encoding="utf-8") as f:
            f.write(html)

        time.sleep(0.3)
        print("【已完成】", gene)

        with open("finished.txt", "a", encoding="utf-8") as f:
            f.write(gene + "\n")

        # 加大爬取间隔，防止被禁
        time.sleep(random.randint(20, 30))
    
    page.close()
    context.close()
    browser.close()

with sync_playwright() as playwright:
    run(playwright)

```



### 解析基因页面获得位点列表
使用BeautifulSoup进行解析

```python
import os
from bs4 import BeautifulSoup

# 基因基本信息
def gene_base_info(gene):
    htmlFile = open("html/gene/{gene}.empty.html".format(gene=gene), "r", encoding="utf-8")
    htmlHandle = htmlFile.read()
    soup = BeautifulSoup(htmlHandle, "html.parser")
    div_ele = soup.find_all("div", role="alert")
    negative = False
    for d in div_ele:
        if "We do not have any information for this gene" in d.text:
            negative = True
    if negative:
        outputList = [gene, "无此基因记录", "", "", "", "", ""]
    else:
        geneInfoTable = soup.find_all("table")[0]
        tds = geneInfoTable.find_all("td")
        tdDict = {}
        for n in range(len(tds)):
            if n % 2 == 0:
                tdDict[tds[n].text] = tds[n+1].text
        ncbi = tdDict["NCBI Gene"]
        if "Ensembl Gene" in tdDict:
            embl = tdDict["Ensembl Gene"].split(" (")[0]
        else:
            embl = "-"
        grch37 = grch38 = "-"
        if "Location" in tdDict:
            if "GRch37" in tdDict["Location"]:
                grch37 = tdDict["Location"].split(" (GRch37)")[0].replace("Chr", "chr")
            if "GRch38" in tdDict["Location"]:
                if "GRch37" in tdDict["Location"]:
                    grch38 = tdDict["Location"].split(" (GRch37)")[1].rstrip(" (GRch38)").replace("Chr", "chr")
                else:
                    grch38 = tdDict["Location"].split(" (GRch38)")[0].replace("Chr", "chr")
        if "Ensembl Transcript" in tdDict:
            embl_t = tdDict["Ensembl Transcript"].split(" (")[0]
        else:
            embl_t = "-"
        if "RefSeq" in tdDict:
            refseq = tdDict["RefSeq"].split(" (")[0]
        else:
            refseq = "-"
        outputList = [gene, ncbi, embl, grch37, grch38, embl_t, refseq]
    htmlFile.close()
    return outputList

# alterations
def alterations_info(gene):
    htmlFile = open("html/gene/{gene}.alterations.html".format(gene=gene), "r", encoding="utf-8")
    htmlHandle = htmlFile.read()
    soup = BeautifulSoup(htmlHandle, "html.parser")
    rows = soup.find_all("div", role="row")
    outputList = []
    for i in rows[1:]:
        divs = i.find_all("div")
        outputList.append([gene, divs[0].text, divs[1].text, divs[2].text])
    htmlFile.close()
    return outputList

# diagnostic
def diagnostic_info(gene):
    htmlFile = open("html/gene/{gene}.diagnostic.html".format(gene=gene), "r", encoding="utf-8")
    htmlHandle = htmlFile.read()
    soup = BeautifulSoup(htmlHandle, "html.parser")
    rows = soup.find_all("div", role="row")
    outputList = []
    for i in rows[1:]:
        divs = i.find_all("div")
        level = divs[0].span.i.get("class")[-1].lstrip("level-")
        output = [gene, level, divs[2].text, divs[3].text]
        outputList.append(output)
    htmlFile.close()
    return outputList

# therapeutic
def therapeutic_info(gene):
    htmlFile = open("html/gene/{gene}.therapeutic.html".format(gene=gene), "r", encoding="utf-8")
    htmlHandle = htmlFile.read()
    soup = BeautifulSoup(htmlHandle, "html.parser")
    rows = soup.find_all("div", role="row")
    outputList = []
    for i in rows[1:]:
        divs = i.find_all("div")
        level = divs[0].span.i.get("class")[-1].lstrip("level-")
        output = [gene, level, divs[2].text, divs[3].text, divs[4].text]
        outputList.append(output)
    htmlFile.close()
    return outputList

# prognostic
def prognostic_info(gene):
    htmlFile = open("html/gene/{gene}.prognostic.html".format(gene=gene), "r", encoding="utf-8")
    htmlHandle = htmlFile.read()
    soup = BeautifulSoup(htmlHandle, "html.parser")
    rows = soup.find_all("div", role="row")
    outputList = []
    for i in rows[1:]:
        divs = i.find_all("div")
        level = divs[0].span.i.get("class")[-1].lstrip("level-")
        output = [gene, level, divs[2].text, divs[3].text]
        outputList.append(output)
    htmlFile.close()
    return outputList

# fda
def fda_info(gene):
    htmlFile = open("html/gene/{gene}.fda.html".format(gene=gene), "r", encoding="utf-8")
    htmlHandle = htmlFile.read()
    soup = BeautifulSoup(htmlHandle, "html.parser")
    rows = soup.find_all("div", role="row")
    outputList = []
    for i in rows[1:]:
        divs = i.find_all("div")
        output = [gene, divs[1].text, divs[2].text, divs[3].text]
        outputList.append(output)
    return outputList

####################
geneList = []
with open("oncokb.all.txt", "r", encoding="utf-8") as g:
    for line in g:
        geneList.append(line.strip())

output1 = open("output_all/output.summary.txt", "w", encoding="utf-8")
output2 = open("output_all/output.anno.txt", "w", encoding="utf-8")
output3 = open("output_all/output.thera.txt", "w", encoding="utf-8")
output4 = open("output_all/output.diagn.txt", "w", encoding="utf-8")
output5 = open("output_all/output.progn.txt", "w", encoding="utf-8")
output6 = open("output_all/output.fda.txt", "w", encoding="utf-8")
htmlList = os.listdir("html/gene")
for gene in geneList:
    print(gene)
    for i in ["empty", "alterations", "diagnostic", "therapeutic", "prognostic", "fda"]:
        if (gene + "." + i + ".html") in htmlList:
            if i == "empty":
                output1.write("\t".join(gene_base_info(gene)) + "\n")
            elif i == "alterations":
                for j in alterations_info(gene):
                    output2.write("\t".join(j) + "\n")
            elif i == "therapeutic":
                for j in therapeutic_info(gene):
                    output3.write("\t".join(j) + "\n")
            elif i == "diagnostic":
                for j in diagnostic_info(gene):
                    output4.write("\t".join(j) + "\n")
            elif i == "prognostic":
                try:
                    for j in prognostic_info(gene):
                        output5.write("\t".join(j) + "\n")
                except:
                    continue
            elif i == "fda":
                try:
                    for j in fda_info(gene):
                        output6.write("\t".join(j) + "\n")
                except:
                    continue
output1.close()
output2.close()
output3.close()
output4.close()
output5.close()
output6.close()
```



### 爬取位点

与爬取基因大差不差，修改一下上面的脚本，然后位点的某些字符特殊处理一下即可。

### 解析位点页面

其实感觉基因页面的信息也够用了，但是实际上位点页面中，评级的信息会比基因页面的信息多，另外位点页面中会包含位点的解析。因此这里就只把需要的突变信息和评级抓出来了。

```python
from bs4 import BeautifulSoup

# 基因基本信息
def variant_base_info(gene, variant):
    # MET	981_1028splice
    # gene = "MET"
    # variant = "981_1028splice"
    oncoFixResult = variantDescription = geneDescription = "-"
    htmlFile = open("html/variant/{gene}_{variant}.empty.html".format(gene=gene, variant=variant), "r", encoding="utf-8")
    htmlHandle = htmlFile.read()
    soup = BeautifulSoup(htmlHandle, "html.parser")
    try:
        geneDescription = soup.find_all("div", attrs={"class": "mb-3"})[0].text.strip()
    except:
        pass
    variantDescriptionList = soup.find_all("span")
    try:
        n = 0
        for vars in variantDescriptionList:
            variant = variant.replace("Fusion", "fusion").replace("Amplification", "Amplification of ").replace("Ter", "*")
            if "PMID" in vars.text:
                n += 1
                if n == 1:
                    variantDescription = vars.text
    except:
        pass
    oncoList = soup.find_all("h5")
    try:
        for onco in oncoList[0].find_all("span"):
            if not "style" in str(onco):
                if onco.text in ["Oncogenic", "Resistance", "Likely Oncogenic", "Likely Neutral", "Inconclusive"]:
                    oncoFixResult = onco.text
    except:
        pass


    return [geneDescription.replace("\n", "\\x0a"), variantDescription.replace("\n", "\\x0a"), oncoFixResult]



####################
input = open("oncokb.variant.txt", "r", encoding="utf-8")
output = open("output_all/oncokb.variant.des.txt", "w", encoding="utf-8")
for line in input:
    lines = line.rstrip().split("\t")
    gene = lines[0]
    variant = lines[1].replace("*", "Ter")
    print(gene, variant)
    geneDes, variantDes, oncoFixResult = variant_base_info(gene, variant)
    output.write("\t".join([gene, variant.replace("Ter", "*"), geneDes, variantDes, oncoFixResult]) + "\n")
output.close()
input.close()
    
```

