---
title: 爬取CKB数据库
tags: coding
---

Jax Ckb是比较常用的肿瘤相关位点注释的数据库。[曾几何时](https://pzweuj.github.io/2018/06/10/ckb-jax.html)数据库还是免费查询的，现在部分基因如EGFR已经转到了高级版本中，不过免费版还是能查若干基因的。

参考文章[python多线程爬取CKB数据库](https://blog.csdn.net/dujidan/article/details/105472604)对CKB数据库的核心基因位点注释进行爬取，修改了其中解析部分，以及补充上位点描述。

```python
import requests
import bs4
from bs4 import BeautifulSoup
import threadpool


# 获得网页源码
def getHTMLText(url):
    try:
        r = requests.get(url, timeout=40)
        r.raise_for_status()
        r.encoding = r.apparent_encoding
        return r.text
    except:
        return ""

# 获得ckb基因对应网页字典
# 例 {"ALK": "https://ckb.jax.org/gene/show?geneId=238"}
def get_gene_id(url):
    html = getHTMLText(url)
    soup = BeautifulSoup(html, "html.parser")
    gene_id_dict = {}
    for a in soup.find_all(name="a", attrs="btn btn-default btn-gene btn-block"):
        gene_name = a.string.replace("\n", "").replace(" ", "")
        ID = a.attrs["href"]
        gene_id = "https://ckb.jax.org" + ID
        gene_id_dict[gene_name] = gene_id
    return gene_id_dict

# 获得位点对应连接
# 例 [["ALK", "A380T", "https://ckb.jax.org/geneVariant/show?geneVariantId=39002"]]
def gene_variant_link(gene, url):
    list_link = []
    url = url
    html = getHTMLText(url)
    soup = BeautifulSoup(html, "html.parser")
    for a in soup.select('a[href^="/geneVariant"]'):
        list_link.append([gene, a.text.replace(" ", "").replace("\n", ""), "https://ckb.jax.org" + a["href"]])
    return list_link

# 解析
def print_gene_variant(links):
    url = links[-1]
    html = getHTMLText(url)
    soup = BeautifulSoup(html, "html.parser")

    des_td = []
    for td in soup.find_all("td"):
        des_td.append(td.text)

    description = "-"
    drug_res = "-"
    for des in range(len(des_td)):
        if "Gene Variant Descriptions" in des_td[des]:
            description = des_td[des + 1]
        if "Associated Drug Resistance" in des_td[des]:
            drug_res = des_td[des + 1]

    if len(soup.find_all("table", attrs={"id": "TranscriptTabTable"})) == 2:
        for tr in soup.find_all("table", attrs={"id": "TranscriptTabTable"})[1].children:
            if isinstance(tr, bs4.element.Tag):
                results = tr.text.replace(" ", "").replace("\n\n", "").strip()
                results = results.replace("\n", "\t")
                if not "GenomeBuild" in results:
                    results = results.split("\tGRCh38/hg38")[0] + "\tGRCh38/hg38"
                    results = "\t".join(links) + "\t" + description  + "\t" + drug_res + "\t" + results
                    print(results)
                    
    if len(soup.find_all("table", attrs={"id": "TranscriptTabTable"})) == 1:
        for tr in soup.find_all("table", attrs={"id": "TranscriptTabTable"})[0].children:
            if isinstance(tr, bs4.element.Tag):
                results = tr.text.replace(" ", "").replace("\n\n", "").strip()
                results = results.replace("\n", "\t")
                if not "GenomeBuild" in results:
                    results = results.split("\tGRCh38/hg38")[0] + "\tGRCh38/hg38"
                    results = "\t".join(links) + "\t" + description + "\t" + drug_res + "\t" + results
                    print(results)

    if len(soup.find_all("table", attrs={"id": "TranscriptTabTable"})) == 0:
        results = "NAN"
        results = "\t".join(links) + "\t" + description + "\t" + drug_res + "\t" + results
        print(results)


url = "https://ckb.jax.org/gene/grid"
gDict = get_gene_id(url)

print("Gene\tVariant\tUrl\tDescriptions\tDrug Resistance\tTranscript\tgDNA\tcDNA\tProtein\tSourceDatabase\tGenomeBuild")
for name, ID in gDict.items():
    url = ID
    gene = name
    list_link = gene_variant_link(gene, url)

    pool = threadpool.ThreadPool(10)
    tasks = threadpool.makeRequests(print_gene_variant, list_link)
    [pool.putRequest(task) for task in tasks]
    pool.wait()
```

最终会得到这种结构
```
CDH1    L583R   https://ckb.jax.org/geneVariant/show?geneVariantId=2454 CDH1 L583R lies within the Cadherin domain 4 of the Cdh1 protein (UniProt.org). L583R results in aberrant maturation and decreased stability of Cdh1 protein, and leads to scattered cell distribution, indicating reduction of Cdh1 adhesive function in cell culture (PMID: 22470475).         NM_004360.4     chr16:g.68822037T>G      c.1748T>G       p.L583R RefSeq  GRCh38/hg38
```

然后需要hg19的话，对chr16:g.68822037T>G进行解析再用liftover就行了。