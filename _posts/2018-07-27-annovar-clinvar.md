---
title: annovar的官方clinvar数据库格式修改教程
tags: default
---
刚刚才发现annovar在7月8日的时候放出了官方是怎么把clinvar转换成annovar格式的流程。

[首先点击这里下载脚本](http://www.openbioinformatics.org/annovar/download/prepare_annovar_user.pl)!
这个脚本其实之前也能下，是用来转换cosmic数据库的。

然后安装小工具包[VT](https://genome.sph.umich.edu/wiki/Vt)
```bash
git clone https://github.com/atks/vt.git
cd vt
make
./vt -h
```

开始测试。
首先去clinvar下载最新的数据库。
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20180701.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20180701.vcf.gz.tbi
```
然后跟着流程来。
```bash
vt decompose clinvar_20180701.vcf.gz -o temp.split.vcf
perl prepare_annovar_user.pl -dbtype clinvar_preprocess2 temp.split.vcf -out temp.split2.vcf
vt normalize temp.split2.vcf \
	-r ~/database/b37/Homo_sapiens.GRCh37.dna.toplevel.fa \
	-o temp.norm.vcf \
	-w 2000000
perl prepare_annovar_user.pl -dbtype clinvar2 temp.norm.vcf -out hg19_clinvar_20180701.txt
```
最后还是按之前那样创建个索引文件吧。
[传送门](https://pzweuj.github.io/2018/04/25/convert-clinvar-to-annovar.html)


[-_-]:123