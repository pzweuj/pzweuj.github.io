---
title: 把clinvar转换成annovar可用的格式
tags: coding
---

annovar可以说是最常用的注释软件了。可是官方的数据库更新很慢，所以，最好是自己更新。
学会了下面的操作，自建数据库用annovar注释也不是问题。

首先我们需要下载最新的clinvar数据库。
点击进入[GRCh37](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/)版本的FTP地址。
我下载的是20180401的Clinvar。
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20180401.vcf.gz
gunzip clinvar_20180401.vcf.gz
```

解压后，用下面的python脚本进行解析。
```python
a = open('clinvar_20180401.vcf', 'r')
b = open('hg19_clinvar_20180401.txt', 'w')

b.write('#chrom	start	end	ref	alt	CLNSIG	CLNDN\n')

for line in a:
	if line.startswith('#'):
		continue
	else:
		lines = line.split('\t')
		chrom = lines[0]
		start = lines[1]
		ref = lines[3]
		alt = lines[4]
		end = str(int(start) + len(ref) - 1)
		clinvar = lines[7]

		clins = clinvar.split(';')

		for i in clins:
			if i.startswith('CLNDN'):
				CLNDN = i.split('=')[1]
			elif i.startswith('CLNSIG'):
				CLNSIG = i.split('=')[1]
			else:
				continue

		l = [chrom, start, end, ref, alt, CLNSIG, CLNDN]
		s = '\t'.join(l) + '\n'
		b.write(s)
b.close()
a.close()
```

然后就可以放进annovar的数据库中，注释的时候填入clinvar_20180401就可以了。
当然，这还不够快，接下来是建立索引文件。
[这里](https://github.com/pzweuj/practice/blob/master/perl/makeAnnovarIndex/makeAnnovarIndex.pl)有一个建立annovar索引的perl程序。
运行的姿势是
```
perl makeAnnovarIndex.pl hg19_clinvar_20180401.txt 1000
```
这样就会生成一个 hg19_clinvar_20180401.txt.idx 索引文件。

[T_T]:###############