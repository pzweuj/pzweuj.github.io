---
title: 截断引物的需求
tags: coding
---

得到一个需求，截断fastq里的引物序列，为什么是截断而不是去除，目的是为了交付数据时不被同行得到引物同时读长看起来也没那么短。一开始想象中是很简单的，只要替换fastq中的引物序列就好了，但是后来发现有些目标区域其实也含有引物的序列，如果直接替换了会导致这些目标区域变成缺失之类的判读。

然后新增了一个要求，在截断的同时需要保留adapter。有一个做法是先除adapter，然后再5‘和3’截断10bp，再把adapter加回去。但是我觉得这样的难点在于adapter的质量值也得加回去，比较复杂。

所以，最终的写法是，先把序列和质量值的两行zip起来，然后搜索里面的adapter序列，再把adapter前去掉10bp，如果搜索不到adapter，就直接3‘去10bp。5’的10bp可以直接去除。

弱智写法，大佬绕道。比较吃内存。
```python
## 读入
import gzip

rawFastq = gzip.open("xxx_1.fq.gz", "rb")
cleanFastq = open("xxx_1.clean.fq.gz", "w")

## 先把adapter放这里
r1adpt = "AAAGAGCATTCAAAGTGTCAA"
r2adpt = "GATCGTCGGACTGTAGAACT"
adpt = r1adpt # 这里可以设置一个传参是r1还是r2

## 因为序列中adapter可能不是全长，所以定了一个10bp的长度来找adapter
## 得到所有10bp的情况，太短的话怕特异性差
a = 0
adptList = []
while (m + 10) <= len(adpt):
	adptBin = adpt[m: m+10]
	m += 1
	adptList.append(adptBin)

## 将每一行的类型分类
c = 0
name = []
seq = []
direct = []
qual = []
for line in rawFastq:
	c += 1
	k = c % 4
	if k == 1:
		name.append(line)
	elif k == 2:
		seq.append(line)
	elif k == 3:
		direct.append(line)
	elif k == 0:
		qual.append(line)
	else:
		continue
rawFastq.close()

## 重点来了，处理seq和qual
seq_qual = zip(seq, qual)
seq_new = []
qual_new = []
for sq in seq_qual:
	s = 0
	
	## 先判断有没有adapter序列，判断和处理分开，不然有点乱
	for adpts in adptList:
		if adpts in sq[0]:
			s += 1
	
	## 没有找到adapter时去除前10后10
	if s == 0:
		seq_new.append(sq[0][10:][:-10] + "\n")
		qual_new.append(sq[1][10:][:-10] + "\n")
	
	## 找到的情况
	else:
		for adpts in adptList:
			if adpts in sq[0]:
				
				# 直接剪切adapter
				adptSplit = sq[0][10:].split(adpts)
				# 保存下来剪剩下多少
				cutLength = len(adptSplit[0][:-10])
				# 剪掉adapter前面的10bp然后把adapter和后面的加回去
				seq_new.append(adptSplit[0][:-10] + adpts + adptSplit[1])
				# 对应的质量值也是这样剪
				qual_new.append(sq[1][10:][:cut] + sq[1][10:][cut+10:])
				# 然后第一次遇到对的adapter就可以跳出循环了，避免多个adapter又切几次
				break

# 写入
fi = zip(name, seq_new, direct, qual_new)
for finals in fi:
	cleanFastq.write(finals[0])
	cleanFastq.write(finals[1])
	cleanFastq.write(finals[2])
	cleanFastq.write(finals[3])
	
cleanFastq.close()
print "task done"
```




[-_-]:love jing