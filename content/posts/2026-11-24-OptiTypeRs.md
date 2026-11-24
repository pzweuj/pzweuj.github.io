---
title: OptiTypeRs，用 Rust 重构 OptiType
tags: software
---

2021 年我写过一篇 [HLA 分型软件 OptiType](/posts/Optitype) 的使用记录。那时候如果不想用 Docker，跑通一次要挨个配 Python 依赖、RazerS3、samtools，再找一个求解器塞进 `config.ini` 里。折腾完这一大圈，其实也就为了从两份 FASTQ 里拿到 HLA-A、HLA-B、HLA-C 的两字段分型（比如 `A*01:01`）。

这次我把 OptiType 的核心用 Rust 重写了一遍，起名叫 [OptiTypeRs](https://github.com/pzweuj/OptiTypeRs)。它保留了 [原版 OptiType](https://github.com/FRED-2/OptiType) 将 A、B、C 放在一起联合求解的整数规划模型，但换掉了外部比对和 Python 那套求解流程——比对在程序内部完成，整数规划则静态链接了 HiGHS。现在运行时只需要一个 `optitype-rs` 二进制文件外加一份本地参考库，输出依然是 HLA-A、HLA-B、HLA-C 的两字段分型。

---

## 两条 read 怎么投票

HLA 区域序列高度相似，一条 read 同时比对上好几个等位基因太正常了。如果先随便挑一个命中或者简单按支持数排序，很容易把本来共享的模糊证据全算到某一个型别头上。OptiType 的思路是把 A、B、C 的候选一起塞进整数规划模型，联合选出能解释这些 reads 的组合，并用惩罚项抑制为了多覆盖几条 reads 而强行增加杂合型的倾向。

Rust 版改用整条 read 的编辑距离做比对，默认阈值卡在 read 长度的 3%，并保留每条 read 的所有最佳命中。双端数据按 read ID 配对，取 R1 和 R2 最佳命中的交集。举个例子：R1 同时比中 `A*01:01` 和 `HLA-H`，R2 只比中 `A*01:01`，取交集后这条片段就能给 `A*01:01` 投上一票；如果取完交集还留着 `HLA-H`，那在后续分型计分时就不会把它算进去。

---

## 只报 A、B、C，为什么库里还要放 HLA-H？

这个问题是在拿 NA11995 做测试时踩出来的。

最早我图省事，建库只收了 A、B、C 三个位点的序列。结果跑 NA11995 时，HLA-A 多报出一个 `A*68:293Q`。把那 15 对和它索引区完全匹配的 reads 拎出来回头检查，发现 **15 对全都能匹配 `HLA-H`**。库里没有 `HLA-H`，这些 reads 没地方去，自然就全投给了 `A*68:293Q`。

更有意思的是，光补上外显子注释完整的竞争序列还不够。下一轮验证又冒出了 `A*29:05`，多覆盖的 15 对 reads 里有 12 对的两端都能匹配一条 `HLA-U` 记录——而那条记录刚好缺少完整的外显子注释，之前建库时被过滤掉了。

所以现在的 format 2 索引干脆把同源的 class I 位点和假基因一并放进去参与竞争，只要有实际序列，哪怕外显子注释不全也会纳入比对。改完之后再跑 NA11995，主结果顺利回到 `A*01:01 / A*01:01`，B 和 C 位点也和 [GM11995 的细胞系分型记录](https://www.cellosaurus.org/CVCL_9600) 对上了。

参考库直接基于 [IPD-IMGT/HLA](https://github.com/ANHIG/IMGTHLA) 的发布版本构建：RNA 用拼接后的 exon 2 和 exon 3，DNA 还会带上相邻有注释的内含子。建库时会把 release 版本号、源 commit 和校验和一并记下来。以后数据库升级了，旧版本目录直接留着，想复现历史结果时 `--db` 指回去就行。

---

## 怎么跑

在 [Releases](https://github.com/pzweuj/OptiTypeRs/releases) 下载 Linux x86-64 压缩包解压，把 `optitype-rs` 放到 `PATH` 里就能直接用，Linux 和 WSL 都支持。建好参考库之后，样本分型可以完全离线跑。

下面固定用 `3.66.0` 版本数据库，方便和仓库里记录的验证结果对照：

```bash
# 1. 构建参考库
optitype-rs db build --release 3.66.0 --out refs

# 2. DNA 双端分型
optitype-rs run --dna \
  -i sample_R1.fastq.gz sample_R2.fastq.gz \
  --db refs/3.66.0 -o results --details
```

跑 RNA 数据把 `--dna` 换成 `--rna`，单端测序就只传一份 FASTQ。

结果落在输出目录的 `result.json` 里。我自己看结果一般先看 `genotype`，紧接着看 `ambiguous_loci` 和 `alternatives`。其中 `reads` 是模型覆盖的加权片段数，加上 `--details` 还能把完整的等位基因组和比对统计吐出来。

---

## 主结果有时候只是并列第一

看分型不能只盯着 `genotype` 第一行，有些样本在现有数据下本来就分不开。

拿 NA12878（从 WES BAM 的几个 HLA 区间抽出来的 5,951 对 reads）举例，OptiTypeRs 给出的 HLA-A 主结果是 `A*01:218 / A*11:01`，而 [公开参考分型（表 4）](https://pmc.ncbi.nlm.nih.gov/articles/PMC11948951/) 是 `A*01:01 / A*11:01`。第一眼看像是报偏了，但顺手看一眼 `alternatives` 就会发现，这两套组合覆盖的片段数一模一样（都是 282），目标函数分数也完全相同。手头这批 reads 在对应区域上根本无法在 `A*01:218` 和 `A*01:01` 之间做区分，主结果只是排在前面的那个并列第一。

所以我自己在下游做汇总时，会把 `ambiguous_loci` 一起带上。想多看几组并列模型的话可以加 `--enumerate`；如果结果里看到 `best_ties_truncated` 为 `true`，说明还有更多并列组合超出上限没列出来。

---

## 回归与现状

用旧版参考库和旧过滤条件做回归测试时，手头的 DNA 和 RNA 样本都保留了原版 OptiType 的两字段分型。现在换上了新的参考库、比对逻辑以及并列命中的处理方式，支持数和主叫结果自然可能出现变动，上面那两个例子主要就是用来定位到底是哪一步改了结果，具体记录都在仓库的 [Validation 说明](https://github.com/pzweuj/OptiTypeRs#validation) 里。

目前这些验证用的都是富集提取后的 HLA reads，整套全基因组 FASTQ 端到端跑下来耗时多少，后面有空再单独测一轮。
