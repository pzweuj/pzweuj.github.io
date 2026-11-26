---
title: HaploGrouperRS，用 Rust 重写单倍群分型工具 HaploGrouper
tags: software
---

Y 染色体和线粒体（mtDNA）单倍群分型，最经典的思路之一就是拿着一条已知的系统发育树和各分支的突变位点表，一路从根节点往下累加支持和反对的突变来打分。Rootara之前在用的 [HaploGrouper](https://gitlab.com/bio_anth_decode/haploGrouper) 就是这套路子，只不过原版是 Python 写的。

在把 [Rootara](/posts/Rootara_v2) 后端改成全 Rust 的过程中，为了把 Y 染色体和线粒体单倍群判定也一并塞进单二进制工具链里，我把 HaploGrouper 用 Rust 重写了一遍，开源在 [HaploGrouperRS](https://github.com/pzweuj/HaploGrouperRS)。目前已经更新到 `0.2.0`。

---

## 直接替换原版命令

HaploGrouperRS 支持读取明文 VCF、gzip 和 BGZF 压缩的 VCF。日常使用基本不需要重新学参数，把原来的 `python haploGrouper.py` 直接替换成 `haplogrouper-rs`（Windows 下是 `haplogrouper-rs.exe`）即可，原有的长短选项、目录和相对路径习惯都保留了下来：

```bash
haplogrouper-rs -v sample.vcf.gz \
  -t data/chrY_isogg2019_tree.txt \
  -l data/chrY_isogg2019-decode1_loci_b37.txt \
  -c Y -m -o calls.tsv
```

- **Y 染色体参考数据**：仓库里直接带了基于 ISOGG 2019（含 deCODE 修订位点）的 GRCh37 和 GRCh38 参考表（详情见 [Y-DATA.md](https://github.com/pzweuj/HaploGrouperRS/blob/main/data/Y-DATA.md)）。跑 GRCh38 数据时把位点表换成 `chrY_isogg2019-decode1_loci_b38.txt` 就行，多染色体 VCF 记得通过 `-c` 指定精确的染色体名（比如 `Y` 或 `chrY`）。
- **线粒体与自定义树**：接口完全通用，做 mtDNA 时把 `-t` 和 `-l` 指向上游的 PhyloTree Build 17 参考表就行。

---

## 打分规则和细节修复

分型判定沿用了原版的单倍型惯例：遇到包含多个等位基因的 `GT` 字段时，只取第一个等位基因；随后从系统发育树根节点出发，沿着每条路径累加支持突变数（support）和反对突变数（opposition），净得分为 `support - opposition`。

最高分的候选必须在末端分支上拿到实际支持；如果出现得分并列，先按反对数最少排序，再退回到最近共同祖先（MRCA）收束；要是没有任何候选满足条件，就直接报根节点。

拿上游自带的 HG00096 测试数据对过一遍结果：mtDNA 判定为 `H16a1`（净得分 45），deCODE GRCh37 Y 染色体参考判定为 `R1b1a1b1a1a2c1a1f1`（净得分 952），和原版完全一致。另外在实现上顺手处理了 VCF 末尾 indel 拖尾碱基的兼容问题，遇到非法输入也会老老实实返回非零退出码。

---

## `0.2.0` 的内存压缩

`0.2.0` 版本专门针对大样本 VCF 和得分矩阵输出（`-x`）重写了存储结构：

- 所有样本基因型一致的位点只存一个等位基因；
- 选中样本数大于等于 16 的二等位位点，直接压缩成每个样本 1 bit；
- 用 `-f` 做参考序列推断时，按行流式扫描 FASTA；
- 输出 `-x` 大矩阵时只开辟一个 4 MiB 的打分缓冲，溢出部分借用系统临时文件流式写出（大概需要 `16 bytes × 输出节点数 × 样本数` 的临时磁盘空间）；如果只跑分型不输出矩阵，连临时文件都不会产生。

在 Windows x64 release 构建下用上游真实队列 VCF 做了一组实测对比，四组运行输出文件的 SHA-256 校验与前一版完全字节级一致，内存峰值直接压下来一大截：

| 数据集 | 样本数 | 输出模式 | 早期 Rust 版 | `0.2.0` |
| --- | ---: | --- | ---: | ---: |
| GRCh37 Y, ISOGG 2019 deCODE | 1,233 | 仅分型结果 | 50.36 MiB | **29.12 MiB** |
| GRCh37 Y, ISOGG 2019 deCODE | 1,233 | 分型 + 得分矩阵 | 119.63 MiB | **32.40 MiB** |
| Phase 3 mtDNA, Build 17 (rCRS) | 2,534 | 仅分型结果 | 16.87 MiB | **10.88 MiB** |
| Phase 3 mtDNA, Build 17 (rCRS) | 2,534 | 分型 + 得分矩阵 | 227.02 MiB | **14.96 MiB** |

尤其是两千多样本跑 mtDNA 加得分矩阵那一项，峰值直接从 227 MiB 掉到了不到 15 MiB。

---

## 下载与许可说明

Windows、Linux x64 和 macOS 的预编译二进制都可以直接从 [GitHub Releases](https://github.com/pzweuj/HaploGrouperRS/releases) 下载，ISOGG Y 染色体参考数据单独打在 `haplogrouper-rs-y-data.tar.gz` 里，解压即可。源码编译直接 `cargo build --release --locked` 就行。

使用前留意一下许可边界：HaploGrouperRS `0.2.0` 新增的代码采用限制商业用途的源码可用协议（上游 HaploGrouper 及早期贡献保持 Apache-2.0），随附的 ISOGG 数据遵循 `CC BY-NC-SA 3.0`。如果是发文章用，记得引用 Jagadeesan 等人的原版论文（*Bioinformatics* 37, 570–572）。
