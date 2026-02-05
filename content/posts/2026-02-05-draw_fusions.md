---
title: SchemaBio draw_fusion融合绘图
tags: software
---

改造了[Arriba](https://github.com/suhrig/arriba)的[draw_fusions](https://github.com/suhrig/arriba/blob/master/draw_fusions.R)，现在对接其他SV程序更友好了。

原代码以GPLv3开源，因此[SchemaBio/draw_fusion](https://github.com/SchemaBio/draw_fusion)也以GPLv3开源。

除了参数的更新外，还优化了代码结构，调整了巨大的if else循环，并且提高了可读性。


## 安装

建议使用docker

```bash
docker pull ghcr.io/schemabio/draw_fusion:latest
```

## 运行

现在不用以Arriba的结果格式作为输入，而是输入所有的信息

```bash
Rscript draw_fusions.R \
  --annotation=annotation.gtf.gz \    # 必填：GTF注释文件
  --output=fusion.pdf \               # 必填：输出PDF
  --gene1=TMPRSS2 \                   # 基因1名称
  --gene2=ERG \                       # 基因2名称
  --contig1=chr2 \                    # 染色体1 (chr2 或 2)
  --contig2=chr2 \                    # 染色体2
  --breakpoint1=42526250 \            # 断点1位置
  --breakpoint2=29447242              # 断点2位置
```

是的，为了提升运行效率，现在的注释文件兼容**.gtf.gz**，即使用bgzip压缩后的格式。

```bash
bgzip annotation.gtf
tabix -p gff annotation.gtf.gz
```

绘制的图和arriba的是一致的。

![plot](https://user-images.githubusercontent.com/2459065/77710742-149a9100-6f8c-11ea-8c83-3dc9a57d9379.png)


## 文档

更多的输入参数请参考[仓库页面](https://github.com/SchemaBio/draw_fusion)。