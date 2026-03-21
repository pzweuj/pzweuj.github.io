---
title: TransVar2API - HGVS变异注释工具的Web服务化
tags: software
---

[TransVar](https://github.com/zwdzwd/transvar)是一个强大的HGVS变异注释工具，可以在不同坐标系统之间进行转换和注释。虽然TransVar本身提供了在线版本，但是我希望能接入我的[BioTools](https://use.biotools.space)，于是我开发了[TransVar2API](https://github.com/pzweuj/TransVar2API)，将其封装为Web服务和API。

## HGVS与TransVar

HGVS（Human Genome Variation Society）命名规范是描述遗传变异的标准。比如 `PIK3CA:p.E545K` 表示PIK3CA基因第545位谷氨酸被赖氨酸替代。这种命名直观，但在进行实际分析时，需要转换为基因组坐标。

TransVar正是解决这个问题的工具，支持：
- **panno**: 蛋白水平注释
- **canno**: cDNA水平注释
- **ganno**: 基因组水平注释

## TransVar2API

我把TransVar封装成了FastAPI服务，主要特点：

1. **多数据库支持**: RefSeq、Ensembl、GENCODE、UCSC、CCDS
2. **双版本支持**: hg38 (GRCh38) 和 hg19 (GRCh37)
3. **Web界面**: 简洁的界面，支持单个和批量注释
4. **RESTful API**: 方便集成到自动化流程
5. **Docker部署**: 一键部署，数据库内置

### 部署

```bash
# 克隆项目
git clone https://github.com/pzweuj/TransVar2API.git
cd TransVar2API

# 构建并启动
docker-compose build
docker-compose up -d
```

构建时间约20-30分钟，包含参考基因组和注释数据库的下载。启动后访问 `http://localhost:8000` 即可使用。

### API使用

单个变异注释：

```bash
curl -X POST http://localhost:8000/api/annotate \
  -H "Content-Type: application/json" \
  -d '{
    "variant": "PIK3CA:p.E545K",
    "refversion": "hg38",
    "mode": "panno",
    "databases": ["refseq", "ensembl"]
  }'
```

批量注释：

```bash
curl -X POST http://localhost:8000/api/batch_annotate \
  -H "Content-Type: application/json" \
  -d '{
    "variants": ["PIK3CA:p.E545K", "EGFR:p.L858R"],
    "refversion": "hg38",
    "mode": "panno",
    "databases": ["refseq"]
  }'
```

### 部署到HuggingFace Spaces

本项目支持一键部署到HuggingFace Spaces：

1. Fork项目到你的GitHub
2. 在HuggingFace创建新的Space (Docker类型)
3. 连接GitHub仓库即可

由于HuggingFace Spaces提供免费的容器运行环境，非常适合作为公开的注释服务。

这是我部署的[在线版本](https://pzweuj-transvarweb.hf.space/)的[API参考文档](https://pzweuj-transvarweb.hf.space/docs)。

当然，你也可以在[BioTools](https://use.biotools.space/tools/transvar)中获得一站式的多工具体验。

## 常见变异格式

| 模式 | 格式示例 | 说明 |
|------|----------|------|
| panno | `PIK3CA:p.E545K` | 蛋白变异 |
| canno | `NM_006218.4:c.1633G>A` | cDNA变异 |
| ganno | `chr3:g.178921852G>A` | 基因组变异 |

