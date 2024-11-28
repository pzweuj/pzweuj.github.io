---
title: GnomAD的nhemi注释
tags: coding
---

需求是注释出GnomAD的nhemi字段，但翻了一下[dbNSFP](https://usf.app.box.com/s/py6pfknr4h6464do2dw322oe2ux09hpd)和[VEP](https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html)，都没有提供这个内容，又得自己下载数据库搞了。

根据提示，nhemi即半合子计数，取GnomAD数据中的AC_XY值就好。这个值给chrX和chrY的突变提供即可。


### 数据库

UCSC里有一个[现成的注释表](https://hgw1.soe.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=varRep&hgta_track=gnomadVariantsV4.1&hgta_table=gnomadGenomesVariantsV4_1&hgta_doSchema=describe+table+schema)，版本是GnomAD v4.1的，足够新。我尝试下载了一下，发现足足有200+G，因此放弃了（如果不嫌弃的话这个版本倒是不错的，值直接了当列好了）。

由于我只需要chrX和chrY，那么就只下载[GnomAD官网](https://gnomad.broadinstitute.org/data)提供的vcf即可。

```bash
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chrX.vcf.bgz
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chrY.vcf.bgz
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chrX.vcf.bgz
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chrY.vcf.bgz
```

甲方的需求是要把exomes和genomes的值相加，所以我都下载了，然后再进行处理。

### 处理

这个文件最终的用途只是注释nhemi，因此这里将多余的字段删除，并计算exomes+genomes的总和，使用python

```python
# coding=utf-8

import gzip
import re

def extract_ac_xy(info_field):
    """从INFO字段中提取AC_XY值"""
    match = re.search(r'AC_XY=(\d+)', info_field)
    return int(match.group(1)) if match else 0

def process_vcf_files(exome_file, genome_file, output_file, chrom_name="chrY"):
    # 存储位置和AC_XY值的字典
    variants = {}
    
    # 处理exome文件
    with gzip.open(exome_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            pos_key = f"{fields[0]}_{fields[1]}_{fields[3]}_{fields[4]}"  # CHROM_POS_REF_ALT
            ac_xy = extract_ac_xy(fields[7])
            variants[pos_key] = {'exome': ac_xy, 'genome': 0, 'pos': fields[1], 
                               'ref': fields[3], 'alt': fields[4]}

    # 处理genome文件
    with gzip.open(genome_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            pos_key = f"{fields[0]}_{fields[1]}_{fields[3]}_{fields[4]}"
            ac_xy = extract_ac_xy(fields[7])
            if pos_key in variants:
                variants[pos_key]['genome'] = ac_xy
            else:
                variants[pos_key] = {'exome': 0, 'genome': ac_xy, 'pos': fields[1], 
                                   'ref': fields[3], 'alt': fields[4]}

    # 写入输出文件
    with open(output_file, 'w') as out:
        # 写入头部信息
        out.write('##fileformat=VCFv4.2\n')
        out.write(f'##INFO=<ID=Exomes_AC_XY,Number=1,Type=Integer,Description="Allele count for {chrom_name} in exomes">\n')
        out.write(f'##INFO=<ID=Genomes_AC_XY,Number=1,Type=Integer,Description="Allele count for {chrom_name} in genomes">\n')
        out.write(f'##INFO=<ID=Gnomad_AC_XY,Number=1,Type=Integer,Description="Combined allele count for {chrom_name}">\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        # 写入变异信息
        for pos_key, data in sorted(variants.items()):
            chrom = chrom_name
            total_ac = data['exome'] + data['genome']
            info = f'Exomes_AC_XY={data["exome"]};Genomes_AC_XY={data["genome"]};Gnomad_AC_XY={total_ac}'
            out.write(f'{chrom}\t{data["pos"]}\t.\t{data["ref"]}\t{data["alt"]}\t.\t.\t{info}\n')

# chrY
process_vcf_files('gnomad.exomes.v4.1.sites.chrY.vcf.bgz', 'gnomad.genomes.v4.1.sites.chrY.vcf.bgz', 'gnomad.v4.1.sites.chrY.combined_AC_XY.vcf')

# chrX
process_vcf_files('gnomad.exomes.v4.1.sites.chrX.vcf.bgz', 'gnomad.genomes.v4.1.sites.chrX.vcf.bgz', 'gnomad.v4.1.sites.chrX.combined_AC_XY.vcf', "chrX")
```

chrX处理完有一点顺序问题，使用bedtools再处理一下，不知道chrY有没有问题，也顺便处理一下

```bash
bedtools sort -i gnomad.v4.1.sites.chrX.combined_AC_XY.vcf -header > gnomad.v4.1.sites.chrX.combined_AC_XY.sort.vcf
bedtools sort -i gnomad.v4.1.sites.chrY.combined_AC_XY.vcf -header > gnomad.v4.1.sites.chrY.combined_AC_XY.sort.vcf
```

### liftover

是的，还是在用GRCh37/hg19的，需要liftover一下，不去管是否会造成的数据损失了。首先下载hg38转hg19的chain file

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
gunzip hg38ToHg19.over.chain.gz
```

为了一步到位，先把chrX和chrY合起来

```bash
cp gnomad.v4.1.sites.chrX.combined_AC_XY.sort.vcf gnomad.v4.1.sites.combined_AC_XY.vcf
awk 'NR>1 && !/^#/' gnomad.v4.1.sites.chrY.combined_AC_XY.sort.vcf >> gnomad.v4.1.sites.combined_AC_XY.vcf
```

然后可以使用gatk/picard对vcf进行liftover

```bash
gatk LiftoverVcf \
    -I gnomad.v4.1.sites.combined_AC_XY.vcf \
    -O gnomad.v4.1.sites.combined_AC_XY.hg19.vcf \
    -C hg38ToHg19.over.chain \
    --REJECT rejected_variants.vcf \
    -R hg19.fa
```


最后使用bgzip压缩，建立索引（用gatk可能可以直接输出为.vcf.gz？我没试。）

```bash
bgzip gnomad.v4.1.sites.combined_AC_XY.hg19.vcf
tabix -p vcf gnomad.v4.1.sites.combined_AC_XY.hg19.vcf.gz
```



### 注释

还是使用VEP进行注释

```bash
vep \
    -i input.vcf -o output.vcf \
    --custom file=gnomad.v4.1.sites.combined_AC_XY.hg19.vcf.gz,short_name=GnomHemi,format=vcf,type=exact,coords=0,fields=Exomes_AC_XY%Genomes_AC_XY%Gnomad_AC_XY
```



