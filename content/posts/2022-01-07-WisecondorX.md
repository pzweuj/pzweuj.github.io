---
title: WisecondorX的使用
tags: software
---



使用[WisecondorX](https://github.com/CenterForMedicalGeneticsGhent/WisecondorX)来进行NIPT分析。

## 安装

可通过pip或者conda来安装WisecondorX，注意pip安装并不会同时安装依赖的R包。

通过pip安装
```bash
pip install -U git+https://github.com/CenterForMedicalGeneticsGhent/WisecondorX
```

通过conda安装
```bash
conda install -f -c conda-forge -c bioconda wisecondorx
```

在docker hub中找到一个可用的镜像
```bash
docker pull hindrek/wisecondorx:1.0.0
```

## 使用

wisecondorX的输入文件是bam或者cram，官方建议使用bowtie2进行比对，并且在比对前，**不要进行过滤**，因为程序需要用到低质量的reads信息。

这里使用bwa比对后，再使用sambamba进行去重（去重步骤未知是否必须）。

```bash
bwa mem -t 8 hg38.fa test.fq.gz | samtools view -bSh - | samtools sort -@ 8 - -o test.bam
sambamba test.bam test.rmdup.bam -r -p -t 8
```

使用wisecondorX将bam转换为npz格式
```bash
WisecondorX convert test.rmdup.bam test.npz
```

这样，在获得多个样本的npz文件后，将正常样本的npz文件放置于同一文件夹下建立基线。wisecondorX建议使用50~500个样本建立参考基线。
```bash
WisecondorX newref reference/*.npz reference.npz --nipt --cpus 8
```

使用基线对其他样本进行预测
```bash
WisecondorX predict test.npz reference.npz output/test --bed --plot
```

重点关注结果的test_statistics.txt文件，文件中列出了chr1-chrX等23个染色体的Z值。Z值>3的情况下是三体高风险。同时，因为加入了--plot参数运行，wisecondorX输出图片结果。





