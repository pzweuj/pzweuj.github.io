---
title: 修改并更新singularity的sif镜像
tags: default
---

我想修改singularity的镜像sif文件中的一个脚本，然后形成新的镜像sif。问了下GPT，可以使用下面的操作。

### 构建沙盒


首先，需要从现有的.sif文件中提取内容。可以使用Singularity的singularity build命令来将.sif文件转换成可修改的目录或Sandbox。

```bash
singularity build --sandbox /path/to/sandbox old.sif
```

### 修改代码

进入到沙盒中，按需修改代码及文件。

```bash
cd /path/to/sandbox
vi script.sh
```

### 测试沙盒

可以通过shell进入沙盒进行测试，和进入sif是一样的

```bash
singularity shell --writable /path/to/sandbox
```

### 生成sif

将修改后的沙盒打包为新的sif，还有删掉临时文件

```bash
singularity build new.sif /path/to/sandbox
rm -rf /path/to/sandbox
```





