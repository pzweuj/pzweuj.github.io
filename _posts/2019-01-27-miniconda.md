---
title: 这才是conda的正确姿势
tags: default
---


Conda 是一个开源的软件包管理系统和环境管理系统，用于安装多个版本的软件包及其依赖关系，并在它们之间轻松切换。


曾几何时，我一度沉迷conda，因为用conda安装各种软件实在是太方便了。不过，后来是因为conda污染了我的环境，所以被我弃用。现在算是为conda正一下名，因为之前用的姿势不太对。

### 安装miniconda
做生信并不需要装anaconda，装miniconda就行了。

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

### 创建环境
这一步很重要！就是最关键的，把每个流程创建一个环境，这样不同流程就换不同的环境。

如创建一个wes环境
```bash
conda create -n wes <然后就是需要的软件，也可以激活环境之后再往里面装软件>

# 激活环境
source activate wes
# 退出环境
source deactivate
```

我觉得这样的缺点是，很多软件和环境实际上已经搭载过，但是还是要重新搭建。
所以其实我个人还是比较喜欢直接的去找到软件来下载而不是通过conda安装的。

### 附上conda添加清华镜像的方法

更新：
清华tuna镜像已经恢复。

删除过去设置的源，恢复默认：

第一种方法：修改~目录下的.condarc。

第二种方法
```bash
conda config --remove-key channels
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda
```

设置为清华源
```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/
```


[-_-]:继续努力