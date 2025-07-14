---
title: docker镜像转换为singularity sif
tags: default
---

生信集群里一般使用[singularity](https://sylabs.io/docs/)来运行容器，有一些本地现成的docker镜像，需要提前转换为singularity的sif格式文件，这样方便迁移和管理。


## 注册表拉取
适用于环境可以同时使用docker和singularity。如果是从Docker Hub (或其他 OCI 兼容的容器注册表) 中拉取，是可以直接转换为singularity镜像的。如果镜像已经存在本地docker缓存里，singularity会优先从本地docker缓存进行构建。

```bash
singularity build img.sif docker://img:tag
```

## 离线环境

适用只拥有singularity环境，但没有docker环境的情况。首先在docker环境下，将docker镜像转换为tar。

```bash
docker save -o img.tar img:tag
```

然后将tar通过scp复制到singularity环境中。

```bash
scp img.tar user@hpc.example.com:/path/to/your/directory/
```

最后在singularity环境里转换tar为sif。
```bash
singularity build img.sif docker-archive://img.tar
```
