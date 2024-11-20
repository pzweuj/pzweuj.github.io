---
title: Singularity 安装和使用
tags: software
---

### 介绍

![sylabs](https://sylabs.io/assets/svg/singularity-logo.svg)



[Singularity](https://sylabs.io/singularity/)是和docker相似的容器化软件。由于在部署Cromwell调度WDL流程时，运行docker失败。看了下文档说是集群无法支持docker，需要使用singularity。

Singularity的sif格式支持直接从docker image中转换过来，因此流程的迁移还是比较方便的。考虑到学习成本，我大概还是只会使用docker来build镜像，再用singularity来部署。在集群中使用singularity，日常普通任务时使用docker。



### 安装singularity

使用root权限安装，安装文档[参阅](https://github.com/hpcng/singularity/blob/master/INSTALL.md)。

安装依赖，我用的是centos，如果没有找到就先使用国内的源进行更新。
```bash
# yum clean all && yum makecache
yum groupinstall -y 'Development Tools'
yum install -y epel-release
yum install -y golang libseccomp-devel squashfs-tools cryptsetup
```

然后使用rpm安装，安装时由于国内网络原因，无法连接GO服务器进行校对，会报超时错误，但是貌似没有影响。
```bash
yum install -y rpm-build
export VERSION=3.8.0
wget https://github.com/hpcng/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
rpmbuild -tb singularity-${VERSION}.tar.gz
rpm -ivh ~/rpmbuild/RPMS/x86_64/singularity-${VERSION}-1.el7.x86_64.rpm
rm -rf ~/rpmbuild singularity-${VERSION}*.tar.gz
```



### 使用singularity

创建镜像
```bash
singularity -d build ubuntu.sif docker://ubuntu:latest
```



拉取镜像

```bash
singularity pull sub://ubuntu:latest
# 从docker hub中pull
singularity pull docker://ubuntu:latest
# pull并保存镜像文件
singularity pull ubuntu.latest.sif sub://ubuntu:latest
```



使用镜像

```bash
# shell方法
singularity shell ubuntu.latest.sif
# -B 参数与docker run的-v参数类似
singularity shell -B /data:/data ubuntu.latest.sif
# run方法
singularity run ubuntu.latest.sif
# exec方法
singularity exec --containall --bind ${cwd}:${docker_cwd} ubuntu.latest.sif /bin/bash script.sh
```

singularity shell的一个特点是从什么用户进去，外部权限就会是什么用户的，也就是说，在container内新建一个文件夹，对应的外部文件夹权限是属于当时创建container的外部用户的。另外，shell方法只要使用了exit，container就会退出且删除，相当于docker的run --rm。使用-B绑定多个目录时，可用","分隔。

查看运行的container，相当于docker ps -l。另外，如果使用分别用普通用户和root用户的instance start创建，在查询list时也是分开的，就是说普通用户instance list是查不到root用户的instance list的。
```bash
singularity instance list
singularity instance start <name>
singularity instance stop <name>
```

### Cromwell config

cromwell的config文件中，对于singularity配置的简单写法如下

```json
include required(classpath("application"))

backend {
	default: singularity
	providers: {
		singularity {
			actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
			config {
				run-in-background = true
				runtime-attributes = """
					String? docker
				"""
				submit-docker = """
					singularity exec --containall --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}
				"""
			}
		}
	}
}
```

合并slurm集群的配置可看[Cromwell的文档](https://cromwell.readthedocs.io/en/stable/tutorials/Containers/)。不要参考Cromwell的github仓库中的例子，那个例子有问题。

由于我将所有docker image都通过singularity转格为sif保存，我是这样改的

```json
backend {
  default = slurm

  providers {
    slurm {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = """
        Int runtime_minutes = 600
        Int cpus = 2
        Int requested_memory_mb_per_core = 8000
        String? docker
        """

        submit = """
          sbatch --wait \
            -J ${job_name} \
            -D ${cwd} \
            -o ${out} \
            -e ${err} \
            -t ${runtime_minutes} \
            ${"-c " + cpus} \
            --mem-per-cpu=${requested_memory_mb_per_core} \
            --wrap "/bin/bash ${script}"
        """
        
        submit-docker = """
          export SINGULARITY_CACHEDIR=/path/to/your/singularity_cache
          CACHE_DIR=$SINGULARITY_CACHEDIR
          DOCKER_NAME=$(sed -e 's/[^A-Za-z0-9._-]/_/g' <<< ${docker})
          IMAGE=$CACHE_DIR/$DOCKER_NAME.sif
          if [ ! -f $IMAGE ]; then
            singularity pull $IMAGE docker://${docker}
          fi
          sbatch \
            --wait \
            -J ${job_name} \
            -D ${cwd} \
            -o ${cwd}/execution/stdout \
            -e ${cwd}/execution/stderr \
            -t ${runtime_minutes} \
            ${"-c " + cpus} \
            --mem-per-cpu=${requested_memory_mb_per_core} \
            --wrap "singularity exec --containall --bind ${cwd}:${docker_cwd} $IMAGE ${job_shell} ${docker_script}"
        """

        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\d+).*"
      }
    }
  }
}
```



