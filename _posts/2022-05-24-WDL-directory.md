---
title: WDL流程中导入文件夹
tags: default
---

在集群中使用WDL流程语言，会通过singularity镜像导入环境。但是，一些巨大的数据库不可能封装进镜像中。之前在使用annovar等注释软件时，我都是不打包为镜像而是在各个节点中都安装上软件的，然后通过在WDL中写入数据库的路径来达到对应效果。

类似的task如下：
```WDL
task annovar {
	input {
		File vcf
	}
	
	String humandb = "/path/to/humandb"
	
	command <<<
		annovar ~{vcf} ~{humandb}
	>>>

	output {
		File xxx = "xxx"
	}
}
```

其实可以把整个数据库文件夹使用tar打包，然后通过**File**这个包导入镜像中。需要注意的是WDL在跨储存器时使用的是**复制**，所以最好保证这个tar文件与运行目录在同一储存器中。

```WDL
task annovar {
	input {
		File vcf
	}

	File humandb = "/path/to/humandb.tar"
	String humandb_dirname = basename(humandb, '.tar')

	command <<<
		tar -xf ~{humandb} -C ~/
		humandb_path = $(cd ~/~{humandb_dirname}; pwd)
		annovar ~{vcf} "${humandb_path}"
	>>>

	runtime {
		docker: "annovar/annovar:latest"
	}
}
```


WDL在[开发中的版本](https://github.com/openwdl/wdl/blob/main/versions/development/SPEC.md#file-directory-and-optional-outputs)有提及以后会支持**Directory**导入及输出，但是目前是尚不支持的。







