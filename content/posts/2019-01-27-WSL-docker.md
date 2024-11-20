---
title: WSL怎么用docker
tags: default
---

[以前写的docker安装流程](https://pzweuj.github.io/2018/07/23/docker-install.html)


WSL里装了docker，但是运行的时候会报
```shell
docker: Cannot connect to the Docker daemon at unix:///var/run/docker.sock. Is the docker daemon running?.                See 'docker run --help'. 
```

解决办法：

运行

```shell
docker -H localhost:2375 images
```

可以把这句添加到环境中
```shell
export DOCKER_HOST=localhost:2375
```

but! 这个时候还是不能用的，是因为这还是在windows下，需要装windows的客户端。
对右下角图标右键点setting。然后勾上这个，就可以啦！


![wsl_docker](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/wsl-docker.jpg)





[-_-]:继续努力