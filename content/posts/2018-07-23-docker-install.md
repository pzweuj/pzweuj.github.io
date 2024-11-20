---
title: 安装docker
tags: software
---


**WSL2只需要安装新版docker软件并在软件设置中勾选WSL即可，此篇作废**



>Docker 是一个开源的应用容器引擎，让开发者可以打包他们的应用以及依赖包到一个可移植的容器中，然后发布到任何流行的 Linux 机器上，也可以实现虚拟化。容器是完全使用沙箱机制，相互之间不会有任何接口。

点击进入[docker](https://www.docker.com/)的官网。

以下是在ubuntu中全新安装docker的方式：
```bash
# set up repository
sudo apt-get update
sudo apt-get install \
	apt-transport-https \
	ca-certificates \
	curl \
	software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo apt-key fingerprint 0EBFCD88
sudo add-apt-repository \
	"deb [arch=amd64] https://download.docker.com/linux/ubuntu \
	$(lsb_release -cs) \
	stable"

# install docker ce
sudo apt-get update
sudo apt-get install docker-ce

# test
sudo docker run hello-world
```

[-_-]:昨晚梦到井