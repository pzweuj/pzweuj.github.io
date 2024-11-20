---
title:  Cromwell GUI
tags: default
---

在搭建好WDL流程，以及使用[Cromwell](https://cromwell.readthedocs.io/en/stable/)来运行之后，理所当然的，就会想用GUI界面来替代命令行。现在比较有名的开源生信GUI是[Galaxy](https://galaxyproject.org/)，但是并不支持WDL和Cromwell。

在github找了一下，找到几个工具。

## diy-cromwell-server
[diy-cromwell-server](https://github.com/FredHutch/diy-cromwell-server)，最近还在更新中。需要联系作者建立账户，太麻烦，放弃。


## cromwellDashboard
[cromwellDashboard](https://github.com/cran/cromwellDashboard)，最后更新是在2018年，比较久远。

用R安装

```R
install.packages("cromwellDashboard")
```

在后台运行cromwell server
```bash
java -jar cromwell-62.jar server
```

运行
```R
library(cromwellDashboard)
runCromwellDashboard(url="127.0.0.1:8000", version="v62")
```

默认是用5412端口，界面比较简单，但是只能用于查看任务，其他所有功能都用不了。
![cdb](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/cromwellDashboard.png)


## wdl-workspace
[wdl-workspace](https://github.com/epam/wdl-workspace)，最后更新于2019年。

安装可使用npm从源码安装
```bash
npm install
npm run build
```

但是我装的时候报错，因此直接用了提供的docker。
```bash
docker pull lifescience/wdl-workspace:develop
docker run -p <LOCAL PORT>:80 -d lifescience/wdl-workspace:develop
```

这个做的比较不错，界面清爽，而且交互逻辑比较清晰，任务运行时能输出实时log。但是导入WDL脚本会有问题。无法识别import的依赖内容，导致流程不能提交（可能和使用docker版有关，没有挂载无法找到对应文件），另外没有保存流程的功能。后面还是试试从源码安装。

![wdlworkspace](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/wdlworkspace.png)


## cromwell-client
[cromwell-client](https://github.com/antonkulaga/cromwell-client)，今年还在更新中。

从docker安装
```bash
docker pull quay.io/comp-bio-aging/cromwell-web:0.3.1
```

使用
```bash
docker run -p 8001:8001 quay.io/comp-bio-aging/cromwell-web:0.3.1
```

能成功进入界面，但是在Update workflows时提示连接拒绝，连不上cromwell server。暂时用不上，后面再看看。

![cromwellweb](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/cromwellweb.png)


## cromwell-frontend
[cromwell-frontend](https://github.com/BiRG/cromwell-frontend)，最后更新于2019年。

安装
```bash
git clone https://github.com/BiRG/cromwell-frontend.git
pip install -r cromwell-frontend/cromwell_frontend/requirements.txt
```

尝试使用失败。后面再试试。




