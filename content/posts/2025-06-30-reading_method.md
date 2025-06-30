---
title: 目前的自托管看书方案
tags: default
---

先说明一下这个看书方案的需求：

1. 自托管；
2. 跨端同步看书进度；
3. 支持epub、txt；
4. 外网访问。

其实这个需求挺简单的，但是一直没有找到比较合适的客户端软件。在GitHub上能找到很多开源的阅读器项目，但是他们一般都将看书进度作为一个付费功能。

我现在的方案是：

[Talebook](https://github.com/talebook/talebook) + Webdav + [阅读](https://github.com/gedoor/legado)。


## Talebook

首先，Talebook（基于Calibre）用于管理所有的电子书，我在NAS中进行部署，所有的电子书都存在自己的本地。电子书的资源可以在[Z-library](https://z-library.cc/)中找。这个服务只是对电子书进行管理、刮削等，其实在这个方案里可有可无。

我通过docker-compose部署

```yaml
services:

  # main service
  talebook:
    restart: always
    image: talebook/talebook
    volumes:
      - ./data:/data
    ports:
      - "<your_port>:80"
      # - "8443:443"
    environment:
      - PUID=1000
      - PGID=1000
      - TZ=Asia/Shanghai
    depends_on:
      - douban-rs-api

  # optional, for meta plugins
  # please set "http://douban-rs-api" in settings
  douban-rs-api:
    restart: always
    image: ghcr.io/cxfksword/douban-api-rs
```

## Webdav

然后，本来我是想通过NAS自行打开一个webdav服务的，但是想了一下，这样做还需要暴露公网，而且单独为了一个看书服务而打开webdav，不太好。因此，我选择了更奇葩的路线。

我在NAS中计划任务定时同步本地的电子书到OneDrive中，然后把OneDrive挂载到[Koofr](https://koofr.eu/)里。Koofr目前提供10GB免费的网盘，而且可以直连，同时可以将挂载到其中的OneDrive通过Webdav共享出来。

这样，NAS也无需暴露于公网。

--------------------

如果NAS已经暴露了Webdav服务在公网中（比如通过[OpenList](https://github.com/OpenListTeam/OpenList)），就无需绕上面的弯路。


## 阅读

阅读这个软件我之前的使用模式一般都是导入书源然后搜索、放入书架、看书这个模式。之前我想过能不能将Talebook的服务暴露于公网然后形成一个阅读的书源，但是没有找到对应的解决方案。

因此，现在使用的是阅读APP里提供的增加远程书籍功能，这个功能可以读取Webdav中的文件。

同时，阅读也带有通过Webdav来备份进度的设置选项，这样一来，我的看书方案就走通了。



