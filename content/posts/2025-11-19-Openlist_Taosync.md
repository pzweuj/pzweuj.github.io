---
title: 照片备份方案
tags: default
---

今天整理一下目前的照片备份方案。方案的核心目标是可以将手机的照片备份到NAS中，同时可以从NAS中**加密**备份到网盘里。

## 手机到NAS

这一步采用的是在NAS中部署[Immich](https://immich.app/)的方案，以下是我的Immich的docker compose部署配置，重点是在 my_immich.env.txt 文件中设置好 `UPLOAD_LOCATION`

```yaml
name: immich

services:
  immich-server:
    container_name: immich_server
    image: ghcr.io/immich-app/immich-server:v2.1.0
    volumes:
      - ${UPLOAD_LOCATION}:/usr/src/app/upload
      - /etc/localtime:/etc/localtime:ro
    env_file: my_immich.env.txt
    ports:
      - 2283:2283
    depends_on:
      - redis
      - database
    restart: always
    healthcheck:
      disable: false

  immich-machine-learning:
    container_name: immich_machine_learning
    image: ghcr.io/immich-app/immich-machine-learning:v2.1.0
    volumes:
      - ./model-cache:/cache
    env_file: my_immich.env.txt
    restart: always
    healthcheck:
      disable: false

  redis:
    container_name: immich_redis
    image: docker.io/valkey/valkey:8-bookworm
    healthcheck:
      test: redis-cli ping || exit 1
    restart: always

  database:
    container_name: immich_postgres
    image: ghcr.io/immich-app/postgres:14-vectorchord0.3.0-pgvectors0.2.0
    environment:
      # postgres 识别不了自定义命名的env file，软链一个.env解决
      POSTGRES_PASSWORD: ${DB_PASSWORD}
      POSTGRES_USER: ${DB_USERNAME}
      POSTGRES_DB: ${DB_DATABASE_NAME}
      POSTGRES_INITDB_ARGS: '--data-checksums'
    volumes:
      - ${DB_DATA_LOCATION}:/var/lib/postgresql/data
    restart: always

volumes:
  model-cache:
```

然后，在手机中安装Immich的客户端，配置好服务器地址，即可达成照片从上机上传到NAS，并且自动整理的这个环节。



## NAS到网盘

增加网盘，目的是本地储存炸了，还能找回数据。当前使用夸克网盘。如果不需要加密的话，直接使用NAS自带的网盘工具就可以将目录同步到网盘中。不过，考虑自己的隐私问题，我需要使用加密的方案。

### 部署Openlist

我在NAS中部署[Openlist](https://doc.oplist.org/)来挂载夸克网盘，然后在夸克网盘里创建一个空目录，再在Openlist里配置一个[Crypt类型](https://doc.oplist.org/guide/drivers/crypt)的储存。

这样，往这个储存点上传的文件，就会被加密上传到夸克网盘里（切记保存好密码）。

还需要配置一个本地储存端点，这个是给Taosync用的。



整理：

| 储存类型 | 作用                   |
| -------- | ---------------------- |
| 夸克     | 作为夸克的主要入口     |
| Crypto   | 作为夸克的加密保存入口 |
| 本地储存 | 留给Taosync的同步入口  |

下面是我的docker compose配置，核心是需要挂载上需要作为本地端点的路径（即Immich的数据路径）。

```yaml
services:
  openlist:
    image: 'openlistteam/openlist:v4.1.7'
    container_name: openlist
    user: '0:0' # Please replace `0:0` with the actual user ID and group ID you want to use to run OpenList.
    volumes:
      - './data:/opt/openlist/data'
      - '${UPLOAD_LOCATION}:/immich_lib'
    ports:
      - '5244:5244'
    environment:
      - UMASK=022
    restart: unless-stopped
```



### 部署Taosync

接下来考虑的是如何把NAS的文件同步到Openlist中。这里使用[Taosync](https://github.com/dr34m-cn/taosync)，参考[这篇文章](https://dr34m.cn/2024/07/newpost-57/)进行配置。

```yaml
version: "3.8"
services:
  taoSync:
    image: dr34m/tao-sync:latest
    container_name: taoSync
    ports:
      - "8023:8023"
    volumes:
      - ./data:/app/data
    restart: always
```

配置Taosync任务，将Openlist中的【本地储存】路径，定时**增量**备份到【Crypto】路径中。

