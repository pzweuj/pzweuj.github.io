---
title: 自托管部署Stirling-PDF
tags: software
---

我在对PDF进行操作时，一般流程是通过Google搜索PDF处理，然后用在线的工具。譬如[ilovepdf](https://www.ilovepdf.com/zh-cn)，又譬如[SmallPDF](https://smallpdf.com/cn)。使用这些在线工具的风险是数据与第三方服务器进行交互，另外也存在只能在线使用的劣势。

[Stirling-PDF](https://www.stirlingpdf.com/)是一个开源的PDF处理工具，可以本地自部署，包含了API供程序调用。

## 部署

这里使用docker compose进行部署。由于我是部署在云服务器中，会在公网中暴露，因此必须打开登录模式。**反代啥的根据自己的实际情况自行配置**。

```yaml
services:
  stirling-pdf:
    image: stirlingtools/stirling-pdf:latest
    container_name: Stirling-PDF
    ports:
      - '8080:8080'                             # 前面的端口按需调整
    volumes:
      - ./trainingData:/usr/share/tessdata:rw   # OCR 语言支持
      - ./extraConfigs:/configs:rw              # setting文件生成的路径
      - ./customFiles:/customFiles:rw           # 可以对网页进行一些元素调整
      - ./logs:/logs:rw
    environment:
      - DOCKER_ENABLE_SECURITY=true             # 启用内部安全功能
      - SECURITY_ENABLE_LOGIN=true              # 启用登录
```

构建对应的docker-compose.yml文件后，输入下面的命令启动项目

```bash
docker compose up -d
```

在成功启动项目后，等待大概10分钟，程序才会初始化完成，然后才能通过\<ip\>:\<端口\>打开。

初始的管理员用户名和密码分别是**admin**、**stirling**，可以在网页端设置中重新设定。

### 启动失败爬坑

在使用登录模式进行部署时，Stirling-PDF会下载下面这个主程序，由于墙的存在，很可能会下载失败

官方方案是将基础镜像改为fat版本

```bash
# 例如，将下面的镜像
stirlingtools/stirling-pdf:0.44.2

# 改为
stirlingtools/stirling-pdf:0.44.2-fat
```

**也**可以使用下面的方法：

```bash
https://github.com/Stirling-Tools/Stirling-PDF/releases/download/v${version}/Stirling-PDF-with-login.jar
```

可以根据部署的版本，自行[使用镜像](https://greasyfork.org/zh-CN/scripts/412245-github-enhancement-high-speed-download)进行下载。


具体的，先确认Stirling-PDF的容器ID，进入到容器中，然后自行下载

```bash
# 确认Stirling-PDF的ID
docker ps

# 进入到容器中
docker exec -it <容器ID> /bin/bash
cd /

# 自行下载
wget https://github.com/Stirling-Tools/Stirling-PDF/releases/download/v${version}/Stirling-PDF-with-login.jar -O app-security.jar

# 将这个文件进行软链
rm app.jar
ln -s app-security.jar app.jar
exit
```

建议重启项目

```bash
docker compose restart
```

## API调用

网页端的使用没啥好说的，这里是同API调用。如果是自己用，就直接用管理员账户就好了。如果是还有其他用户，Stirling-PDF是做了用户角色分类的，这里使用一个拥有API权限的用户，在设定面板中可以看到这个用户的API密钥。

可以通过Swagger-UI查看和测试Stirling-PDF的API，具体的是访问**\<ip\>\<端口\>/swagger-ui/index.html**。

我这里以将test.docx文件转为pdf为例，查询到API结构是

```bash
curl -X 'POST' \
  'http://<ip>:<port>/api/v1/convert/file/pdf' \
  -H 'accept: */*' \
  -H 'Content-Type: multipart/form-data' \
  -H "X-API-KEY: your-api-key-here" \
  -F 'fileInput=@test.docx;type=application/vnd.openxmlformats-officedocument.wordprocessingml.document' --output test.pdf
```
