---
title: 抱抱脸部署OpenWebUI
tags: software
---

我因为自部署的[OpenWebUI](https://openwebui.com/)占用内存过高，因此转而在Vercel中部署[Lobe-Chat](https://pzweuj.github.io/posts/LobeChat)。使用了Lobe-Chat三个月，感觉还是不如OpenWebUI清爽。

当然，这个清爽是表面上的，实际OpenWebUI非常重。不过，[大善人抱抱脸](https://huggingface.co/)提供了免费的2C **16GB Ram**的免费应用实例，足够OpenWebUI用的。


## 项目部署

项目部署非常简单，在抱抱脸中，创建一个新的Space。

1. 创建新的Space(Public)
2. Space SDK选择Docker
3. 使用Blank模板即可
4. Space硬件选择免费的`CPU Basic`
5. 创建好Space后，自行在Space中增加文件`Dockerfile`

Dockerfile的内容为

```dockerfile
FROM ghcr.io/open-webui/open-webui:main
```

即可构建出OpenWebUI。

首次进入OpenWebUI的页面，会提示创建管理员账户，创建好账户后，进入到管理员设置中

将**允许新用户注册**关掉！

然后即可部署自己的模型端点和模型爽用。

## 持久化

从查询到的资料看，抱抱脸部署的项目在服务端是有数据储存的。但是，抱抱脸的免费实例存在**Sleep Mode**，即当你的 Space 长时间无人访问（通常约 45–60 分钟无请求），HF 会自动将其 休眠（sleep） 以节省资源，然后下次有人访问时，会**重新构建**并启动容器（cold start）。

为了持久化，我们需要配置一个外部数据库`DATABASE_URL`，这里可以使用[supabase](https://supabase.com/)提供的免费数据库，也可以使用其他，自便咯。


示例，在项目的Setting的Variables and secrets中，添加Secrets
```env
DATABASE_URL="postgresql://postgres.[YOUR-DBID]:[YOUR-PASSWORD]@aws-1-eu-north-1.pooler.supabase.com:6543/postgres?pgbouncer=true"
```

配置好后，可以通过更新Space中任意的文件，让它自行触发一次重新构建。

## 自定义域名

抱抱脸部署的项目，会在自己的Space下，比如我部署的OpenWebUI的域名是

```
https://huggingface.co/spaces/pzweuj/open-webui
```

同时，抱抱面会给项目生成一个域名，通过域名访问就是完整的web端，不用看到抱抱脸的框框。比如我上面这个项目，对应的域名是

```
https://pzweuj-open-webui.hf.space
```

但是我们这里讨论的是自定义域名。抱抱面的自定义域名是一个高级功能，需要开通订阅才能使用。在Cloudflare中直接用CNAME指向上面的域名是会报错的。

因此，我们需要绕一下路，使用Cloudflare的Worker路由来实现。

### Worker建立

1. 在Cloudflare左侧标签中选择“计算和AI” -> “Workers和Pages”
2. 然后在右上机点击“创建应用程序”
3. 然后使用“从Hello World！开始”模板
4. 可以给这个worker设定一个名称，比如我设定的是“hf-chat-proxy”
5. 进入到worker中，点击右上角的“编辑代码”


下面是代码，注意把我的项目域名改为自己的
```javascript
export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    // 替换为你实际的 HF Space 地址
    const target = "https://pzweuj-open-webui.hf.space";
    
    // 构造新的请求 URL
    const newUrl = new URL(request.url);
    newUrl.hostname = new URL(target).hostname;
    newUrl.protocol = "https:";

    // 转发请求（保留原始 headers）
    const modifiedRequest = new Request(newUrl, {
      method: request.method,
      headers: request.headers,
      body: request.body,
      redirect: "manual",
    });

    let response = await fetch(modifiedRequest);

    // 修改响应头，避免 CORS 或安全问题（可选）
    const newResponse = new Response(response.body, response);
    newResponse.headers.set("Access-Control-Allow-Origin", "*");

    return newResponse;
  },
};
```

6. 在设置中，把这个worker的路由设置到自己的域名里

比如

```
chat.example.com/*
```

这样就完成自定义域名了。

