---
title: 批量下载NCCN指南
tags: coding
---

使用playwright批量下载NCCN英文指南。~~后续我想对指南建立向量知识库，然后使用DeepSeek来批量整理其中的靶向用药位点信息。~~注意，该行为违反NCCN的最终用户许可，千万不要做。下面也只是一个测试代码，我也不知道有没有用🤪。

使用下面的代码前，首先需要注册NCCN的账户。

照例，为了反反爬，会用到(stealth.min.js)[https://github.com/requireCool/stealth.min.js]

## 爬取代码

下面是爬取的python代码，修改自己的账户密码。为了避免失败，分两阶段进行，第一阶段只查询pdf的网址并保存，然后在第二阶段再进行下载。如果已获得第一阶段文件，完全可以直接进行第二阶段。

可以查看[示例代码](https://github.com/pzweuj/practice/blob/master/python/NCCN/nccn_batch_download.py)。

这里的一个坑是，playwright打开pdf url时，会默认打开为pdf viewer。需要拦截请求为下载pdf。
