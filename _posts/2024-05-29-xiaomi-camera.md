---
title: 小米摄像头储存到nas
tags: coding
---

买了一个小米云台摄像机2，在买之前问了客服是能自动将录像储存到局域网下开启了smb的nas里的。

这样可以完成我的两个需求：

1，摄像机接入米家，可随时查看；
2，至少保存1个月的录像到nas，可以追溯。

在米家中完成设置之后，发现会在选择的储存目录下创捷一个这样的结构

```
.
└── xiaomi_camera_videos
	└── 78xx7x0axxxd
		├── 2024052708
		├── 2024052709
		├── 2024052710
		└── 2024052711
			└── 58M28S_1716778708.mp4
```

这里的 **78xx7x0axxxd** 是摄像机设备的mac地址，然后在这个路径下面的日期命名目录，每个文件夹存着一个小时的录像文件。

每个mp4文件大概时间是1分钟。在我目前的清晰度和画面设置下，每个mp4文件大小约12.5MB。也就是说我想保存一个月的录像，大概需要540GB的空间。

然后这个mp4文件的命名大概是“58M28S”指的就是“2024年05月27日 11:58:28”。可以从后面这个“1716778708” unix时间戳来计算。


```python
import datetime
import pytz

# 给定的时间戳
timestamp = 1716770788

# 转换为UTC可读的日期和时间
utc_time = datetime.datetime.fromtimestamp(timestamp, pytz.utc)

# 转换为东8区时间
eastern_8_time = utc_time.astimezone(pytz.timezone('Asia/Shanghai'))

print("东8区时间:", eastern_8_time)
```

我目前在米家中设置了只保存一个月的录像，到时看看如果没有自动清理，就写个脚本挂到nas的crontab中定时清理超时的录像文件就好。
