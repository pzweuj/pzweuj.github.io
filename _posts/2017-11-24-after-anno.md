---
title: 关于注释之后怎么进行基本的筛选
tags: default
---
>因为最近很冷！不想码字，所以随便更新一点点

用annovar注释出来的文件。
我们可以放进excel里面看。
首先应该去关注一下clinvar有没有注释出致病（pathogenic）的位点。
但是！因为clinvar不是很准确，所以我们只能用来作为一个参考的标准。

第二步，筛选一下1000g_all的突变频率。一般以0.001也就是0.1%为准。
因为突变率太高的话，就说明这个突变在人群中是常见的，并不是罕见的变异，没有参考的价值。
同时，可以筛选EXac_eas的频率（表示东亚人），当然其他区域的人筛选其他的。

第三步，去除同义突变，我们要的是没有研究过的以及非同义突变。这样才有意义。

第四步，对剩下的进行与临床表型的匹配。这时候可以借助一些软件。
比如Exomiser。也可以用一些数据库帮助判断。比如国内的奇恩生物。

最后，筛选剩下大概五六个位点，可以去人工判断。
借助IGV直接去查看测序的相关区域。

得到最最准确的位点。
拿去一代验证。

[T_T]:我的电源快来！