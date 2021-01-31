---
title: oculus quest2激活
tags: default
---

为了玩[Half-Life: Alyx](https://store.steampowered.com/app/546560/HalfLife_Alyx/)，在美亚买了一台oculus quest2，目前装备未到货，但是在之前看到很多说激活麻烦的，先来预习一下。
看了几篇文章，激活设备的难点主要在于魔法上网，那么其实对我来说就没有难点。主要是看用什么样的方式连进魔法WiFi。

###  连魔法路由
我本来是有个华硕的ac68u刷了梅林的，可以轻松在路由器层面实现魔法，但是后面升级到WiFi6的时候换成了小米的ax1800，把ac68u卖掉了。然后魔法的事情留给了PC来做，毕竟家里除了我也没有别的Google重度依赖者。

然后我查了一下，发现ax1800也是可以开启ssh然后装clash来实现魔法的。首先参考[这篇文章](https://www.right.com.cn/forum/forum.php?mod=viewthread&tid=4032490&extra=page%3D1%26filter%3Dtypeid%26typeid%3D44)，开启ax1800的ssh。然后使用再ssh进去装一个shell下的[clash](https://github.com/juewuy/ShellClash)。再自己找个机场就ok。

然而问题又来了，我的ax1800固件已经自动升级到了1.0.378，比上面可以开启ssh的版本号都要新，我估计是已经把这个漏洞堵上了（纯猜测，没试验）。这时可以考虑把路由固件降级，再开启ssh。但是，我毕竟只是需要让oculus quest2用上魔法，并不是需要重新拥有一个魔法路由，这样做反而**把简单问题复杂化**了。因此直接否决以上整个路线。

### 连魔法手机
想了想其实可以通过手机先架好梯子，然后再开热点给oculus quest2，实现同样的效果。问题是我用其他设备（switch），发现直接连上还是不能到外面世界的，猜测应该是需要设置代理。这时候的问题就转变为要找到手机的ip地址。然后想了想其实直接用电脑不就好了，前面提到过[我的电脑](https://pzweuj.github.io/2020/11/12/ax200.html)其实也是连的网线，同时也有个WiFi6网卡。

### 连电脑
我的电脑现在是通过clash实现魔法的，把clash设置里的Allow Lan打开，同时Clash上会显示当前的网卡的内网ip，都不用自己去cmd用ipconfig查了。然后win10自带了热点开关，打开开关再把手机连上，同时再在手机wifi连接里设置手动代理，主机名选电脑的ip，端口填clash设置的端口，就可以了。测试也成功了，到时oculus quest2来了直接按这个操作也就可以了，我就不信它的wifi连接不能设代理。

