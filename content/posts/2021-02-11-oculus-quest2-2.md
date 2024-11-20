---
title: oculus quest2 流程与感想
tags: default
---

[前情提要](https://pzweuj.github.io/2021/01/31/oculus-quest2.html)

## 到货
从美亚入手了oculus quest2含税包邮大概是2232元人民币，速度还挺快的，1月31日下单，2月5日就到了。到手后就着手激活，但是之前考虑的激活流程无法激活。查了下大概是clash还是走不了UDP通道还是啥来着。原以为很轻松，然而我在这里居然卡了一个小时（为了找个合适的方法）。

最后我用了[Netch](https://github.com/NetchX/Netch)来进行魔法。和clash不同，Netch是创建了一个虚拟网卡，通过虚拟网卡来实现魔法（必须安装TAP-windows）。按照Github页面上的这篇[教程](https://github.com/NetchX/Netch/blob/master/docs/Quickstart.zh-CN.md)，就能安装完成。完成后还是同样的，把自己的机场填进去，然后最重要的一步是代理必须选择TUN/TAP那个模式。然后打开win10的热点，再在网络适配器中选择tap的网卡，右键属性然后共享，下拉选择热点的网卡。

这时在quest2里面就可以成功连接啦。

## 开发者
接下来，需要在手机中装上oculus app。然后登录自己的FB账户再绑定设备。同时，为了后面便于安装应用，需要打开开发者模式。先到[这个网站](https://developer.oculus.com/)进行登录自己的FB账户，然后就可以创建一个属于自己的开发者团队了。后面为了将设备与电脑连接，还需要安装这个[USB驱动](https://developer.oculus.com/downloads/package/oculus-adb-drivers/)。接下来就是在手机的oculus app中，打开设备的开发者模式。

## SideQuest
sideQuest实际上是相当于往电脑上装的手机应用市场之类的软件，可以通过sideQuest往设备中安装第三方应用。具体安装方法就参考[sideQuest](https://sidequestvr.com/setup-howto)官网好了。就这样，设置好了之后，我们自己找一些安卓apk，比如clash，在设备连接到电脑时，把apk拖到sideQuest中，就直接安装了。这时就可以在设备自己里面实现魔法，没必要用电脑热点了。

## 无线串流
买一体机肯定是为了无线串流，不可能用有线的。现在最好的串流应用只能是[virtual desktop](https://www.oculus.com/experiences/quest/2017050365004772/?locale=zh_HK)。建议20刀直接买好了，效果比免费的ALVR好的不是一点半点。这个可以通过paypal或者在某鱼中找人代购，会给你发兑换码，然后在手机的oculus app里兑换就好了。同时在设备和电脑中安装virtual desktop，只要两者是在同一个网段里，就能直接找到。这里有个小bug，就是如果在设备中开了魔法，需要在魔法应用里设置virtual desktop不走魔法通道，不然可能会无法连接。

另外，电脑端的oculus app也要装上，可能是有些驱动或者依赖项的原因，一开始我没装这个，就没法把steam里的半衰期爱丽丝串流到quest2上，每次打开都报错。

## 感想
VR真提供了一种异常沉浸的感觉。本来我以为会很晕，实际上，因为移动是通过传送的方式，所以并没有晕。佩戴时，需要反复调整头盔找一个最清晰的角度再进行后面的操作，不然会很累的。

玩半衰期，在序章那个突然一个巨大的铺网线（？）的机器人从头顶出现，真的很震撼，后面打怪物也十分吓人（所以我把修改器的无敌打开了，至少有了安全感）。

virtual desktop就这个软件打开看到一片星空，也很震撼。归根结底，就是**身临其境**四个字。