---
title: 科学上网，shadowsocks+BBR
tags: software
---
>首先，你需要一个服务器。

推荐[vultr](https://www.vultr.com/?ref=7251924)的VPS，最低一个月2.5刀（好像卖完了），还有一个月5刀的。

然后弄好服务器之后。接下来就是装科学上网的软件。

首先是大名鼎鼎的shadowsocks（低调。。。）推荐使用大神[teddysun](https://github.com/teddysun)写的一键安装的版本。建议安装的时候选择python版本。
```
# 下载安装脚本
wget --no-check-certificate -O shadowsocks-all.sh https://raw.githubusercontent.com/teddysun/shadowsocks_install/master/shadowsocks-all.sh

# 设置权限
chmod +x shadowsocks-all.sh

# 运行
./shadowsocks-all.sh 2>&1 | tee shadowsocks-all.log
```
接下来就是按提示进行设置。
```
# 卸载
./shadowsocks-all.sh uninstall
```

然后接下来装BBR。BBR是一个加速网络的算法。
```
# 下载脚本
wget --no-check-certificate https://github.com/teddysun/across/raw/master/bbr.sh

# 设置权限
chmod +x bbr.sh

# 运行
./bbr.sh
```
服务器会重启。

然后是安装各种客户端。
[windows](https://github.com/shadowsocks/shadowsocks-windows/releases/download/4.0.6/Shadowsocks-4.0.6.zip)

[macOS X](https://github.com/shadowsocks/ShadowsocksX-NG/releases/download/v1.6.1/ShadowsocksX-NG.1.6.1.zip)

[android](https://github.com/shadowsocks/shadowsocks-android/releases/download/v4.2.5/shadowsocks-nightly-4.2.5.apk)

[iOS](https://itunes.apple.com/us/app/shadowsocks/id665729974?ls=1&mt=8)

运行的客户端，输入各种配置参数。然后更新一下PAC，就可以啦。
建议用PAC模式，这样不会代理国内网站，不会消耗太多流量。

最后感谢[clowwindy](https://github.com/clowwindy)

[^_^]:研究了很久还是觉得伸手党最简单
