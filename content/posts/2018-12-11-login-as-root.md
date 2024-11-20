---
title: root登录ubuntu的图形界面
tags: default
---

今天设置ubuntu的terminal，然后发现打不开了，解决方案就是回去恢复默认设置。
但是问题是terminal都进不去，怎么恢复默认。

第一步
----
安装另外的terminal模拟器。这一步的目的是为了设置能从root进入图形界面。

第二步
----
在root用户下，修改以下，记得保存。
```bash
#设置密码
sudo passwd root
```
```bash
vi /usr/share/lightdm/lightdm.conf.d/50-unity-greeter.conf
# 增加两行
greeter-show-manual-login=true 
all-guest=false
```
```bash
vi /etc/pam.d/gdm-autologin
# 注释掉 auth required pam_succeed_if.so user != root quiet_success
```
```bash
vi /etc/pam.d/gdm-password
# 注释掉 auth required pam_succeed_if.so user != root quiet_success
```
```bash
vi /root/.profile
# 将mesg n || true修改成tty -s&&mesg n || true
```

第三步
----
注销，然后以root登录。
在文件管理中进入目录usr/share/applications。
找到gnome-terminal。右键属性，把command修改为gnome-terminal --preferences。
然后注销回到本来的用户，进入目录usr/share/applications打开gnome-terminal即可直接进入设置。
把之前导致打不开的设置修改回来，就ok了。




[-_-]:前景