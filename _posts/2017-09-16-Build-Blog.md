---
title: github建立博客！
tags: default
---

> 事情是这样的，前几天，我想着利用github page建一个博客。由于零基础，我google了很多教程。~~但是！都不好用！因为他们都过时了！~~我发现，现在用github建博客非常简单！根本不像那些教程做的那么复杂！当然，如果你要很炫很折腾的那种，就去看那些教程。下面来说说我的教程。

第一步，因为是github博客，一个github账户是必须有的。点击[github](https://github.com/)进去，注册一个账户。有了账户之后，点击如下位置，新建一个仓库（repository）。
![new-repository](https://github.com/pzweuj/pzweuj.github.io/raw/master/downloads/images/github-new-repository.png)

第二步，按下面新建一个仓库，仓库名是username.github.io。注意username改成自己的名字！
![create-repository](https://github.com/pzweuj/pzweuj.github.io/raw/master/downloads/images/github-create-new-repository.png)

第三步，安装[github desktop](https://desktop.github.com/)，安装过程略，装好之后登录你的github账户，ctrl+shift+O，然后选择你的username.github.io，clone到本地。

第四步，到[jekyll themes](http://jekyllthemes.org/)挑一个自己喜欢的模板。（在这里顺便感谢[jekyll](https://github.com/jekyll/jekyll)，写了个那么好的轮子！）选到自己的心水模板，点Homepage可以跳转到他的github仓库，这时候，点击右上角那个绿绿的**Clone or download**，然后点击download zip。之后，把这个压缩包解压到你本地的username.github.io文件夹里，就可以了。在这里推荐一下我正在用的[_SimpleGray_](https://github.com/mytharcher/SimpleGray)。

第五步，修改_config.yml这个文件，把内容改成你自己的东西。然后往_posts文件夹里扔你的文章，利用github desktop同步到github仓库，就可以发表了。注意文章的标题必须是：
```
yyyy-mm-dd-title.md
```
[^_^]:真正意义上的第一篇文章，献给我的第一个读者$#。