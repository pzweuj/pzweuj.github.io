---
title: 华法林剂量计算器
tags: coding
---

很久以前写了一个基于IWPC公式的[华法林剂量计算函数](https://pzweuj.github.io/posts/IWPC)，现在基于此，实现一个线上的计算器项目。


## 信息收集

为了让项目内容更丰度，增加更多的计算模型。前期搜集信息基于此篇文章：

[模型引导的华法林精准用药：中国专家共识（2022 版）](https://xadxyylib.yuntsg.com/ueditor/jsp/upload/file/20240322/1711076013715062015.pdf)

文章中提及了五个模型的计算公式：

1. IWPC模型
2. Gage模型
3. 湘雅模型
4. 苜蓿草模型
5. Biss模型

我已将上述模型的计算方法都实现到了项目中。

但是，文章中IWPC模型的公式存在问题，遗漏了平方计算（IWPC文章中对此公式的描述是该得分等于周剂量的平方根）。

## 项目使用

![warfarin_cal](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/warfarin_cal.png)



项目基于NextJs开发，已部署到Vercel。可以通过下面的链接在线使用。

[warfarin-dosage-calculator.vercel.app](https://warfarin-dosage-calculator.vercel.app/)



如无法访问，也可以使用Tauri封装的本地版本。

[warfarin-dosage-calculator_0.0.3_x64-setup.exe](https://github.com/pzweuj/Warfarin-Dosage-Calculator/releases/download/v0.0.3/warfarin-dosage-calculator_0.0.3_x64-setup.exe)







