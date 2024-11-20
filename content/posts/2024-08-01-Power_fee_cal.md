---
title: 我需要转到峰谷电价？
tags: coding
---

开始柴米油盐，需要抠抠搜搜！居民用电是默认阶梯电价的，需要转换到峰谷电价需要申请，让电网来更换峰谷电表。

查了下自己的用电量，当我只保留冰箱、净水机、NAS、各种待机设备时，每天大概耗电3.5 kWh（保留设备都正常在线，然后出门溜达一天算出来的。妈啊怎么这么多）。正常的话，每天大概用电10 kWh（也不知道正不正常）。

以一个月30天计算（当然，周末的用电曲线是不同的，但是影响不会太大），我每月用电量大概是300 kWh。再经过研究自己的作息，发现我在10-12,14-19这个用电高峰根本就不在家！现在可以算一下究竟是阶梯电价划算还是峰谷电价划算了。

## 阶梯计费

以南方电网 - 广州的电价为例（注意时效，夏天和非夏天电价不同）。

|  档位  | 电价（分/kWh） | 阶梯（kWh） |
| :----: | :------------: | :---------: |
| 第一档 |     58.89      |  [0, 260]   |
| 第二档 |     63.89      | (260, 600]  |
| 第三档 |     88.89      |  (600, +∞)  |


```python
# 阶梯
def step_billing(month_total, step1=58.89, step2=63.89, step3=88.89, step2_cutoff=260, step3_cutoff=600):
    step2_amount = step3_amount = 0
    if month_total > step3_cutoff:
        step3_amount = month_total - step3_cutoff
        billing = step3_amount * step3 + (step3_cutoff - step2_cutoff) * step2 + step2_cutoff * step1
    elif month_total > step2_cutoff:
        step2_amount = month_total - step2_cutoff
        billing = step2_amount * step2 + step2_cutoff * step1
    else:
        billing = month_total * step1
    return billing
```

计算所得是178.67元。


## 峰谷计费

| 档位 | 电价（分/kWh） |     时间段      |
| :--: | :------------: | :-------------: |
| 峰电 |     99.50      |   10-12,14-19   |
| 平电 |     58.89      | 8-10,12-14,19-0 |
| 谷电 |     22.92      |       0-8       |


但是！峰谷计费仍需要叠加一个阶梯！

|  档位  | 电价（分/kWh） | 阶梯（kWh） |
| :----: | :------------: | :---------: |
| 第一档 |       0        |  [0, 200]   |
| 第二档 |      5.00      | (200, 400]  |
| 第三档 |     30.00      |  (400, +∞)  |


```python
# 阶梯
def peekValley_billing(peak, flat, valley, peak_f=99.50, flat_f=58.89, valley_f=22.92):
    billing_base = peak * peak_f + flat * flat_f + valley * valley_f
    month_total = peak + flat + valley
    billing_append = step_billing(month_total, 0, 5, 30, 200, 400)
    return billing_base + billing_append
```

我大概估算了一下我每个月这三个时间档的用电量，峰电：36kWh；平电：135kWh；谷电：129kWh。

计算所得是149.89元。

## 结论

也就是说，按照峰谷计费比阶梯计费每个月可以少给30块钱电费，可以多吃两个猪脚饭！

但是！**夏天和非夏天电价不同。**一般夏天模式是5-10月；非夏天模式是11-4月，而且非夏天模式的电费相对更高（因为阶梯差额小）。申请了改计费模式后，一年内是不能再申请的，因此，还需要综合评估一下非夏天模式的用电量。所以我决定验证够一年再决定要不要切换。

注：电车的电表独立报装，不影响原本家庭电表的阶梯。

