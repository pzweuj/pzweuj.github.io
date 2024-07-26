---
title: Pandas的SettingWithCopyError定位到行
tags: coding
---

在使用pandas对数据表df进行操作时，由于对表进行了筛选过滤，但并没有去生成新表，容易触发pandas的SettingWithCopyError。

有人会设置warning把这个警告直接屏蔽掉，但是这个实际是个错误【ERROR】。

虽然表现和警告⚠️一样，流程能运行，但是实际上结果是错误的，因为原始表并没有发生改变（相当于只是在excel表里做了筛选，但是没有把筛选过滤掉的结果drop掉，没有形成新表）。

因此，这个警告实际是个错误，会导致最后输出的结果实际上是没被处理的！

一般来说，使用loc来索引，将值改变能避免这个问题
```python
df.loc[:, "col"] = "xxx"
```

然而，SettingWithCopyError并不会告诉你发生问题的位置，当代码又臭又长之后，定位困难。下面是解决方案，给代码套个try except：

```python
# 改变链式错误的回报模式
pd.options.mode.chained_assignment = 'raise'
try:
	<你的代码>
except pd.core.common.SettingWithCopyError as e:
	print(e)
```

