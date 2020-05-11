---
title: 随机森林
tags: coding
---

使用sklearn

```python
from sklearn.ensemble import RandomForestClassifier
import pandas as pd

# 创建树
model = RandomForestClassifier(n_estimators=100, bootstrap=True, max_depth=4)

# 导入训练集
train = pd.read_table("train.txt", sep="\t", header=0, col_index=False)

# 处理，一般需要把文本处理为数值，这里我把良性处理为0，恶性处理为1
# 切片
ytrain = train["Class"]
xtrain = train.iloc[:, 1:]

# 训练
clsf = model.fit(xtrain, ytrain)

# 导入测试组
test = pd.read_table("test.txt", sep="\t", header=0, col_index=False)

# 预测
predict_clsf = clsf.predict(test)
```

### 参考
[实战：用Python实现随机森林](https://segmentfault.com/a/1190000017320801)