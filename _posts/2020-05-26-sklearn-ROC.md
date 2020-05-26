---
title: sklearn ROC与AUC曲线
tags: coding
---

关于ROC与AUC曲线，[这篇文章](https://blog.csdn.net/u013385925/article/details/80385873)写的比较详细，看完基本有个大致的了解了。

这里写一写sklearn画这个曲线。

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import pandas as pd
import matplotlib.pyplot as plt

# 导入训练组测试组
train = pd.read_table("train.txt", sep="\t", header=0, index_col=False)
test = pd.read_table("test.txt", sep="\t", header=0, index_col=False)

xtrain = train.iloc[:, 1:]
ytrain = train["Class"]
xtest = test.iloc[:, 1:]
ytest = train["Class"]

# 模型
clf = RandomForestClassifier()
clf = clf.fit(xtrain, ytrain)

# 画ROC
predict_test = clf.predict(xtest)
prob_test = clf.predict_proba(xtest)
predict_test_value = prob_test[:, 1]
fpr, tpr, thresholds = roc_curve(ytest, predict_test_value)
roc_auc = auc(fpr, tpr)
plt.title("ROC Curve")
plt.plot(fpr, tpr, "b", label="AUC = %0.2f" % roc_auc)
plt.legend(loc="lower right")
plt.plot([0, 1], [0, 1], "r--")
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel("True Positive Rate")
plt.xlabel("False Positive Rate")
plt.show()
```

画出来是这样的
![ROC](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/81on227_5E_274M.png)


参考：
[AUC，ROC我看到的最透彻的讲解](https://blog.csdn.net/u013385925/article/details/80385873)

[python-sklearn中RandomForestClassifier函数以及ROC曲线绘制](https://blog.csdn.net/hjxu2016/article/details/78337308)