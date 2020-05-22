---
title: sklearn 鸢尾花
tags: coding
---

[上次](https://pzweuj.github.io/2020/03/05/weka.html)使用weka对鸢尾花数据集进行了一个预测，这次使用sklearn试一下。其实sklearn貌似内置了iris的数据，可以直接用，但是为了直观一点，还是自己去下载数据。然后为了方便之后带到其他情况中，分成两个数据来处理。

照旧，数据集在[这里](https://archive.ics.uci.edu/ml/datasets/Iris)可以找到。

这里处理成一个比较方便的格式。
格式如下：
```
Class	sepal_length	sepal_width	petal_length	petal_width
Iris-setosa	5.1	3.5	1.4	0.2
Iris-setosa	4.9	3	1.4	0.2
Iris-setosa	5	3.6	1.4	0.2
Iris-setosa	5.4	3.9	1.7	0.4
Iris-setosa	4.6	3.4	1.4	0.3
Iris-setosa	5	3.4	1.5	0.2
Iris-setosa	4.4	2.9	1.4	0.2
Iris-setosa	5.4	3.7	1.5	0.2
Iris-setosa	4.8	3.4	1.6	0.2
Iris-setosa	4.3	3	1.1	0.1
Iris-setosa	5.7	4.4	1.5	0.4
Iris-setosa	5.4	3.9	1.3	0.4
Iris-setosa	5.1	3.5	1.4	0.3
```
随意弄一些作为测试组。

---------------

试一试用ComplementNB算法来训练
```python
from sklearn.naive_bayes import ComplementNB
import pandas as pd
import joblib

# 处理成数字
def normalization(arff):
	df = pd.read_table(arff, sep="\t", header=0, index_col=False)
	for col in df.columns.values.tolist():
		if col == "Class":
			for i in df["Class"].index:
				if df.loc[i, "Class"] == "Iris-setosa":
					df.loc[i, "Class"] = 0
				elif df.loc[i, "Class"] == "Iris-versicolor":
					df.loc[i, "Class"] = 1
				else:
					df.loc[i, "Class"] = 2
	return df
	
model = ComplementNB()
train = normalization("train.txt")
xtrain = train.iloc[:, 1:]
ytrain = train["Class"]
ytrain = ytrain.astype("int")

clf = model.fit(xtrain, ytrain)

test = normalization("test.txt")
xtest = test.iloc[:, 1:]
predict_clf = list(clf.predict(xtest))

```

结果并没有特别好，但是流程是走通了。多换几个模型试一试就可以了。