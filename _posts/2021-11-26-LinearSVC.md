---
title: SKlearn LinearSVC
tags: coding
---



记录一些用得上的。



### 数据初始化

使用pandas对数据进行处理，将分类转换为数字。

数据表样式是

| Sample | Class      | Feature1 | Feature2 | Feature3 | Feature4 |
| ------ | ---------- | -------- | -------- | -------- | -------- |
| S1     | Pathogenic | 0        | 2        | 8        | 1        |
| S2     | Benign     | 6        | 2        | 5        | 0        |
| S3     | Other      | 9        | 3        | 7        | 1        |

计算log2CPM值，但为了避免0值使用的加1方法可能会[导致倍数不准](https://support.bioconductor.org/p/107719/)。


```python
import pandas as pd
from math import log

def normalizationDF(table, clsList):
	df = pd.read_table(table, sep="\t", header=0, index_col=0)
	for i in df.index:
		types = df.loc[i, "Class"]
		for c in range(len(clsList)):
			df.loc[i, "Class"] = c
	dropRow = []
	for i in df.index:
		if not isinstance(df.loc[i, "Class"], int):
			dropRow.append(i)
	df = df.drop(dropRow)
	colList = df.columns.values.tolist()
	colSumDict = {}
	for i in df.index:
		gSum = 0
		for col in colList:
			if col != "Class":
				gSum += df.loc[i, col]
		colSumDict[i] = gSum
	for i in df.index:
		for col in colList:
			if col !=  "Class":
				# 使用log2CPM
				df.loc[i, col] = log((df.loc[i, col] / colSumDict[i] * 1000000 + 1), 2)
	return df
	
# 数据矩阵获得，不分析Other
classList = ["Benign", "Pathogenic"]
train = normalizationDF("train.txt", classList)
test = normalizationDF("test.txt", classList)
```



### 模型建立

使用SKlearn的LinearSVC模型。

```python
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.pipeline import make_pipeline

def LinearSVCTrain(xtrain, ytrain):
	model = make_pipeline(StandardScaler(), LinearSVC(random_state=0, tol=1e-4, max_iter=1000, C=1.0))
	clf = model.fit(xtrain, ytrain)
	return clf
	
# 获得结果
xtrain = train.iloc[:, 1:]
ytrain = train["Class"]
clf = LinearSVCTrain(xtrain, ytrain)
```



### 交叉验证

有时候需要交叉验证。
```python
from sklearn.model_selection import cross_val_score

def crossValScorePipe(x, y, cv=10):
	model = make_pipeline(StandardScaler(), LinearSVC(random_state=0, tol=1e-4, max_iter=1000, C=1.0))
	scores = cross_val_score(model, x, y, cv)
	return scores

# 打印结果
print(crossValScorePipe(xtrain, ytrain))
```



### 模型保存及读取

将模型保存下来，方便下次对其他测试数据使用。

```python
import joblib

def model_save(clf, model):
	joblib.dump(clf, model)
	
def model_load(model):
	clf = joblib.load(model)
	return clf
	
# 保存
model_save(clf, "LinearSVC.model")
# 读取
clf = model_load("LinearSVC.model")
```



### 贡献度统计

统计模型中，每个feature的贡献度，并画出前20个feature的柱状图。
```python
import matplotlib.pyplot as plt
import seaborn as sns

def feature_importance(x, model, clsName="linearsvc"):
	feature_names = x.columns.values.tolist()
	coefs = model.named_steps[clsName].coef_.flatten()
	zipped = zip(feature_names, coefs)
	df = pd.DataFrame(zipped, columns=["Feature", "Value"])
	df["Abs_value"] = df["Value"].apply(lambda x: abs(x))
	df["Colors"] = df["Value"].apply(lambda x: "green" if x > 0 else "red")
	df = df.sort_values("Abs_value", ascending=False)
	return df

def feature_plot(df, picName):
	fig, ax = plt.subplots(1, 1, figsize=(15, 12))
	sns.barplot(x="Feature", y="Value", data=df.head(20), palette=df.head(20)["Colors"])
	ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=20)
	ax.set_title("Top 20 Features", fontsize=25)
	ax.set_ylabel("Coef", fontsize=22)
	ax.set_xlabel("Feature Name", fontsize=20)
	plt.savefig(picName)

#  画图及保存
f_imp = feature_importance(xtrain, clf)
f_imp.to_csv("LinearSVC.FeaturesRanking.txt", sep="\t", index=False)
feature_plot(f_imp, "LinearSVC.Top20Features.png")
```



### 模型测试

使用模型对测试数据进行预测。
```python
def modelPredict(xtest, clf, clsList):
	predict_clf = list(clf.predict(xtest))
	predict_transform = []
	for pc in predict_clf:
		predict_transform.append(clsList[pc])
	name_test = xtest.index
	sampleName = list(name_test)
	zipped = list(zip(sampleName, predict_transform))
	return zipped

# 测试
xtest = test.iloc[:, 1:]
ytest = test["Class"]
zTest = modelPredict(xtest, clf, classList)
results = open("LinearSVC.results.txt", "w", encoding="utf-8")
for z in zTest:
	results.write(z[0] + "\t" + z[1] + "\n")
results.close()
```
