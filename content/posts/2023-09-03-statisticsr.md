---
title: python几种统计方法
tags: coding
---

记录一下几种统计学方法的python写法



## fisher精确检验

odds ratio：两个事件的相关度，大于1为正相关，小于1为负相关，等于1为不相关；

p-value：置信度

置信区间：odds ratio在此区间内置信度为p-value

```python
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats import norm
import math

def fisher_own(dataframe, groupAname, groupBname, yates=True):
    groupA_muts = dataframe.loc["突变", groupAname]
    groupA_unmuts = dataframe.loc["未突变", groupAname]
    groupB_muts = dataframe.loc["突变", groupBname]
    groupB_unmuts = dataframe.loc["未突变", groupBname]

    # 构建二联表
    observed = [[groupA_muts, groupB_muts], [groupA_unmuts, groupB_unmuts]]

    # 如果列联表中包含0值，需要Yates修正，将所有值增加一个平均数
    if yates:
        while 0 in [groupA_muts, groupA_unmuts, groupB_muts, groupB_unmuts]:
            ave = float(sum([groupA_muts, groupA_unmuts, groupB_muts, groupB_unmuts]) / 4)
            groupA_muts = groupA_muts + ave
            groupA_unmuts = groupA_unmuts + ave
            groupB_muts = groupB_muts + ave
            groupB_unmuts = groupB_unmuts + ave
            observed = [[groupA_muts, groupB_muts], [groupA_unmuts, groupB_unmuts]]
    odds_ratio, p_value = fisher_exact(observed)
    
    # 计算0.95置信水平下的置信区间，存在0值时，理论上fisher精确检验无法计算置信区间
    # 改用wilson方法进行计算
    alpha = 0.05
    if 0 in [groupA_muts, groupA_unmuts, groupB_muts, groupB_unmuts]:
        n = groupA_muts + groupA_unmuts + groupB_muts + groupB_unmuts
        p_hat = (groupA_muts + 0.5) / (groupA_muts + groupB_muts + 1)
        z_alpha_half = norm.ppf(1 - alpha / 2)
        denominator = 1 + (z_alpha_half ** 2) / n
        centre_adjusted_probability = p_hat + (z_alpha_half ** 2) / (2 * n)
        adjusted_standard_error = math.sqrt((p_hat * (1 - p_hat) + (z_alpha_half ** 2) / (4 * n)) / n)
        lower_bound = centre_adjusted_probability - adjusted_standard_error * z_alpha_half / denominator
        upper_bound = centre_adjusted_probability + adjusted_standard_error * z_alpha_half / denominator

    else:
        try:
            log_odds_ratio = math.log(odds_ratio)
            se_log = math.sqrt((1 / groupA_muts) + (1 / groupA_unmuts) + (1 / groupB_muts) + (1 / groupB_unmuts))
            z_alpha_half = abs(math.erf(alpha / math.sqrt(2)))
            lower_bound = math.exp(log_odds_ratio - z_alpha_half * se_log)
            upper_bound = math.exp(log_odds_ratio + z_alpha_half * se_log)
        except:
            lower_bound = upper_bound = ""

    return odds_ratio, p_value, lower_bound, upper_bound

data = {"组A": {"突变": 34, "未突变": 0}, "组B": {"突变": 10, "未突变": 20}}
df = pd.DataFrame(data)
output = fisher_own(df, "组A", "组B")
```



## Bonferroni校正

Bonferroni校正是一种用于多重比较中控制错误率的方法，当我们进行多个统计假设检验时，以防止因执行大量假设检验而增加，在对多个位点逐个进行了fisher精确检验等统计后，将p值与新的显著性水平进行比较。

```python
def bonferroni_fix(p_value_dict, alpha=0.05):
    bonfierroni_pass = {}
    n = len(p_value_dict)
    new_alpha = alpha / n
    for i in p_value_dict:
        if p_value_dict[i] < new_alpha:
            bonfierroni_pass[i] = True
        else:
            bonfierroni_pass[i] = False
    return bonfierroni_pass
```



## 卡方检验

期望值 = （行总和 x 列总和）/ 总样本数

卡方值：计算观察与期望值的差异，越大越显著，0值时不显著

自由度：可用独立变化的数据个数，一般是行数减去1再乘以列数减去1

P值：观察的数据与原假设的矛盾程度，如果小于显著性水平（一般是0.05）则拒绝原假设，并认为观察到的结果是显著的

贡献度：某些情况下不同因素对结果的贡献程度，越大则贡献度越大，负值表示减少了整体差异
```python
import pandas as pd
import numpy as np
from scipy.stats import chi2

def chi2_own(dataframe):
    observed = dataframe.T.values
    
    # 计算期望频次
    row_totals = observed.sum(axis=1)
    col_totals = observed.sum(axis=0)
    total = observed.sum()
    expected = np.outer(row_totals, col_totals) / total

    # 计算卡方值
    chi_squared = np.sum((observed - expected) ** 2 / expected)

    # 计算自由度和显著性水平
    degrees_of_freedom = (observed.shape[0] - 1) * (observed.shape[1] - 1)
    p_value = 1 - chi2.cdf(chi_squared, degrees_of_freedom)

    # 计算贡献度
    with np.errstate(divide="ignore", invalid="ignore"):
        contribution = (observed / total) * (np.log(observed / expected))
        # nan值改为0
        contribution[np.isnan(contribution)] = 0
    
    df_contri = pd.DataFrame(contribution).T
    df_contri.index = dataframe.index
    df_contri.columns = dataframe.columns
    
    return chi_squared, p_value, degrees_of_freedom, df_contri


data = {"组A": {"位点1": 34, "位点2": 0, "位点3": 29, "位点4": 10}, "组B": {"位点1": 4, "位点2": 10, "位点3": 20, "位点4": 10}}
df = pd.DataFrame(data)
output = chi2_own(df)
```



## 蒙特卡洛模拟
个别位点仅在单个样本中突变，计算所得卡方值比较大，因此使用Monte Carlo模拟数据，并计算一个卡方值的阈值，认为95分位数时，显著，实际操作中卡方值以此为上限。这里使用了上面的chi2_own函数。

```python
import random
import numpy as np

def monte_carlo_sim(act_data_list):
    # observed_data = np.array(act_data_list)
    num_simulations = 10000
    chi2_values = []
    for n in range(num_simulations):
        # 模拟所得的虚拟样本
        simulated_data = random.choice(act_data_list)
        chi_squared, p_value, _, _ = chi2_own(simulated_data)
        if not str(chi_squared) == "nan":
            chi2_values.append(chi_squared)

    # 模拟过程中观察到的卡方值分布中的 95 分位数
    percentile_95 = np.percentile(chi2_values, 95)
    print("观察卡方值分布中的 95 分位数:", percentile_95)
```



## ANOVA多因素方差分析
当进行多组比较时，显然不适用单位点的假设检验，因此这里需要直接构建一个多位点的矩阵

```python
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from scipy.stats import kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd

def anova_own(dataframe):
    df = dataframe.T
    anova_results = f_oneway(*[df[group] for group in df.columns])
    anova_statistic, anova_pvalue = anova_results

    # 如果方差分析的p值低于显著性水平，进行多重比较校正
    alpha = 0.05
    if anova_pvalue < alpha:
        tukey_results = pairwise_tukeyhsd(df.values.flatten(), np.repeat(df.columns, df.shape[0]))

    # 计算每个位点对于分组差异的贡献度
    mean_per_group = df.mean()
    total_mean = df.values.mean()
    contribution = (mean_per_group - total_mean).abs() / total_mean

    # 输出位点的p值
    p_values = {}
    for column in df.columns:
        group_data = [df[column][group] for group in df.T.columns]
        _, p_value = kruskal(*group_data)
        p_values[column] = p_value
    return anova_results, contribution, p_values


data = {
    "组A": {"位点1": 34, "位点2": 10, "位点3": 0, "位点4": 24},
    "组B": {"位点1": 0, "位点2": 20, "位点3": 13, "位点4": 77},
    "组C": {"位点1": 6, "位点2": 15, "位点3": 1, "位点4": 5},
    "组D": {"位点1": 1, "位点2": 0, "位点3": 2, "位点4": 1}
}
df = pd.DataFrame(data)
```



## Shapiro-Wilk test

正态分布检验 p≥0.05正态分布，p＜0.05非正态分布，对每个组的数据进行Shapiro-Wilk检验。

```python
from scipy.stats import shapiro

def shapiro_wilk_own(dataframe):
    groupDict = {}
    for group in dataframe.columns:
        groupDict[group] = shapiro(dataframe[group].values)
    return groupDict
```



## Levene's test

方差齐性检验 p≥0.05方差齐，p＜0.05方差齐。

```python
from scipy.stats import levene

def levene_own(dataframe):
    dataGroupList = []
    for group in dataframe.columns:
        dataGroupList.append(dataframe[group].values)
    lev_results = levene(*dataGroupList)
    return lev_results
```







