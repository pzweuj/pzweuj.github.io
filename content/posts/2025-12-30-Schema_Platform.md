---
title: 2026新坑开挖-开源NGS临床分析平台
tags: default
---



项目地址：[schema-platform](https://github.com/SchemaBio/schema-platform)

项目包含“遗传全外显子分析平台”以及“实体瘤泛癌种分析平台”（暂不考虑血液瘤，比较复杂）。

以[Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0)许可开源。

**这个项目我预估会使用2026年全年的业余时间来打磨。**



## 前端

项目前端整体采用NextJs进行构建，数据的对接使用[Parquet](https://parquet.apache.org/)格式

整体采用类Github的风格

![sp_1](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/schema_platform_1.png)

UI上，会找点人脉继续聊聊，从使用方的感受来进一步打磨

![sp_2](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/schema_platform_2.png)



## 后端

整体采用Go语言建立后端，提高效率，降低资源损耗。



## 分析流程

采用[Nextflow](https://www.nextflow.io/)建立。所有流程和软件均采用docker进行封装，用docker-in-docker的方式运行。开源版本的数据库使用公共数据库，同时我自建一个订阅数据库（我的盈利点）。

当前的思路是，用户可以在这个前端自行导入，通过配置VEP参数达到加入自己注释的效果

![sp_3](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/schema_platform_3.png)



### 定制化报告

每个分析任务最终都会产生一个Parquet文件（包含所有的样本信息、检出报表、文件路径等内容），用户可以自行以Restful API构建自己的报告生成端点，对接到系统中。

![sp_4](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master/content/data/images/schema_platform_4.png)


