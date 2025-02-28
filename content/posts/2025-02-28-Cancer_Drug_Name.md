---
title: 癌症靶向药物名称
tags: default
---

## 整理流程

这个库整理自**[新型抗肿瘤药物临床应用指导原则(2024)](http://www.nhc.gov.cn/yzygj/s7659/202501/c1b6a86e2e6040aca75fdf89bc382184/files/e8eb67d6cd0444bbbde66da45d956d3c.pdf)**，使用通义千问 Qwen 2.5VL 72B识别对中文靶向药名进行提取，提取后进行了人工校对（这个原则中的中文名称也是有问题的，主要在音译后中文的不统一，没有按照获批上市的名称作为正确名称），可能存在纰漏。

接下来我使用了Gemini 2.0 Flash和通义千问 Qwen 2.5 Max将中文翻译为英文，对于两个模型翻译一致的结果，直接PASS，对不一致的结果进行了人工校对（Google）。【注，没有联网搜索】

## 核对过程

有一些奇怪的结果，譬如**瑞普替尼**和**瑞派替尼**，Gemini和通义千问都认为英文名称都是Ripretinib，我原本以为是中文音译问题，但是查了一下，这俩又好像是两种药。最后的结论是，用于**肺癌**的**瑞普替尼**的英文名是**Repotrectinib**，而用于**消化道肿瘤**的**瑞派替尼**的英文名是**Ripretinib**。

## 赛后结算

在全部147种靶向药中，双方翻译一致的有110种，不一致的是37种。

在这37种药物种，与我最终校正的结果对比，通义千问答对了13个，正确率35.14%；Gemini答对了6个，正确率16.22%。当然，如果放到所有147种药物看，通义千问正确率是83.67%，Gemini是78.91%。


## 整理结果

整理所得的药物表格在下面，同时提供一个[Excel表格](https://github.com/pzweuj/pzweuj.github.io/raw/refs/heads/master/content/data/document/cancer_drug_name_2024.xlsx)。


| 中文名称 | 癌种 | 英文名称 |
| :-: | :-: | :-: |
| 恩曲替尼 | 泛实体瘤 | Entrectinib |
| 思沃利单抗 | 泛实体瘤 | Envafolimab |
| 拉罗替尼 | 泛实体瘤 | Larotrectinib |
| 帕博利珠单抗 | 泛实体瘤 | Pembrolizumab |
| 普特利单抗 | 泛实体瘤 | Pucotenlimab |
| 斯鲁利单抗 | 泛实体瘤 | Serplulimab |
| 替雷利珠单抗 | 泛实体瘤 | Tislelizumab |
| 安罗替尼 | 骨与软组织肿瘤 | Anlotinib |
| 地舒单抗 | 骨与软组织肿瘤 | Denosumab |
| 依维莫司 | 骨与软组织肿瘤 | Everolimus |
| 他泽司他 | 骨与软组织肿瘤 | Tazemetostat |
| 阿得贝利单抗 | 呼吸系统肿瘤 | Adebrelimab |
| 阿法替尼 | 呼吸系统肿瘤 | Afatinib |
| 阿来替尼 | 呼吸系统肿瘤 | Alectinib |
| 阿美替尼 | 呼吸系统肿瘤 | Almonertinib |
| 安罗替尼 | 呼吸系统肿瘤 | Anlotinib |
| 阿替利珠单抗 | 呼吸系统肿瘤 | Atezolizumab |
| 贝福替尼 | 呼吸系统肿瘤 | Befotertinib |
| 贝莫苏拜单抗 | 呼吸系统肿瘤 | Benmelstobart |
| 贝伐珠单抗 | 呼吸系统肿瘤 | Bevacizumab |
| 伯瑞替尼 | 呼吸系统肿瘤 | Bozitinib |
| 布格替尼 | 呼吸系统肿瘤 | Brigatinib |
| 卡瑞利珠单抗 | 呼吸系统肿瘤 | Camrelizumab |
| 卡马替尼 | 呼吸系统肿瘤 | Capmatinib |
| 塞瑞替尼 | 呼吸系统肿瘤 | Ceritinib |
| 克唑替尼 | 呼吸系统肿瘤 | Crizotinib |
| 达拉非尼 | 呼吸系统肿瘤 | Dabrafenib |
| 达可替尼 | 呼吸系统肿瘤 | Dacomitinib |
| 度伐利尤单抗 | 呼吸系统肿瘤 | Durvalumab |
| 恩沙替尼 | 呼吸系统肿瘤 | Ensartinib |
| 恩曲替尼 | 呼吸系统肿瘤 | Entrectinib |
| 依奉阿克 | 呼吸系统肿瘤 | Envonalkib |
| 厄洛替尼 | 呼吸系统肿瘤 | Erlotinib |
| 依维莫司 | 呼吸系统肿瘤 | Everolimus |
| 伏美替尼 | 呼吸系统肿瘤 | Furmonertinib |
| 吉非替尼 | 呼吸系统肿瘤 | Gefitinib |
| 谷美替尼 | 呼吸系统肿瘤 | Glumetinib |
| 埃克替尼 | 呼吸系统肿瘤 | Icotinib |
| 伊匹木单抗 | 呼吸系统肿瘤 | Ipilimumab |
| 伊鲁阿克 | 呼吸系统肿瘤 | Iruplinalkib |
| 依沃西单抗 | 呼吸系统肿瘤 | Ivonescimab |
| 洛拉替尼 | 呼吸系统肿瘤 | Lorlatinib |
| 纳武利尤单抗 | 呼吸系统肿瘤 | Nivolumab |
| 奥希替尼 | 呼吸系统肿瘤 | Osimertinib |
| 帕博利珠单抗 | 呼吸系统肿瘤 | Pembrolizumab |
| 派安普利单抗 | 呼吸系统肿瘤 | Penpulimab |
| 普拉替尼 | 呼吸系统肿瘤 | Pralsetinib |
| 重组人血管内皮抑制素 | 呼吸系统肿瘤 | Recombinant Human Endostatin |
| 瑞普替尼 | 呼吸系统肿瘤 | Repotrectinib |
| 瑞齐替尼 | 呼吸系统肿瘤 | Rezivertinib |
| 瑞厄替尼 | 呼吸系统肿瘤 | Rilertinib |
| 赛沃替尼 | 呼吸系统肿瘤 | Savolitinib |
| 塞普替尼 | 呼吸系统肿瘤 | Selpercatinib |
| 斯鲁利单抗 | 呼吸系统肿瘤 | Serplulimab |
| 信迪利单抗 | 呼吸系统肿瘤 | Sintilimab |
| 舒格利单抗 | 呼吸系统肿瘤 | Sugemalimab |
| 舒沃替尼 | 呼吸系统肿瘤 | Sunvozertinib |
| 特泊替尼 | 呼吸系统肿瘤 | Tepotinib |
| 替雷利珠单抗 | 呼吸系统肿瘤 | Tislelizumab |
| 特瑞普利单抗 | 呼吸系统肿瘤 | Toripalimab |
| 曲美替尼 | 呼吸系统肿瘤 | Trametinib |
| 安奈克替尼 | 呼吸系统肿瘤 | Unecritinib |
| 阿比特龙 | 泌尿系统肿瘤 | Abiraterone |
| 阿帕他胺 | 泌尿系统肿瘤 | Apalutamide |
| 阿昔替尼 | 泌尿系统肿瘤 | Axitinib |
| 达罗他胺 | 泌尿系统肿瘤 | Darolutamide |
| 维迪西妥单抗 | 泌尿系统肿瘤 | Disitamab vedotin |
| 恩扎卢胺 | 泌尿系统肿瘤 | Enzalutamide |
| 依维莫司 | 泌尿系统肿瘤 | Everolimus |
| 仑伐替尼 | 泌尿系统肿瘤 | Lenvatinib |
| 纳武利尤单抗 | 泌尿系统肿瘤 | Nivolumab |
| 奥拉帕利 | 泌尿系统肿瘤 | Olaparib |
| 培唑帕尼 | 泌尿系统肿瘤 | Pazopanib |
| 帕博利珠单抗 | 泌尿系统肿瘤 | Pembrolizumab |
| 瑞维鲁胺 | 泌尿系统肿瘤 | Rezvilutamide |
| 索拉非尼 | 泌尿系统肿瘤 | Sorafenib |
| 舒尼替尼 | 泌尿系统肿瘤 | Sunitinib |
| 替雷利珠单抗 | 泌尿系统肿瘤 | Tislelizumab |
| 特瑞普利单抗 | 泌尿系统肿瘤 | Toripalimab |
| 伏罗尼布 | 泌尿系统肿瘤 | Vorolanib |
| 达拉非尼 | 皮肤肿瘤 | Dabrafenib |
| 伊马替尼 | 皮肤肿瘤 | Imatinib |
| 帕博利珠单抗 | 皮肤肿瘤 | Pembrolizumab |
| 普特利单抗 | 皮肤肿瘤 | Pucotenlimab |
| 索立德吉 | 皮肤肿瘤 | Sonidegib |
| 特瑞普利单抗 | 皮肤肿瘤 | Toripalimab |
| 曲美替尼 | 皮肤肿瘤 | Trametinib |
| 妥拉美替尼 | 皮肤肿瘤 | Tunlametinib |
| 维莫非尼 | 皮肤肿瘤 | Vemurafenib |
| 阿贝西利 | 乳腺癌 | Abemaciclib |
| 西达本胺 | 乳腺癌 | Chidamide |
| 达尔西利 | 乳腺癌 | Dalpiciclib |
| 恩替司他 | 乳腺癌 | Entinostat |
| 依维莫司 | 乳腺癌 | Everolimus |
| 伊尼妥单抗 | 乳腺癌 | Inetetamab |
| 拉帕替尼 | 乳腺癌 | Lapatinib |
| 奈拉替尼 | 乳腺癌 | Neratinib |
| 哌柏西利 | 乳腺癌 | Palbociclib |
| 帕博利珠单抗 | 乳腺癌 | Pembrolizumab |
| 帕妥珠单抗 | 乳腺癌 | Pertuzumab |
| 帕妥珠曲妥珠单抗(皮下注射) | 乳腺癌 | Pertuzumab and Trastuzumab (Subcutaneous Injection) |
| 吡咯替尼 | 乳腺癌 | Pyrotinib |
| 瑞波西利 | 乳腺癌 | Ribociclib |
| 戈沙妥珠单抗 | 乳腺癌 | Sacituzumab Govitecan |
| 特瑞普利单抗 | 乳腺癌 | Toripalimab |
| 曲妥珠单抗 | 乳腺癌 | Trastuzumab |
| 德曲妥珠单抗 | 乳腺癌 | Trastuzumab deruxtecan |
| 恩美曲妥珠单抗 | 乳腺癌 | Trastuzumab Emtansine |
| 贝伐珠单抗 | 生殖系统肿瘤 | Bevacizumab |
| 卡度尼利单抗 | 生殖系统肿瘤 | Cadonilimab |
| 氟唑帕利 | 生殖系统肿瘤 | Fluzoparib |
| 尼拉帕利 | 生殖系统肿瘤 | Niraparib |
| 奥拉帕利 | 生殖系统肿瘤 | Olaparib |
| 帕米帕利 | 生殖系统肿瘤 | Pamiparib |
| 索卡佐利单抗 | 生殖系统肿瘤 | Socazolimab |
| 赛帕利单抗 | 生殖系统肿瘤 | Zimberelimab |
| 安罗替尼 | 头颈部肿瘤 | Anlotinib |
| 卡瑞利珠单抗 | 头颈部肿瘤 | Camrelizumab |
| 西妥昔单抗 | 头颈部肿瘤 | Cetuximab |
| 仑伐替尼 | 头颈部肿瘤 | Lenvatinib |
| 尼妥珠单抗 | 头颈部肿瘤 | Nimotuzumab |
| 纳武利尤单抗 | 头颈部肿瘤 | Nivolumab |
| 帕博利珠单抗 | 头颈部肿瘤 | Pembrolizumab |
| 普拉替尼 | 头颈部肿瘤 | Pralsetinib |
| 塞普替尼 | 头颈部肿瘤 | Selpercatinib |
| 索拉非尼 | 头颈部肿瘤 | Sorafenib |
| 昔雷利珠单抗 | 头颈部肿瘤 | Tislelizumab |
| 特瑞普利单抗 | 头颈部肿瘤 | Toripalimab |
| 阿帕替尼 | 消化系统肿瘤 | Apatinib |
| 阿替利珠单抗 | 消化系统肿瘤 | Atezolizumab |
| 阿伐替尼 | 消化系统肿瘤 | Avapritinib |
| 贝伐珠单抗 | 消化系统肿瘤 | Bevacizumab |
| 卡度尼利单抗 | 消化系统肿瘤 | Cadonilimab |
| 卡瑞利珠单抗 | 消化系统肿瘤 | Camrelizumab |
| 西妥昔单抗 | 消化系统肿瘤 | Cetuximab |
| 维迪西妥单抗 | 消化系统肿瘤 | Disitamab vedotin |
| 多纳非尼 | 消化系统肿瘤 | Donafenib |
| 度伐利尤单抗 | 消化系统肿瘤 | Durvalumab |
| 恩沃利单抗 | 消化系统肿瘤 | Envafolimab |
| 依维莫司 | 消化系统肿瘤 | Everolimus |
| 呋喹替尼 | 消化系统肿瘤 | Fruquintinib |
| 伊马替尼 | 消化系统肿瘤 | Imatinib |
| 仑伐替尼 | 消化系统肿瘤 | Lenvatinib |
| 尼妥珠单抗 | 消化系统肿瘤 | Nimotuzumab |
| 纳武利尤单抗 | 消化系统肿瘤 | Nivolumab |
| 帕博利珠单抗 | 消化系统肿瘤 | Pembrolizumab |
| 佩米替尼 | 消化系统肿瘤 | Pemigatinib |
| 普特利单抗 | 消化系统肿瘤 | Pucotenlimab |
| 雷莫西尤单抗 | 消化系统肿瘤 | Ramucirumab |
| 瑞戈非尼 | 消化系统肿瘤 | Regorafenib |
| 瑞派替尼 | 消化系统肿瘤 | Ripretinib |
| 斯鲁利单抗 | 消化系统肿瘤 | Serplulimab |
| 信迪利单抗 | 消化系统肿瘤 | Sintilimab |
| 索拉非尼 | 消化系统肿瘤 | Sorafenib |
| 舒格利单抗 | 消化系统肿瘤 | Sugemalimab |
| 舒尼替尼 | 消化系统肿瘤 | Sunitinib |
| 索凡替尼 | 消化系统肿瘤 | Surufatinib |
| 替雷利珠单抗 | 消化系统肿瘤 | Tislelizumab |
| 特瑞普利单抗 | 消化系统肿瘤 | Toripalimab |
| 曲妥珠单抗 | 消化系统肿瘤 | Trastuzumab |
| 德曲妥珠单抗 | 消化系统肿瘤 | Trastuzumab deruxtecan |
| 阿基仑赛 | 血液肿瘤 | Axicabtagene ciloleucel |
| 贝林妥欧单抗 | 血液肿瘤 | Blinatumomab |
| 硼替佐米 | 血液肿瘤 | Bortezomib |
| 维布妥昔单抗 | 血液肿瘤 | Brentuximab Vedotin |
| 卡瑞利珠单抗 | 血液肿瘤 | Camrelizumab |
| 卡非佐米 | 血液肿瘤 | Carfilzomib |
| 西达本胺 | 血液肿瘤 | Chidamide |
| 达雷妥尤单抗 | 血液肿瘤 | Daratumumab |
| 达沙替尼 | 血液肿瘤 | Dasatinib |
| 度维利塞 | 血液肿瘤 | Duvelisib |
| 氟马替尼 | 血液肿瘤 | Flumatinib |
| 吉瑞替尼 | 血液肿瘤 | Gilteritinib |
| 格菲妥单抗 | 血液肿瘤 | Glofitamab |
| 戈利昔替尼 | 血液肿瘤 | Golidocitinib |
| 伊布替尼 | 血液肿瘤 | Ibrutinib |
| 伊基奥仑赛 | 血液肿瘤 | Idecabtagene vicleucel |
| 伊马替尼 | 血液肿瘤 | Imatinib |
| 奥加伊妥珠单抗 | 血液肿瘤 | Inotuzumab ozogamicin |
| 艾伏尼布 | 血液肿瘤 | Ivosidenib |
| 伊沙佐米 | 血液肿瘤 | Ixazomib |
| 来那度胺 | 血液肿瘤 | Lenalidomide |
| 林普利塞 | 血液肿瘤 | Linperlisib |
| 莫格利珠单抗 | 血液肿瘤 | Mogamulizumab |
| 尼洛替尼 | 血液肿瘤 | Nilotinib |
| 奥妥珠单抗 | 血液肿瘤 | Obinutuzumab |
| 奥雷巴替尼 | 血液肿瘤 | Olverembatinib |
| 奥布替尼 | 血液肿瘤 | Orelabrutinib |
| 派安普利单抗 | 血液肿瘤 | Penpulimab |
| 维泊妥珠单抗 | 血液肿瘤 | Polatuzumab |
| 泊马度胺 | 血液肿瘤 | Pomalidomide |
| 瑞基奥仑赛 | 血液肿瘤 | Relmacabtagene Autoleucel |
| 瑞帕妥单抗 | 血液肿瘤 | Ripertamab |
| 利妥昔单抗 | 血液肿瘤 | Rituximab |
| 罗培干扰素α-2b | 血液肿瘤 | Ropeginterferon alfa-2b |
| 芦可替尼 | 血液肿瘤 | Ruxolitinib |
| 塞利尼索 | 血液肿瘤 | Selinexor |
| 司妥昔单抗 | 血液肿瘤 | Siltuximab |
| 信迪利单抗 | 血液肿瘤 | Sintilimab |
| 特立妥单抗 | 血液肿瘤 | Teclistamab |
| 沙利度胺 | 血液肿瘤 | Thalidomide |
| 替雷利珠单抗 | 血液肿瘤 | Tislelizumab |
| 维奈克拉 | 血液肿瘤 | Venetoclax |
| 泽布替尼 | 血液肿瘤 | Zanubrutinib |
| 泽沃基奥仑赛 | 血液肿瘤 | Zevorcabtagene Autoleucel |
| 赛帕利单抗 | 血液肿瘤 | Zimberelimab |
| 泽贝妥单抗 | 血液肿瘤 | Zuberitamab |


## 叠甲

注意，由于我没有去核对Gemini与通义千问一致的结果，因此很可能存在谬误。后续会将这个库应用于文献、指南等文件中的用药信息快速提取。

