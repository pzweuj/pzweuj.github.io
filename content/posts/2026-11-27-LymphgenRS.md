---
title: LymphgenRS，单文件跑完 NCI LymphGen 2.0 淋巴瘤分型
tags: software
---

去年写过一篇 [使用 LymphPlexPy 对 DLBCL 进行分型](/posts/LymphPlexPy)，当时是为了省掉 R 环境和复杂的输入整理，把 36 基因的 LymphPlex 重构成了一个吞 JSON 的 Python 工具。不过真要跑经典的 **NCI LymphGen 2.0**（Wright et al., *Cancer Cell* 2020），官方那套 R 脚本带着将近 200 MB 的参考矩阵用起来还是偏重；而且在临床靶向 panel 上，每次为了跑它去凑齐 `Sample.annot.txt`、`Mutation flat.txt`、`Copy flat.txt` 那 6 张表实在繁琐。

干脆用 Rust 把完整的 LymphGen 2.0 重写了一遍，开源在 [LymphgenRS](https://github.com/pzweuj/LymphgenRS)。核心想法就两个：把官方参考队列直接压进单个二进制文件里实现零依赖运行，同时提供和 LymphPlexPy 一样顺手的 JSON 接口，支持 **BN2**、**EZB**（含 `EZB-MYC+` / `EZB-MYC-` 分层）、**MCD**、**N1**、**ST2**、**A53**、各类复合亚型以及 `Other`。

---

## 单文件嵌入与位图加速

LymphGen 2.0 和固定公式的分类器不太一样，它是“看菜吃饭”的动态建模——拿到一个新样本后，要根据它实际检测了哪些基因、有没有做 `BCL2` / `BCL6` FISH 或拷贝数检测，在 574 例 NCI 官方参考队列（163,427 个特征）上现场重训练层级朴素贝叶斯模型。

为了既跑得快又不用到处挂数据目录，LymphgenRS 做了两处工程优化：

- **编译期嵌入参考矩阵**：用仓库里的 `tools/pack_reference.py` 把原本约 200 MB 的 NCI 训练数据无损压缩成 3.4 MB 左右的 `data/lymphgen_ref.bin.gz`，编译时通过 `include_bytes!` 直接打进可执行文件里，运行时单二进制就能直接跑。
- **576-bit 位图重拟合**：把 574 例样本矩阵编码成 `[u64; 9]`（576-bit）的位图，利用硬件 `popcnt` 指令加速 2×2 Fisher 精确检验和条件概率统计，单次动态重训练与预测耗时压到了几十到一百多毫秒。

算法版本默认 `--algo-version 2.0`（也兼容 `1.0`），在 NCI 官方示例数据集上做过逐字符比对，输出和官方 `Result.txt` 完全一致。

---

## 面向 Panel 的三态 JSON 输入

以前写 LymphPlexPy 时输入只分 `0` 和 `1`。但在真实靶向 panel 里，“测了是阴性”和“panel 里压根没包含这个基因”完全是两回事，把未检测基因强行填成 `0` 会拉偏贝叶斯后验概率。所以 LymphgenRS 的 JSON 输入直接改成了三态编码：

- `1`（或 `true`）：阳性改变，纳入已检测集合
- `0`（或 `false`）：野生型 / 阴性，纳入已检测集合
- `null` / `"none"` / `"NA"`：未检测或无信息，直接排除出检测集合（例如缺失 `BCL2` 或拷贝数时，模型会自动切到 `NoBCL2CGH` 等对应子模型）

几类常见特征的 key 约定也很直白：

- **普通体细胞突变**：用 Gene Symbol 或 Entrez ID，如 `"TP53": 1`、`"CD79B": 1`、`"NOTCH1": 0`。如果手上是一份只列出了阳性位点的 WES 稀疏 JSON，可以在命令行加上 `--panel full`，把未列出的参考基因统一视作野生型 `0`。
- **MYD88 双位点**：因为 L265P 和非 L265P 在模型里权重不同，分别用 `"MYD88_L265P"` 和 `"MYD88_Other"` 表示（不带下划线的写法也认）。
- **FISH 易位与 `EZB-MYC` 分层**：`"BCL2"` 和 `"BCL6"` 直接参与 LymphGen 主模型计算；当预测结果包含 `EZB` 时，程序再根据 `"MYC"` 状态做二次分层——`1` 标成 `EZB-MYC+`，`0` 标成 `EZB-MYC-`，`null` 则保持 `EZB`。在 DNA + FISH 流程下，这一步只做后处理标注，不改变前面 6 个核心亚型的贝叶斯概率。
- **局灶与染色体臂 CNV**：基因层面支持 `_AMP`、`_GAIN`、`_HETLOSS`（或 `_LOSS` / `_DEL`）、`_HOMDEL` 后缀（比如 `"TP53_HETLOSS": 1`、`"CDKN2A_HOMDEL": 1`）；臂级改变支持 `1p..22q`、`Xp..Yq`、`1chrom..24chrom`（也认 `1q+`、`17p-` 这类简写）。

单样本 JSON 示例（`examples/single_sample.json`）：

```json
{
  "MYD88_L265P": 1,
  "MYD88_Other": 0,
  "CD79B": 1,
  "PIM1": 1,
  "PRDM1": 1,
  "NOTCH1": 0,
  "EZH2": 0,
  "BCL2": null,
  "BCL6": 0,
  "MYC": null
}
```

多样本直接包一层 `{ "Sample1": { ... }, "Sample2": { ... } }` 就能批量跑。

---

## 怎么跑

可以从 [GitHub Releases](https://github.com/pzweuj/LymphgenRS/releases) 直接下载 Linux（`x86_64`）、Windows（`x86_64`）或 macOS（`arm64` / `x86_64`）的预编译二进制，也可以本地 `cargo build --release` 编译。

日常几种跑法：

```bash
# 默认 LymphGen 2.0，输出 JSON
lymphgen-rs predict examples/single_sample.json

# 多样本预测并输出 TSV，跳过 574 例参考集对照评估以获取最快速度
lymphgen-rs predict examples/multi_sample.json --format tsv --no-compare

# 只列出阳性突变的稀疏 JSON，按全外显子背景将未列出基因视作野生型 0
lymphgen-rs predict positive_only.json --panel full

# 从标准输入读取 JSON
cat examples/single_sample.json | lymphgen-rs predict -

# 兼容 NCI 官方 6 文件 flat 格式
lymphgen-rs flat \
  --sample-annot Sample.annot.txt \
  --mut-flat "Mutation flat.txt" \
  --mut-genes mut.genelist.txt \
  --copy-flat "Copy flat.txt" \
  --cgh-genes CGH.genelist.txt \
  --arm-file Arm.file.txt
```

不确定某个特殊 key 或染色体臂怎么拼写时，直接跑 `lymphgen-rs list-genes` 就能把支持的清单全打出来。

如果需要直接嵌进别的 Rust 服务里，也可以作为库调用：

```rust
use lymphgen_rs::{predict_from_json_str, JsonPredictOptions, ReferenceData};

fn main() -> anyhow::Result<()> {
    let rf = ReferenceData::global();
    let json = r#"{"MYD88_L265P": 1, "MYD88_Other": 0, "CD79B": 1, "PIM1": 1, "BCL2": null, "BCL6": 0}"#;
    let result = predict_from_json_str(rf, json, &JsonPredictOptions::default())?;
    println!("Prediction: {}", result.predictions[0].subtype_prediction);
    Ok(())
}
```

---

## 许可

代码采用 **BSD 3-Clause License**（三条款 BSD 许可）开源，算法本身来自 Wright 等人 2020 年发表在 *Cancer Cell* 的 LymphGen 2.0 研究。
