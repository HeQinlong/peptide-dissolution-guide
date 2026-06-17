# 🧪 多肽溶解指导

一个轻量级的命令行工具，根据多肽序列在 pH 7 下的电荷性质，评估其溶解特性并给出具体的溶解方法建议。分类逻辑完全遵循 **Sigma-Aldrich 合成肽处理与储存实验方案**。

---

## 功能特性

- ✅ 验证标准单字母氨基酸序列，自动过滤非法字符
- ⚡ 计算 pH 7 下的净电荷（酸性残基 D/E、碱性残基 K/R，含 N/C 末端贡献）
- 📊 输出疏水残基比例与带电残基比例
- 🔬 将多肽分类为：疏水/不带电肽、酸性肽、碱性肽或中性肽
- 💊 给出具体溶解溶剂建议（水、碳酸氢铵、乙酸、ACN / DMF / DMSO）
- ⚠️ 对含 Cys（C）或 Met（M）的肽给出 DMSO 不相容警告

---

## 使用方法

**运行环境要求：** Python 3.6+，无需安装任何第三方依赖。

```bash
python peptide_solubility.py
```

在提示符下输入单字母氨基酸序列（如 `ACDEFGHIK`），输入 `q` 或 `quit` 退出。

### 示例输出

```
=======================================================
 多肽溶解性分析报告
=======================================================
 序列（长度 9）: ACDEFGHIK
 酸性残基（D/E）数量   : 2
 碱性残基（K/R）数量   : 1
 疏水残基数量          : 3
 pH 7 净电荷           : -1
 带电残基比例（含末端）: 55.6%
 疏水残基比例          : 33.3%

 ▶ 肽的类型 : 带电肽
 ▶ 电荷类型 : 酸性肽（净电荷为负）

 ▶ 溶解方法建议：
   加入少量 0.1 M 碳酸氢铵溶解肽，再用水稀释至所需浓度。
   确保肽溶液最终 pH ≈ 7，必要时调节 pH。
=======================================================
```

---

## 分类判断逻辑

| 判断条件 | 肽的类型 | 推荐溶剂 |
|---|---|---|
| 疏水残基比例 > 50%，或带电比例 < 10% | 疏水/不带电肽 | ACN / DMF / DMSO |
| 净电荷 < 0 | 酸性肽 | 0.1 M NH₄HCO₃（碳酸氢铵） |
| 净电荷 > 0 | 碱性肽 | 25% 乙酸 |
| 净电荷 = 0，带电比例 > 25% | 中性极性肽 | 0.1 M NH₄HCO₃（碳酸氢铵） |
| 净电荷 = 0，带电比例 10%–25% | 中性肽 | 25% 乙酸 |

> 💡 重要提示：建议先尝试用水溶解；若难溶，可冻干后再尝试其他溶剂。

---

## 参考文献

1. Sigma-Aldrich. *合成肽的处理和储存实验方案* \[EB/OL\]. \[2026-06-16\]. https://www.sigmaaldrich.cn/CN/zh/technical-documents/protocol/protein-biology/protein-and-nucleic-acid-interactions/peptide-solubility
2. Pommié C, Levadoux S, Sabatier R, et al. IMGT standardized criteria for statistical analysis of immunoglobulin V-REGION amino acid properties\[J\]. *Journal of Molecular Recognition*, 2004, 17(1): 17–32.

---

## 许可证

MIT License
