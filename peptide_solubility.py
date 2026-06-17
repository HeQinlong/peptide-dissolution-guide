#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
多肽溶解性分析工具
根据 Sigma-Aldrich 指南评估多肽电荷性质并给出溶解方法建议
参考: https://www.sigmaaldrich.cn/CN/zh/technical-documents/protocol/protein-biology/
      protein-and-nucleic-acid-interactions/peptide-solubility
"""

# ──────────────────────────────────────────────
# 氨基酸分类数据
# ──────────────────────────────────────────────

# 疏水性残基（用于判断是否为疏水肽）（非精确热力学分类）
# 核心来源：
#   - Kyte & Doolittle (1982) J.Mol.Biol. 157:105-132（A/C/F/I/L/M/V 正疏水值）
#   - Pierce/Thermo Fisher 肽溶解指南（新增 P/W/Y 作为溶解性疏水残基）
#   - Merck Millipore 氨基酸参考表（F/W/Y 芳香族疏水）
# 注：C 在 Kyte-Doolittle 中疏水，但含游离巯基，DMSO 中不稳定（另有警告）
#     P 在 K-D 中略偏亲水(-1.6)，但溶解性实践中普遍归为疏水
#     W/Y 在 K-D 中偏中性，但因苯环结构在溶解指南中归为疏水
HYDROPHOBIC_AA = set("ACFILMPVWY")

# 带电残基：酸性（-1）和碱性（+1）
ACIDIC_AA  = set("DE")   # Asp, Glu → 净电荷 -1
BASIC_AA   = set("KR")   # Lys, Arg → 净电荷 +1
# His 在 pH<6 时带电 +1，pH≥6 时为 0；本程序按 pH7 计算，故 H 不计入净电荷

# 全部有效单字母氨基酸代码
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def validate_sequence(seq: str) -> str:
    """
    验证并清理输入序列。
    返回大写、去空格的合法序列；遇到非法字符则抛出 ValueError。
    """
    seq = seq.strip().upper().replace(" ", "")
    invalid = set(seq) - VALID_AA
    if invalid:
        raise ValueError(f"序列包含非法字符: {', '.join(sorted(invalid))}")
    if len(seq) == 0:
        raise ValueError("序列不能为空。")
    return seq


def calculate_charge(seq: str) -> dict:
    """
    计算多肽在 pH 7 时的净电荷及相关统计。

    计分规则（参考 Sigma-Aldrich 指南）：
      - 每个酸性残基 D、E → -1
      - 每个碱性残基 K、R → +1
      - N 末端 NH2        → +1
      - C 末端 COOH       → -1
      - H（His）在 pH 7   →  0
      注：Cys（C）游离巯基 pKa ≈ 8，pH7 时通常不去质子化，不计入净电荷。

    返回包含计数和净电荷的字典。
    """
    n_acidic  = sum(1 for aa in seq if aa in ACIDIC_AA)
    n_basic   = sum(1 for aa in seq if aa in BASIC_AA)
    n_hydro   = sum(1 for aa in seq if aa in HYDROPHOBIC_AA)

    # 末端贡献：N 末端 +1，C 末端 -1
    terminal_charge = +1 + (-1)  # = 0，两者相抵

    net_charge = n_basic - n_acidic + terminal_charge  # 末端净贡献为 0

    # 总带电残基数（含末端）用于比例计算
    # 带电残基 = D + E + K + R（pH7 下 H 不带电）
    total_charged = n_acidic + n_basic
    # 加上末端（N 末端 +1、C 末端 -1 各算一个带电基团）
    total_charged_with_termini = total_charged + 2

    total_residues = len(seq)
    charge_ratio = total_charged_with_termini / total_residues  # 带电比例

    hydro_ratio = n_hydro / total_residues  # 疏水残基比例

    return {
        "sequence":      seq,
        "length":        total_residues,
        "n_acidic":      n_acidic,
        "n_basic":       n_basic,
        "n_hydrophobic": n_hydro,
        "net_charge":    net_charge,
        "total_charged": total_charged_with_termini,
        "charge_ratio":  charge_ratio,
        "hydro_ratio":   hydro_ratio,
    }


def classify_peptide(info: dict) -> dict:
    """
    根据电荷信息将多肽分类，并给出溶解建议。

    判断逻辑（依照 Sigma-Aldrich 指南）：
    1. 若疏水残基 > 50%，或带电比例 < 10%
       → 疏水/不带电肽，建议有机溶剂（ACN / DMSO / DMF）。
    2. 若带电比例 ≥ 10%，进入带电肽判断：
       a. 净电荷 < 0  → 酸性肽  → 0.1 M 碳酸氢铵溶解
       b. 净电荷 > 0  → 碱性肽  → 25% 乙酸溶解
       c. 净电荷 = 0  → 中性肽
          - 带电比例 > 25% → 按酸性肽方法（碳酸氢铵）
          - 带电比例 10%–25% → 按碱性肽方法（乙酸）
    """
    cr   = info["charge_ratio"]
    hr   = info["hydro_ratio"]
    nc   = info["net_charge"]
    seq  = info["sequence"]

    # ── 疏水/不带电肽 ──
    has_cys = "C" in seq
    has_met = "M" in seq

    if hr > 0.50 or cr < 0.10:
        peptide_type  = "疏水/不带电肽"
        charge_class  = "—"
        method_lines  = [
            "建议使用有机溶剂：",
            "  1. 乙腈（ACN）",
            "  2. 二甲基甲酰胺（DMF）",
            "  3. 二甲基亚砜（DMSO）",
            "     ⚠ 含 Cys（C）或 Met（M）的肽在 DMSO 中不稳定，请避免使用。" if (has_cys or has_met) else "",
            "  4. 若仍难溶，加离液序列高的化合物（盐酸胍 / 尿素）破坏氢键网络。",
            "注意：先将肽完全溶于有机溶剂，再逐滴缓慢加入缓冲液并轻轻搅拌。",
        ]
    else:
        # ── 带电肽 ──
        peptide_type = "带电肽"

        if nc < 0 or cr>0.25:
            charge_class = "酸性肽（净电荷为负），且/或该肽在pH 7下的总电荷数大于残基总数的 25％"
            method_lines = [
                "加入少量 0.1 M 碳酸氢铵溶解肽，再用水稀释至所需浓度。",
                "确保肽溶液最终 pH ≈ 7，必要时调节 pH。",
            ]
        elif nc > 0 or 0.1<=cr<=0.25:
            charge_class = "碱性肽（净电荷为正），且/或该肽在pH 7下的总电荷数介于残基总数的 10% 至 25％"
            method_lines = [
                "加入少量 25% 乙酸溶解肽，再用水稀释至所需浓度。",
            ]
        else:  # nc == 0，中性肽
            charge_class = "中性肽（净电荷为零）"
            if cr > 0.25:
                method_lines = [
                    "带电残基比例 > 25%，按酸性肽方法处理：",
                    "加入少量 0.1 M 碳酸氢铵溶解肽，再用水稀释至所需浓度。",
                    "确保最终 pH ≈ 7，必要时调节 pH。",
                ]
            else:
                method_lines = [
                    "建议使用有机溶剂：",
                    "  1. 乙腈（ACN）",
                    "  2. 二甲基甲酰胺（DMF）",
                    "  3. 二甲基亚砜（DMSO）",
                    "     ⚠ 含 Cys（C）或 Met（M）的肽在 DMSO 中不稳定，请避免使用。" if (has_cys or has_met) else "",
                    "  4. 若仍难溶，加离液序列高的化合物（盐酸胍 / 尿素）破坏氢键网络。",
                    "注意：先将肽完全溶于有机溶剂，再逐滴缓慢加入缓冲液并轻轻搅拌。",
                ]

    # 过滤空行
    method_lines = [l for l in method_lines if l]

    return {
        "peptide_type": peptide_type,
        "charge_class": charge_class,
        "method":       method_lines,
    }


def print_report(info: dict, result: dict) -> None:
    """格式化打印分析报告。"""
    sep = "=" * 55

    print(f"\n{sep}")
    print("          多肽溶解性分析报告")
    print(sep)

    print(f"  序列（长度 {info['length']}）: {info['sequence']}")
    print(f"  酸性残基（D/E）数量  : {info['n_acidic']}")
    print(f"  碱性残基（K/R）数量  : {info['n_basic']}")
    print(f"  疏水残基数量         : {info['n_hydrophobic']}")
    print(f"  pH 7 净电荷          : {info['net_charge']:+d}")
    print(f"  带电残基比例（含末端）: {info['charge_ratio']:.1%}")
    print(f"  疏水残基比例         : {info['hydro_ratio']:.1%}")
    print()
    print(f"  ▶ 肽的类型  : {result['peptide_type']}")
    if result['charge_class'] != "—":
        print(f"  ▶ 电荷类型  : {result['charge_class']}")
    print()
    print("  ▶ 溶解方法建议：")
    for line in result["method"]:
        print(f"    {line}")
    print(sep)


def analyze_peptide(sequence: str) -> None:
    """主流程：验证序列 → 计算电荷 → 分类 → 输出报告。"""
    try:
        seq    = validate_sequence(sequence)
        info   = calculate_charge(seq)
        result = classify_peptide(info)
        print_report(info, result)
    except ValueError as e:
        print(f"\n[输入错误] {e}")


# ──────────────────────────────────────────────
# 主程序入口
# ──────────────────────────────────────────────

def main():
    print("=" * 55)
    print("  多肽溶解特性和溶解方法分析工具（根据 Sigma-Aldrich 合成肽的处理和储存实验方案）")
    print("  参考文献：")
    print("  [1] 合成肽的处理和储存实验方案[EB/OL]. [2026-06-16]. https://www.sigmaaldrich.cn/CN/zh/technical-documents/protocol/protein-biology/protein-and-nucleic-acid-interactions/peptide-solubility.\
")
    print("  [2] Pommié C，Levadoux S，Sabatier R，et al. IMGT standardized criteria for statistical analysis of immunoglobulin V-REGION amino acid properties[J]. Journal of molecular recognition: JMR，2004，17（1）：17-32. \
")
    print("  重要提示：建议先尝试水和乙酸，若不溶，可以冻干后尝试别的试剂")
    print("  输入单字母氨基酸序列（如 ACDEFGHIK）")
    print("  输入 'Q'、'q' 或 'quit' 退出")
    print("=" * 55)

    while True:
        raw = input("\n请输入多肽序列: ").strip()
        if raw.lower() in ("q","Q", "quit", "exit", ""):
            print("已退出。")
            break
        analyze_peptide(raw)


if __name__ == "__main__":
    main()
