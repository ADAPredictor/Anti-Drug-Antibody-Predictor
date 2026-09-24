"""Ensemble immunogenicity risk score (0-100).

The rule-based model combines mechanistically motivated components:
MHC-I / MHC-II epitope burden, linear B-cell epitope content, aggregation-
prone hydrophobic surface, sequence liability motifs, a modality prior
calibrated on the curated FDA dataset, and sequence novelty versus known
protein drugs.  A scikit-learn logistic-regression trainer can replace the
weights once enough labelled sequences are available (see ``train_model``).
"""

# prior immunogenicity by modality, roughly calibrated on curated FDA data
MODALITY_PRIOR = {
    "murine_mab": 0.90,
    "chimeric_mab": 0.60,
    "humanized_mab": 0.40,
    "fab_pegylated": 0.45,
    "fc_fusion": 0.45,
    "human_mab": 0.25,
    "enzyme": 0.55,
    "cytokine": 0.35,
    "hormone_peptide": 0.25,
    "gene_therapy": 0.50,
    "other": 0.35,
    "unknown": 0.35,
}

WEIGHTS = {
    "mhc1_burden": 0.26,
    "mhc2_burden": 0.14,
    "bcell_epitope": 0.16,
    "aggregation": 0.10,
    "liabilities": 0.09,
    "modality_prior": 0.13,
    "novelty": 0.12,
}

# component value at which the saturation (1.0) is reached
SATURATION = {
    # strong-binder-rate (8-11mer windows, <=500 nM, 8-allele panel) at the
    # p90 of a 10-drug FDA reference panel measured with the real mhcflurry
    # 2.2.1 predictor (rates 0.053-0.217, median 0.096, p90 ~0.20). The
    # component thus measures HLA-binding density vs an approved-drug
    # background; reproducible via tools/calibrate_mhc1_saturation.py.
    # This is a scale calibration, NOT a fit to clinical ADA incidence.
    "mhc1_burden": 0.20,
    "mhc2_burden": 0.15,   # fraction of scanned 15mers that are strong binders
    "bcell_epitope": 0.35,
    "aggregation": 0.25,
    "liabilities": 8.0,    # liability events per 100 residues
}


def _sat(value, cap):
    return max(0.0, min(1.0, value / cap))


def tier(total):
    if total < 30:
        return "low"
    if total < 50:
        return "low_moderate"
    if total < 70:
        return "moderate_high"
    return "high"


TIER_LABELS = {
    "low": "低风险",
    "low_moderate": "中低风险",
    "moderate_high": "中高风险",
    "high": "高风险",
}


# --------------------------------------------------------------------- #
# 临床暴露情境调整（定性经验规则层，非统计拟合）
#
# 证据锚点（均为定性结论，出处见方法页）：
#   * 给药途径：慢性生物制剂中，皮下自我注射在说明书中的 ADA 报告率常高于
#     静脉给药（如抗 TNF 类 SC 与 IV 制剂的标签数据趋势，FDA 2019 指南
#     将给药途径列为免疫原性决定因素）。
#   * 给药频次：高频/长期给药累积暴露增加，ADA 检出概率上升（指南定性结论）。
#   * 合并用药：合并免疫抑制剂（甲氨蝶呤/硫唑嘌呤）在多数单抗说明书数据中
#     一致性地降低 ADA 检出率。
#   * 人群：部分生物制剂在儿童人群中报告的 ADA 率偏高（数据有限，小幅调整）。
#
# 各项为有界序数分，刻意远小于序列分分项尺度；输入"未指定"时不产生调整。
# 明确不与临床 ADA 发生率做统计回归——没有逐患者数据集，不允许编造拟合。
CLINICAL_FACTORS = {
    "route":        {"iv": 0.0, "sc": 2.0, "im": 1.0, "unknown": 0.0},
    "frequency":    {"weekly_or_more": 2.0, "every_2_4_weeks": 1.0,
                     "monthly_or_less": 0.0, "unknown": 0.0},
    "comedication": {"none": 2.0, "immunosuppressant": 0.0, "unknown": 0.0},
    "population":   {"adult": 0.0, "pediatric": 1.0, "unknown": 0.0},
}

CLINICAL_LABELS = {
    "route": ("给药途径", {"iv": "静脉输注", "sc": "皮下注射", "im": "肌肉注射",
                           "unknown": "未指定"}),
    "frequency": ("给药频次", {"weekly_or_more": "每周或更频繁",
                               "every_2_4_weeks": "每 2–4 周",
                               "monthly_or_less": "每月或更少",
                               "unknown": "未指定"}),
    "comedication": ("合并用药", {"none": "无免疫抑制剂",
                                  "immunosuppressant": "合并免疫抑制剂",
                                  "unknown": "未指定"}),
    "population": ("目标人群", {"adult": "成人常规", "pediatric": "儿童人群",
                                "unknown": "未指定"}),
}


def clinical_context(ctx):
    """Qualitative clinical-exposure adjustment. Returns None when no factor
    is specified (or all resolve to "unknown"); otherwise per-factor rows and
    the summed bonus (capped at +8)."""
    if not ctx:
        return None
    rows = []
    total = 0.0
    for key, spec in CLINICAL_FACTORS.items():
        choice = str(ctx.get(key) or "unknown").lower()
        if choice not in spec:
            choice = "unknown"
        pts = spec[choice]
        total += pts
        label, choices = CLINICAL_LABELS[key]
        rows.append({"factor": key, "label": label,
                     "choice": choice,
                     "choice_label": choices.get(choice, choice),
                     "points": pts})
    if all(r["choice"] == "unknown" for r in rows):
        return None
    adj = round(min(total, 8.0), 1)
    return {"inputs": rows, "adjustment": adj}


def liability_rate(features, length):
    """Liability events per 100 residues (methionine weighted 0.5)."""
    events = (features["n_deamidation_NG_NS_NN"]
              + features["n_isomerization_DG_DS_DD"]
              + features["n_acid_clipping_DP"]
              + 0.5 * features["n_methionine"]
              + 0.3 * features["n_glycosylation"])
    return 100.0 * events / length if length else 0.0


def compute_score(features, mhc1_rate, mhc2_rate, bcell_cov,
                  exposed_fraction, modality, novelty):
    """Return (total 0..100, component breakdown).

    ``mhc1_rate`` / ``mhc2_rate``: strong-binder peptides per scanned window.
    ``bcell_cov``: residues covered by linear B-cell epitope segments.
    """
    length = features["length_total"]
    comps = {
        "mhc1_burden": _sat(mhc1_rate, SATURATION["mhc1_burden"]),
        "mhc2_burden": _sat(mhc2_rate, SATURATION["mhc2_burden"]),
        "bcell_epitope": _sat(bcell_cov / length if length else 0,
                              SATURATION["bcell_epitope"]),
        "aggregation": _sat(features["frac_hydrophobic_patch"],
                            SATURATION["aggregation"]),
        "liabilities": _sat(liability_rate(features, length),
                            SATURATION["liabilities"]),
        "modality_prior": MODALITY_PRIOR.get(modality, MODALITY_PRIOR["unknown"]),
    }
    if exposed_fraction is not None:
        # fold structural evidence into the B-cell component
        comps["bcell_epitope"] = min(1.0, 0.6 * comps["bcell_epitope"]
                                     + 0.4 * _sat(exposed_fraction, 0.5))
    if novelty is None:
        comps["novelty"] = None
    else:
        comps["novelty"] = max(0.0, min(1.0, novelty))

    active = {k: v for k, v in comps.items() if v is not None}
    wsum = sum(WEIGHTS[k] for k in active)
    total = sum(WEIGHTS[k] * v for k, v in active.items()) / wsum * 100.0

    breakdown = []
    labels = {
        "mhc1_burden": "MHC-I 限制性 T 细胞表位负荷",
        "mhc2_burden": "MHC-II 限制性辅助性 T 细胞表位",
        "bcell_epitope": "B 细胞线性/表面表位",
        "aggregation": "聚集倾向疏水区",
        "liabilities": "翻译后修饰/降解风险位点",
        "modality_prior": "药物模态先验（FDA 数据）",
        "novelty": "序列新颖性（相对已知药物）",
    }
    for key, val in comps.items():
        breakdown.append({
            "key": key,
            "label": labels[key],
            "value": None if val is None else round(val, 3),
            "weight": round(WEIGHTS[key] / wsum if key in active and wsum else 0, 3),
        })
    return round(total, 1), breakdown
