"""
Neutralizing ADA (NAb) risk heuristics.

Clinical neutralizing activity is not the same as binding-ADA incidence.
This module combines:
  - molecular / modality foreignness (structure & origin proxies)
  - high binding-ADA prior (titer/incidence related risk)
  - literature/neighbor NAb annotations when available

Research heuristic only — not a substitute for cell-based NAb assays.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any
try:
    from typing import Literal
except ImportError:  # Python 3.7
    from typing_extensions import Literal

from .features import FeatureSet

NabLevel = Literal["low", "moderate", "high", "unknown"]

NAB_LABELS = {
    "low": {"zh": "低", "en": "Low"},
    "moderate": {"zh": "中", "en": "Moderate"},
    "high": {"zh": "高", "en": "High"},
    "unknown": {"zh": "未知", "en": "Unknown"},
}

NAB_RULES = [
    {
        "zh": "低中和等级：综合分 < 25（结构偏人源/内源样，且 ADA 先验不高）",
        "en": "Low NAb level: score < 25 (more human/endogenous-like; ADA prior not high).",
    },
    {
        "zh": "中中和等级：综合分 25–54.9",
        "en": "Moderate NAb level: score 25–54.9.",
    },
    {
        "zh": "高中和等级：综合分 ≥ 55（外源结构 + 高 ADA 发生率/文献中和提示）",
        "en": "High NAb level: score ≥ 55 (foreign scaffold + high ADA and/or literature NAb signal).",
    },
]


@dataclass
class NabScore:
    score: float
    nab_level: NabLevel
    nab_level_label: str
    nab_level_label_en: str
    structure_component: float
    ada_incidence_component: float
    literature_component: float | None
    explanation_items: list[dict[str, str]]
    level_rules_bi: list[dict[str, str]]

    def to_dict(self) -> dict[str, Any]:
        return {
            "score": self.score,
            "nab_level": self.nab_level,
            "nab_level_label": self.nab_level_label,
            "nab_level_label_en": self.nab_level_label_en,
            "structure_component": self.structure_component,
            "ada_incidence_component": self.ada_incidence_component,
            "literature_component": self.literature_component,
            "explanation_items": self.explanation_items,
            "level_rules_bi": self.level_rules_bi,
        }


def _clip(x: float, lo: float = 0.0, hi: float = 100.0) -> float:
    return max(lo, min(hi, x))


def _item(zh: str, en: str, note_zh: str, note_en: str) -> dict[str, str]:
    return {"zh": zh, "en": en, "note_zh": note_zh, "note_en": note_en}


def _structure_foreignness(features: FeatureSet) -> tuple[float, list[dict[str, str]]]:
    """
    Proxy for molecular foreignness from sequence homology / Ab format / neighbor class.
    Higher = more likely functional (neutralizing) antibodies matter clinically.
    """
    notes: list[dict[str, str]] = []
    score = 25.0  # baseline for recombinant proteins
    top = features.neighbors[0] if features.neighbors else None
    drug_class = (top.drug_class if top else "").lower()
    inn = (top.inn if top else "").lower()
    seq_signal = False

    ab = features.ab_epitope_calibration or {}
    pc = features.protein_epitope_calibration or {}
    if isinstance(ab, dict) and ab.get("format_exposure") is not None:
        fe = float(ab.get("format_exposure") or 0)
        fmt = ""
        if isinstance(ab.get("ab_format"), dict):
            fmt = str(ab["ab_format"].get("format") or "")
        if fe >= 0.55:
            score += 18
            seq_signal = True
            notes.append(
                _item(
                    zh=f"格式暴露面偏高（{fmt or 'fragment'}，{fe:.2f}）：预存/片段 neoepitope 可抬高中和相关风险。",
                    en=f"Elevated format exposure ({fmt or 'fragment'}, {fe:.2f}): preexisting/fragment neoepitopes can raise NAb-relevant risk.",
                    note_zh="来自序列架构检测，不是晶体对接。",
                    note_en="From sequence architecture detection — not docking.",
                )
            )
    if isinstance(pc, dict) and pc.get("homology"):
        hom = pc["homology"] or {}
        hf = float(hom.get("foreignness") or 0)
        if hf >= 70:
            score += 20
            seq_signal = True
            notes.append(
                _item(
                    zh=f"对人源内源库异源度高（≈{hf:.0f}）：功能中和更需关注。",
                    en=f"High foreignness vs human self pool (≈{hf:.0f}): functional neutralization more relevant.",
                    note_zh="序列比对异源度，非类别标签。",
                    note_en="Sequence homology foreignness, not a class label alone.",
                )
            )
        elif hf <= 15:
            score -= 12
            seq_signal = True
            notes.append(
                _item(
                    zh=f"对人源内源库高度相似（异源≈{hf:.0f}）：中和风险通常较低。",
                    en=f"Near-self vs human pool (foreignness≈{hf:.0f}): neutralizing risk usually lower.",
                    note_zh="ERT 等例外见耐受注释。",
                    note_en="ERT-like exceptions covered by tolerance notes.",
                )
            )

    rules = [
        (
            any(k in drug_class for k in ("enzyme-ert", "enzyme", "asparaginase", "urate")),
            35,
            "酶替代/外源酶类别：中和抗体常直接削弱疗效或引起超敏。",
            "ERT / exogenous enzymes: NAbs often blunt efficacy or drive hypersensitivity.",
        ),
        (
            any(k in drug_class for k in ("peptide-glp1", "peptide", "exendin"))
            or "exenatide" in inn,
            20,
            "非人源肽（如 Exendin-4 类）：结合 ADA 常见，部分可中和生物活性。",
            "Non-human peptides (e.g. exendin-4 class): binding ADA common; a subset can neutralize activity.",
        ),
        (
            "chimeric" in drug_class or "murine" in drug_class,
            25,
            "嵌合/鼠源骨架：抗独特型或抗恒定区抗体更易影响功能。",
            "Chimeric/murine scaffolds: anti-idiotype/constant-region Abs more often affect function.",
        ),
        (
            any(k in drug_class for k in ("bacterial", "toxin", "immunotoxin", "streptokinase", "ides", "imlifidase"))
            or any(k in inn for k in ("imlifidase", "streptokinase", "pegvaliase")),
            40,
            "细菌/毒素来源蛋白：中和与超敏风险通常很高。",
            "Bacterial/toxin-derived proteins: neutralizing and hypersensitivity risk usually high.",
        ),
        (
            any(k in drug_class for k in ("coagulation", "factor viii", "factor-")),
            20,
            "凝血因子：抑制物（中和性 ADA）是关键临床终点。",
            "Clotting factors: inhibitors (neutralizing ADA) are a key clinical endpoint.",
        ),
        (
            any(k in drug_class for k in ("humanized", "fully human", "whole mab", "monoclonal"))
            and "chimeric" not in drug_class,
            -10,
            "人源化/全人抗体：中和性 ADA 相对少见，但仍取决于表位与滴度。",
            "Humanized/fully human mAbs: neutralizing ADA less common, still epitope/titer dependent.",
        ),
        (
            any(k in drug_class for k in ("insulin", "hormone-gh", "somatropin", "human")),
            -15,
            "接近内源激素序列：结合 ADA 多为低滴度、非中和。",
            "Near-endogenous hormones: binding ADA often low-titer and non-neutralizing.",
        ),
    ]

    applied = False
    if not seq_signal:
        for cond, delta, zh, en in rules:
            if cond:
                score += delta
                applied = True
                notes.append(
                    _item(
                        zh=f"结构/类别因子: {zh}",
                        en=f"Structure/class factor: {en}",
                        note_zh="根据最邻近上市药的分子类别做外源性代理，不是晶体结构对接。",
                        note_en="Foreignness proxy from nearest marketed product class — not a docking model.",
                    )
                )
                break
        if not applied:
            notes.append(
                _item(
                    zh="结构/类别因子: 未命中强规则，使用基线外源性分。",
                    en="Structure/class factor: no strong rule hit; using baseline foreignness.",
                    note_zh="缺少明确类别邻居时，中和等级更依赖 ADA 发生率与表位负担。",
                    note_en="Without a clear class neighbor, NAb level leans on ADA incidence and epitope burden.",
                )
            )

    if features.t_epitope_density >= 15:
        score += 8
        notes.append(
            _item(
                zh="MHC-II 表位密度偏高，可能支持更强的 T 依赖应答与更高滴度。",
                en="Elevated MHC-II epitope density may support stronger T-dependent responses and higher titers.",
                note_zh="滴度升高时中和活性概率上升（经验相关，非因果证明）。",
                note_en="Higher titers correlate empirically with neutralizing activity (not proof of causality).",
            )
        )

    return _clip(score), notes


def score_nab(features: FeatureSet, ada_score: float) -> NabScore:
    items: list[dict[str, str]] = []
    structure, struct_notes = _structure_foreignness(features)
    items.extend(struct_notes)

    # Binding ADA incidence / composite score as titer-related proxy
    ada_prior = features.neighbor_ada_pct
    if ada_prior is not None and features.neighbor_trust > 0:
        ada_comp = _clip(float(ada_prior))
        items.append(
            _item(
                zh=f"结合 ADA 发生率先验={ada_comp:.1f}%（用于近似滴度/暴露相关风险）",
                en=f"Binding-ADA incidence prior={ada_comp:.1f}% (proxy for titer/exposure-related risk)",
                note_zh="高发生率常伴随更高滴度亚群；中和更常见于高滴度，但并非一一对应。",
                note_en="High incidence often coexists with higher-titer subsets; neutralization is more common at high titer, not 1:1.",
            )
        )
    else:
        ada_comp = _clip(ada_score)
        items.append(
            _item(
                zh=f"无临床邻居时，用综合 ADA 分={ada_comp:.1f} 作为发生率代理。",
                en=f"No clinical neighbor; using composite ADA score={ada_comp:.1f} as incidence proxy.",
                note_zh="缺少文献发生率时不确定性更大。",
                note_en="Uncertainty is higher without literature incidence.",
            )
        )

    # Literature / curated neighbor NAb signal if present on FeatureSet
    lit = getattr(features, "neighbor_nab_pct", None)
    lit_note_en = getattr(features, "neighbor_nab_note", "") or ""
    lit_note_zh = getattr(features, "neighbor_nab_note_zh", "") or ""
    if lit is not None and features.neighbor_trust > 0:
        lit_comp = _clip(float(lit))
        items.append(
            _item(
                zh=f"文献/注释中和提示≈{lit_comp:.1f}% {lit_note_zh}".strip(),
                en=f"Literature/annotated NAb signal≈{lit_comp:.1f}% {lit_note_en}".strip(),
                note_zh="来自邻近上市药的中和 ADA 注释；覆盖不全。",
                note_en="From neighbor marketed-product NAb annotations; coverage is incomplete.",
            )
        )
    else:
        lit_comp = None
        items.append(
            _item(
                zh="无直接中和 ADA 文献数值；仅用结构外源性 + 结合 ADA 先验推断。",
                en="No direct NAb literature value; inferring from foreignness + binding-ADA prior.",
                note_zh="建议结合细胞法 NAb 实验与滴度数据解读。",
                note_en="Interpret alongside cell-based NAb assays and titer data.",
            )
        )

    if lit_comp is None:
        score = 0.45 * structure + 0.55 * ada_comp
    else:
        score = 0.35 * structure + 0.35 * ada_comp + 0.30 * lit_comp

    score = _clip(score)
    level = _nab_level(score)
    items.append(
        _item(
            zh=f"中和等级综合分={score:.1f} → {NAB_LABELS[level]['zh']} 中和等级",
            en=f"NAb level score={score:.1f} → {NAB_LABELS[level]['en']} NAb level",
            note_zh="预测的是中和可能性倾向，不是确认的中和滴度。",
            note_en="This is a tendency for neutralizing potential, not a confirmed neutralizing titer.",
        )
    )

    return NabScore(
        score=round(score, 2),
        nab_level=level,
        nab_level_label=NAB_LABELS[level]["zh"],
        nab_level_label_en=NAB_LABELS[level]["en"],
        structure_component=round(structure, 2),
        ada_incidence_component=round(ada_comp, 2),
        literature_component=None if lit_comp is None else round(lit_comp, 2),
        explanation_items=items,
        level_rules_bi=list(NAB_RULES),
    )


def _nab_level(score: float) -> NabLevel:
    if score < 25:
        return "low"
    if score < 55:
        return "moderate"
    return "high"
