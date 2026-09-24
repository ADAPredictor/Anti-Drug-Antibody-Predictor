"""
Interpretable ADA-level scoring model.

Combines sequence-derived T/B epitope features with clinical ADA priors from
similar marketed biologics. This is a transparent research heuristic, not a
regulatory-grade predictor.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict
try:
    from typing import Literal
except ImportError:  # Python 3.7
    from typing_extensions import Literal

from .features import FeatureSet

ADALevel = Literal["low", "moderate", "high", "unknown"]

ADA_LEVEL_LABELS_ZH = {
    "low": "低ADA",
    "moderate": "中ADA",
    "high": "高ADA",
    "unknown": "未知",
}
ADA_LEVEL_LABELS_EN = {
    "low": "Low ADA",
    "moderate": "Moderate ADA",
    "high": "High ADA",
    "unknown": "Unknown",
}

ADA_LEVEL_RULES_BI = [
    {
        "zh": "低ADA：综合分 < 20（约对应临床 ADA 发生率偏低）",
        "en": "Low ADA: composite score < 20 (typically lower clinical ADA incidence).",
    },
    {
        "zh": "中ADA：综合分 20–44.9",
        "en": "Moderate ADA: composite score 20–44.9.",
    },
    {
        "zh": "高ADA：综合分 ≥ 45（约对应临床 ADA 发生率偏高）",
        "en": "High ADA: composite score ≥ 45 (typically higher clinical ADA incidence).",
    },
]

ADA_LEVEL_LABELS = ADA_LEVEL_LABELS_ZH
ADA_LEVEL_RULES = tuple(r["zh"] for r in ADA_LEVEL_RULES_BI)


@dataclass
class RiskScore:
    score: float
    tier: ADALevel
    ada_level: ADALevel
    ada_level_label: str
    ada_level_label_en: str
    epitope_component: float
    ada_prior_component: float | None
    blend_weight: float | None
    explanation: list[str]
    explanation_items: list[dict[str, str]]
    level_rules: list[str]
    level_rules_bi: list[dict[str, str]]
    breakdown: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "score": self.score,
            "tier": self.tier,
            "ada_level": self.ada_level,
            "ada_level_label": self.ada_level_label,
            "ada_level_label_en": self.ada_level_label_en,
            "epitope_component": self.epitope_component,
            "ada_prior_component": self.ada_prior_component,
            "blend_weight": self.blend_weight,
            "explanation": self.explanation,
            "explanation_items": self.explanation_items,
            "level_rules": self.level_rules,
            "level_rules_bi": self.level_rules_bi,
            "breakdown": self.breakdown,
        }


def _clip(x: float, lo: float = 0.0, hi: float = 100.0) -> float:
    return max(lo, min(hi, x))


def _item(zh: str, en: str, note_zh: str, note_en: str) -> dict[str, str]:
    return {"zh": zh, "en": en, "note_zh": note_zh, "note_en": note_en}


def _is_near_endogenous(features: FeatureSet) -> bool:
    if not features.neighbors or features.neighbor_similarity < 0.85:
        return False
    for top in features.neighbors[:3]:
        if top.similarity < 0.85:
            continue
        cls = (top.drug_class or "").lower()
        inn = (top.inn or "").lower()
        if any(
            k in cls or k in inn
            for k in ("exendin", "exenatide", "salmon", "toxin", "bacterial")
        ):
            continue
        if any(k in cls for k in ("antagonist", "pegylat", "fusion", "albumin")):
            continue
        if any(
            k in cls
            for k in (
                "hormone-gh",
                "hormone-insulin",
                "insulin",
                "hormone-pth",
                "hormone-igf",
                "somatropin",
                "glucagon",
            )
        ):
            return True
    return False


def _inn_antibody_lineage(inn: str) -> str:
    """WHO-style suffix: chimeric (-ximab), humanized (-zumab), human (-umab)."""
    token = inn.lower().replace("-", " ").split()[0]
    if token.endswith("ximab"):
        return "chimeric"
    if token.endswith("zumab"):
        return "humanized"
    if token.endswith("umab"):
        return "human"
    return "unknown"


def _is_humanlike_mab(features: FeatureSet) -> bool:
    """Humanized / human mAbs: germline frameworks inflate heuristic MHC density."""
    if not features.neighbors or features.neighbor_similarity < 0.92:
        return False
    for top in features.neighbors[:3]:
        if top.method not in {"vh", "vl", "fv"} or top.similarity < 0.92:
            continue
        cls = (top.drug_class or "").lower()
        inn = (top.inn or "").lower()
        lineage = _inn_antibody_lineage(inn)
        if lineage == "chimeric" or any(k in cls for k in ("chimeric", "murine", "mouse")):
            continue
        if lineage in {"humanized", "human"}:
            return True
        if any(
            k in cls
            for k in (
                "monoclonal",
                "antibody",
                "humanized",
                "human",
                "igg",
                "fab",
                "mab",
                "bispecific",
                "adc",
                "conjugate",
            )
        ):
            return True
    return False


def _epitope_scale(features: FeatureSet) -> tuple[float, str]:
    """
    Return (scale, reason) for dampening raw epitope component.
    reason in {endogenous, humanlike, none}
    Clinical ADA is never used here.
    """
    if _is_near_endogenous(features):
        return 0.35, "endogenous"
    if _is_humanlike_mab(features):
        return 0.40, "humanlike"
    return 1.0, "none"


def score_risk(features: FeatureSet) -> RiskScore:
    """Weighted ADA-oriented score in [0, 100] with bilingual annotated notes."""
    items: list[dict[str, str]] = []

    t_density_term = _clip(features.t_epitope_density / 28.0 * 100.0)
    t_score_term = _clip(features.t_mean_top_score / 2.5 * 100.0)
    b_cov_term = _clip(features.b_coverage * 100.0)
    mhc1_term = _clip(getattr(features, "mhc1_epitope_density", 0.0) / 40.0 * 100.0)
    phys = getattr(features, "physchem", None)
    agg_term = 0.0 if phys is None else _clip(phys.agg_risk_hint)

    full_seq_epitope = (
        0.45 * t_density_term
        + 0.20 * t_score_term
        + 0.20 * b_cov_term
        + 0.08 * mhc1_term
        + 0.07 * agg_term
    )

    # Transparent weighted raw terms (0-100 each) + their fixed weights.
    # These compose full_seq_epitope; for antibodies it serves as the
    # uncalibrated reference before CDR/germline/lineage calibration.
    breakdown: Dict[str, Any] = {
        "path": "heuristic",
        "raw_terms": [
            {"key": "mhc2_density", "label_zh": "MHC-II 表位密度（CD4）",
             "label_en": "MHC-II epitope density (CD4)",
             "weight": 0.45, "term": round(t_density_term, 2)},
            {"key": "mhc2_score", "label_zh": "MHC-II 强结合分值",
             "label_en": "MHC-II strong-binding score",
             "weight": 0.20, "term": round(t_score_term, 2)},
            {"key": "b_cell", "label_zh": "B 细胞表位覆盖",
             "label_en": "B-cell epitope coverage",
             "weight": 0.20, "term": round(b_cov_term, 2)},
            {"key": "mhc1_density", "label_zh": "MHC-I 表位密度（CD8）",
             "label_en": "MHC-I epitope density (CD8)",
             "weight": 0.08, "term": round(mhc1_term, 2)},
            {"key": "aggregation", "label_zh": "聚集倾向",
             "label_en": "Aggregation propensity",
             "weight": 0.07, "term": round(agg_term, 2)},
        ],
        "full_seq_epitope": round(full_seq_epitope, 2),
    }

    regions = getattr(features, "antibody_regions", None)
    use_cdr = bool(
        regions is not None
        and features.ab_epitope_calibration is not None
        and features.early_screen_epitope is not None
    )
    use_protein = bool(
        not use_cdr
        and features.protein_epitope_calibration is not None
        and features.early_screen_epitope is not None
    )

    if use_cdr:
        # Predictive Ab: foreignness × T/B × lineage tolerance (no clinical kNN)
        cal = getattr(features, "ab_epitope_calibration", None) or {}
        epitope_component = _clip(float(features.early_screen_epitope))
        breakdown["path"] = "antibody"
        breakdown["base_score"] = round(float(cal.get("base_score", epitope_component)), 2)
        breakdown["risk_uplift"] = round(float(cal.get("risk_uplift") or 0.0), 2)
        breakdown["lineage_inferred"] = cal.get("lineage_inferred", "unknown")
        breakdown["foreignness"] = round(float(cal.get("knn_clinical", 0.0)), 1)
        breakdown["risk_flags_zh"] = list(cal.get("flag_labels_zh") or [])
        breakdown["risk_flags_en"] = list(cal.get("flag_labels_en") or [])
        breakdown["pre_context_score"] = round(
            float(cal.get("pre_context_score") or 0.0), 2)
        breakdown["context_mult"] = float(cal.get("context_mult") or 1.0)
        breakdown["context_cap"] = cal.get("context_cap")
        breakdown["context_note_zh"] = str(cal.get("context_note_zh") or "")
        breakdown["context_note_en"] = str(cal.get("context_note_en") or "")
        breakdown["context_tags"] = list(cal.get("context_tags") or [])

        items.append(
            _item(
                zh=(
                    f"预测·抗体表位={epitope_component:.1f}"
                    f"（异源度≈{float(cal.get('knn_clinical', 0)):.0f}，"
                    f"推断人源化={cal.get('lineage_inferred', '?')}，"
                    f"框架人源距离={float(cal.get('fr_mismatch', 0)):.2f}，"
                    f"非种系负担={float(features.non_germline_burden or 0):.1f}）"
                ),
                en=(
                    f"Predictive Ab epitope={epitope_component:.1f}"
                    f" (foreignness≈{float(cal.get('knn_clinical', 0)):.0f}, "
                    f"inferred lineage={cal.get('lineage_inferred', '?')}, "
                    f"FR-human distance={float(cal.get('fr_mismatch', 0)):.2f}, "
                    f"non-germline burden={float(features.non_germline_burden or 0):.1f})"
                ),
                note_zh=str(
                    cal.get("note_zh")
                    or "异源度→T/B 表位→人源化耐受；临床 ADA 不参与矫正。"
                ),
                note_en=str(
                    cal.get("note_en")
                    or "Foreignness → T/B epitopes → lineage tolerance; clinical ADA not used to correct."
                ),
            )
        )
        items.append(
            _item(
                zh=(
                    f"分解：人源化项={float(cal.get('humanness_term', 0)):.1f}，"
                    f"负担项={float(cal.get('burden_term', 0)):.1f}，"
                    f"电荷项={float(cal.get('charge_term', 0)):.1f}；"
                    f"对照全长 Heuristic={full_seq_epitope:.1f}"
                ),
                en=(
                    f"Breakdown: humanness={float(cal.get('humanness_term', 0)):.1f}, "
                    f"burden={float(cal.get('burden_term', 0)):.1f}, "
                    f"charge={float(cal.get('charge_term', 0)):.1f}; "
                    f"full-length heuristic ref={full_seq_epitope:.1f}"
                ),
                note_zh=regions.method_note_zh if regions else "",
                note_en=regions.method_note_en if regions else "",
            )
        )
        if regions and regions.to_dict().get("cdrs"):
            cdr_txt = ", ".join(
                f"{c['name']}({c['abs_start']}-{c['abs_end']})"
                for c in regions.to_dict()["cdrs"][:6]
            )
            items.append(
                _item(
                    zh=f"识别到的 CDR：{cdr_txt}",
                    en=f"Detected CDRs: {cdr_txt}",
                    note_zh="Heuristic 锚点划分；正式申报请用 IMGT/ANARCI 编号复核。",
                    note_en="Heuristic anchors; verify with IMGT/ANARCI for formal filings.",
                )
            )
        uplift = float(cal.get("risk_uplift") or 0.0)
        flags = list(cal.get("risk_flags") or [])
        if uplift > 0 or flags:
            labels_zh = cal.get("flag_labels_zh") or flags
            labels_en = cal.get("flag_labels_en") or flags
            items.append(
                _item(
                    zh=(
                        f"多信号风险："
                        f"CDR异源性代理={float(cal.get('cdr_foreignness', 0)):.2f}，"
                        f"flags=[{', '.join(labels_zh) or '无'}]，"
                        f"上浮=+{uplift:.0f}"
                        f"（基线={float(cal.get('base_score', epitope_component)):.1f}→"
                        f"{epitope_component:.1f}）"
                    ),
                    en=(
                        f"Multi-signal risk: "
                        f"CDR foreignness proxy={float(cal.get('cdr_foreignness', 0)):.2f}, "
                        f"flags=[{', '.join(labels_en) or 'none'}], "
                        f"uplift=+{uplift:.0f}"
                        f" (base={float(cal.get('base_score', epitope_component)):.1f}→"
                        f"{epitope_component:.1f})"
                    ),
                    note_zh=str(
                        cal.get("risk_note_zh")
                        or "需多项独立信号交叉印证才上浮；耗竭型机制不加分。"
                    ),
                    note_en=str(
                        cal.get("risk_note_en")
                        or "Uplift requires multiple corroborating signals; depleting mAbs get no uplift."
                    ),
                )
            )
        gl = cal.get("germline") or getattr(features, "germline_report", None) or {}
        if gl and gl.get("hits"):
            hit_txt = ", ".join(
                f"{h.get('chain','?')}:{h.get('v_gene','?')}({float(h.get('identity', 0)):.0%},"
                f"mm={h.get('n_mismatch', 0)}/{h.get('n_aligned', 0)})"
                for h in gl["hits"][:4]
            )
            items.append(
                _item(
                    zh=(
                        f"种系 V 最佳命中：{hit_txt}；"
                        f"异源分率={float(gl.get('hetero_frac', 0)):.2f}"
                        f"（CDR={float(gl.get('cdr_hetero_frac', 0)):.2f}，"
                        f"FR={float(gl.get('fr_hetero_frac', 0)):.2f}），"
                        f"错配={int(gl.get('n_mismatch', 0))}"
                        f"（CDR={int(gl.get('n_cdr_mismatch', 0))}，"
                        f"FR={int(gl.get('n_fr_mismatch', 0))}），"
                        f"t_on_hetero={int(cal.get('t_on_hetero_count', 0))}"
                    ),
                    en=(
                        f"Germline V best-hit: {hit_txt}; "
                        f"hetero frac={float(gl.get('hetero_frac', 0)):.2f}"
                        f" (CDR={float(gl.get('cdr_hetero_frac', 0)):.2f}, "
                        f"FR={float(gl.get('fr_hetero_frac', 0)):.2f}), "
                        f"mismatches={int(gl.get('n_mismatch', 0))}"
                        f" (CDR={int(gl.get('n_cdr_mismatch', 0))}, "
                        f"FR={int(gl.get('n_fr_mismatch', 0))}), "
                        f"t_on_hetero={int(cal.get('t_on_hetero_count', 0))}"
                    ),
                    note_zh=str(
                        gl.get("method_note_zh")
                        or "Biopython 最佳命中 vs 捆绑人源 V 基因；非 ANARCI/HMMER。"
                    ),
                    note_en=str(
                        gl.get("method_note_en")
                        or "Biopython best-hit vs bundled human V genes; not ANARCI/HMMER."
                    ),
                )
            )
    elif use_protein:
        pcal = features.protein_epitope_calibration or {}
        epitope_component = _clip(float(features.early_screen_epitope))
        breakdown["path"] = "protein"
        breakdown["protein_class"] = pcal.get("protein_class", "unknown")
        breakdown["foreignness"] = round(
            float(pcal.get("foreignness", pcal.get("class_prior", 0))), 1)
        breakdown["epitope_raw"] = round(
            float(pcal.get("epitope_raw", pcal.get("raw_epitope", full_seq_epitope))), 2)
        breakdown["tolerance_mult"] = round(float(pcal.get("tolerance_mult", 1)), 3)
        breakdown["gate"] = round(float(pcal.get("gate", 1)), 3)
        items.append(
            _item(
                zh=(
                    f"预测·蛋白药={epitope_component:.1f}"
                    f"（类别={pcal.get('protein_class', '?')}，"
                    f"异源度={float(pcal.get('foreignness', pcal.get('class_prior', 0))):.0f}，"
                    f"T/B表位={float(pcal.get('epitope_raw', pcal.get('raw_epitope', full_seq_epitope))):.1f}，"
                    f"耐受系数={float(pcal.get('tolerance_mult', 1)):.2f}，"
                    f"门控={float(pcal.get('gate', 1)):.2f}）"
                ),
                en=(
                    f"Predictive protein={epitope_component:.1f}"
                    f" (class={pcal.get('protein_class', '?')}, "
                    f"foreignness={float(pcal.get('foreignness', pcal.get('class_prior', 0))):.0f}, "
                    f"T/B epitope={float(pcal.get('epitope_raw', pcal.get('raw_epitope', full_seq_epitope))):.1f}, "
                    f"tolerance={float(pcal.get('tolerance_mult', 1)):.2f}, "
                    f"gate={float(pcal.get('gate', 1)):.2f})"
                ),
                note_zh=str(
                    pcal.get("note_zh")
                    or "异源度→T/B 表位→产品类型耐受；临床 ADA 不参与矫正。"
                ),
                note_en=str(
                    pcal.get("note_en")
                    or "Foreignness → T/B epitopes → product-type tolerance; clinical ADA not used to correct."
                ),
            )
        )
    else:
        epitope_raw = full_seq_epitope
        epitope_component = epitope_raw
        scale, scale_reason = _epitope_scale(features)
        breakdown["tolerance_scale"] = scale
        breakdown["tolerance_reason"] = scale_reason
        # Do not use clinical self-match dampening
        if scale_reason == "endogenous":
            epitope_component = _clip(epitope_raw * scale)
            items.append(
                _item(
                    zh=(
                        f"近内源序列：表位分量从 {epitope_raw:.1f} 下调至 "
                        f"{epitope_component:.1f}（耐受校正）"
                    ),
                    en=(
                        f"Near-endogenous: epitope {epitope_raw:.1f} → "
                        f"{epitope_component:.1f} (tolerance correction)"
                    ),
                    note_zh="近内源：全长 MHC 密度易虚高，按耐受下调。",
                    note_en="Near-endogenous: full-length MHC often inflated; dampened by tolerance.",
                )
            )
        elif scale_reason == "humanlike":
            epitope_component = _clip(epitope_raw * scale)
            items.append(
                _item(
                    zh=(
                        f"人源化/全人抗体后备校正：表位 {epitope_raw:.1f} → "
                        f"{epitope_component:.1f}"
                    ),
                    en=(
                        f"Humanlike mAb fallback: epitope {epitope_raw:.1f} → "
                        f"{epitope_component:.1f}"
                    ),
                    note_zh="未走完整 CDR 标定路径时的耐受校正。",
                    note_en="Tolerance correction when full CDR path was not used.",
                )
            )
        else:
            items.append(
                _item(
                    zh=f"全长表位 Heuristic={epitope_component:.1f}（未走抗体/蛋白专用路径）",
                    en=f"Full-length epitope heuristic={epitope_component:.1f} (no Ab/protein path)",
                    note_zh="后备路径。",
                    note_en="Fallback path.",
                )
            )

    if phys is not None:
        items.append(
            _item(
                zh=(
                    f"物化性质：疏水比例={phys.hydrophobic_frac:.1%}，"
                    f"电荷密度={phys.charge_density:+.3f}，"
                    f"聚集相关提示={phys.agg_risk_hint:.1f}"
                ),
                en=(
                    f"Physicochemistry: hydrophobic frac={phys.hydrophobic_frac:.1%}, "
                    f"charge density={phys.charge_density:+.3f}, "
                    f"aggregation-related hint={phys.agg_risk_hint:.1f}"
                ),
                note_zh="疏水/电荷为次要调制，不能替代制剂开发数据。",
                note_en="Hydrophobicity/charge are minor modulators, not developability data.",
            )
        )

    items.append(
        _item(
            zh=(
                f"MHC-II (CD4) 强表位密度={features.t_epitope_density:.2f}/100aa "
                f"(n={features.t_strong_count})"
            ),
            en=(
                f"MHC-II (CD4) strong-epitope density={features.t_epitope_density:.2f}/100aa "
                f"(n={features.t_strong_count})"
            ),
            note_zh="全长密度；抗体请优先看上方 CDR/非种系负担。",
            note_en="Full-length density; for antibodies prefer CDR/non-germline burden above.",
        )
    )
    items.append(
        _item(
            zh=(
                f"MHC-I (CD8) 强表位密度="
                f"{getattr(features, 'mhc1_epitope_density', 0.0):.2f}/100aa "
                f"(n={getattr(features, 'mhc1_strong_count', 0)})"
            ),
            en=(
                f"MHC-I (CD8) strong-epitope density="
                f"{getattr(features, 'mhc1_epitope_density', 0.0):.2f}/100aa "
                f"(n={getattr(features, 'mhc1_strong_count', 0)})"
            ),
            note_zh="对 ADA 为参考。",
            note_en="Secondary for ADA.",
        )
    )
    items.append(
        _item(
            zh=(
                f"B 细胞线性表位覆盖率={features.b_coverage:.1%} "
                f"(n={features.b_epitope_count})"
            ),
            en=(
                f"B-cell linear epitope coverage={features.b_coverage:.1%} "
                f"(n={features.b_epitope_count})"
            ),
            note_zh="线性 B 表位 Heuristic；真实抗体表位多为构象型。",
            note_en="Linear B-cell heuristic; most Ab epitopes are conformational.",
        )
    )

    # Clinical neighbors are display-only (UI neighbor list); never fused into the score
    score = epitope_component
    breakdown["epitope_component"] = round(epitope_component, 2)
    ada_component: float | None = None
    blend_w: float | None = None
    if features.neighbor_ada_pct is not None and features.neighbor_trust > 0:
        ada_component = _clip(float(features.neighbor_ada_pct))

    score = _clip(score)
    level = _ada_level(score)
    label_zh = ADA_LEVEL_LABELS_ZH[level]
    label_en = ADA_LEVEL_LABELS_EN[level]
    items.append(
        _item(
            zh=f"综合分={score:.1f} → {label_zh}",
            en=f"Composite={score:.1f} → {label_en}",
            note_zh="分档见「名词解释」。",
            note_en="Bands in Glossary.",
        )
    )

    return RiskScore(
        score=round(score, 2),
        tier=level,
        ada_level=level,
        ada_level_label=label_zh,
        ada_level_label_en=label_en,
        epitope_component=round(epitope_component, 2),
        ada_prior_component=None if ada_component is None else round(ada_component, 2),
        blend_weight=None if blend_w is None else round(blend_w, 2),
        explanation=[it["zh"] for it in items],
        explanation_items=items,
        level_rules=[r["zh"] for r in ADA_LEVEL_RULES_BI],
        level_rules_bi=list(ADA_LEVEL_RULES_BI),
        breakdown=breakdown,
    )


def _ada_level(score: float) -> ADALevel:
    if score < 20:
        return "low"
    if score < 45:
        return "moderate"
    return "high"
