"""
Empirical context annotations: route of administration & patient population.

These notes do NOT change the main ADA / NAb scores. They summarize marketed-product
patterns as qualitative empiric deltas for interpretation only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .features import FeatureSet


@dataclass
class ContextNote:
    kind: str  # "route" | "population" | "general"
    title_zh: str
    title_en: str
    empiric_delta: str  # e.g. "+5~+15" or "−5~0"
    body_zh: str
    body_en: str
    examples_zh: str = ""
    examples_en: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "title_zh": self.title_zh,
            "title_en": self.title_en,
            "empiric_delta": self.empiric_delta,
            "body_zh": self.body_zh,
            "body_en": self.body_en,
            "examples_zh": self.examples_zh,
            "examples_en": self.examples_en,
        }


def _cls_blob(features: FeatureSet) -> str:
    parts: list[str] = []
    for n in features.neighbors[:5]:
        parts.append(f"{n.inn} {n.drug_class}")
    return " ".join(parts).lower()


def _inn_blob(features: FeatureSet) -> str:
    return " ".join((n.inn or "").lower() for n in features.neighbors[:5])


def _is_antibody_like(features: FeatureSet) -> bool:
    regions = getattr(features, "antibody_regions", None)
    return bool(regions is not None and getattr(regions, "is_antibody_like", False))


def build_context_notes(features: FeatureSet) -> list[ContextNote]:
    """Pick relevant route/population empiric notes from neighbor class."""
    blob = _cls_blob(features)
    inns = _inn_blob(features)
    notes: list[ContextNote] = []

    # One short scope line
    if _is_antibody_like(features):
        notes.append(
            ContextNote(
                kind="general",
                title_zh="序列分 vs 临床",
                title_en="Sequence vs clinical",
                empiric_delta="例外见下",
                body_zh="多数人源/人源化低 ADA 一致；自身免疫 SC、眼内 scFv、耗竭等例外见下。",
                body_en="Most human/humanized mAbs align low; autoimmunity/SC, ocular scFv, depletion below.",
            )
        )
    else:
        notes.append(
            ContextNote(
                kind="general",
                title_zh="序列分 vs 临床",
                title_en="Sequence vs clinical",
                empiric_delta="近内源低；ERT/异源酶高",
                body_zh="近内源激素/凝血偏低；ERT 缺失自身与病原酶偏高。",
                body_en="Near-endogenous hormones/coag low; ERT missing-self and pathogen enzymes high.",
            )
        )

    # Mechanism outliers (only when neighbor matches)
    if any(
        k in inns
        for k in (
            "adalimumab", "golimumab", "certolizumab", "infliximab",
            "brolucizumab", "ranibizumab", "aflibercept",
            "rituximab", "ofatumumab", "obinutuzumab", "ocrelizumab", "alemtuzumab",
        )
    ) or any(k in blob for k in ("anti-tnf", "tnf", "vegf", "cd20", "ophthalm", "scfv")):
        if any(k in inns for k in ("adalimumab", "golimumab", "certolizumab", "infliximab")) or "tnf" in blob:
            notes.append(
                ContextNote(
                    kind="population",
                    title_zh="抗 TNF：人源序列但临床可偏高",
                    title_en="Anti-TNF: human-like sequence, clinical may be high",
                    empiric_delta="表位易低估",
                    body_zh="长期 SC + 自身免疫人群；联用 MTX 常降低 ADA。",
                    body_en="Chronic SC + autoimmune setting; MTX combo often lowers ADA.",
                    examples_zh="Adalimumab",
                    examples_en="Adalimumab",
                )
            )
        if any(k in inns for k in ("brolucizumab", "ranibizumab", "aflibercept")) or any(
            k in blob for k in ("vegf", "ophthalm", "scfv", "intravitreal")
        ):
            notes.append(
                ContextNote(
                    kind="route",
                    title_zh="眼内 / scFv 格式",
                    title_en="Ocular / scFv format",
                    empiric_delta="格式暴露↑",
                    body_zh="scFv 暴露面计入异源度；玻璃体腔反复给药可再抬高临床检出。",
                    body_en="scFv exposure enters foreignness; repeated IVT may further raise clinical rates.",
                    examples_zh="Brolucizumab",
                    examples_en="Brolucizumab",
                )
            )
        if any(
            k in inns
            for k in ("rituximab", "ofatumumab", "obinutuzumab", "ocrelizumab", "alemtuzumab")
        ) or any(k in blob for k in ("cd20", "cd52", "deplet")):
            notes.append(
                ContextNote(
                    kind="population",
                    title_zh="耗竭型：测得 ADA 可偏低",
                    title_en="Depleting: measured ADA can look low",
                    empiric_delta="≠无免疫原性",
                    body_zh="B 细胞耗竭常压低测得结合 ADA；阿仑单抗等仍可很高，需个案。",
                    body_en="B-cell depletion often suppresses measured ADA; alemtuzumab can still be high.",
                    examples_zh="Rituximab vs Alemtuzumab",
                    examples_en="Rituximab vs alemtuzumab",
                )
            )

    # Route — only when class-relevant (IV baseline documented in glossary)
    if _is_antibody_like(features) or any(k in blob for k in ("fusion", "mab", "antibody")):
        notes.append(
            ContextNote(
                kind="route",
                title_zh="皮下 vs 静脉",
                title_en="SC vs IV",
                empiric_delta="+5~+20（相对 IV）",
                body_zh="同分子 SC 常高于 IV；主分不区分途径。",
                body_en="Same molecule often higher ADA on SC than IV; main score is route-agnostic.",
            )
        )

    if any(k in blob for k in ("insulin", "glp1", "peptide", "hormone-gh", "somatropin", "pramlintide")):
        notes.append(
            ContextNote(
                kind="route",
                title_zh="肽/激素反复皮下",
                title_en="Repeated SC peptide / hormone",
                empiric_delta="+0~+15",
                body_zh="慢性 SC 可提高结合抗体检出；近内源多为低滴度非中和。",
                body_en="Chronic SC can raise binding Ab; near-endogenous often low-titer non-NAb.",
            )
        )

    if any(k in blob for k in ("enzyme", "ert", "asparaginase", "urate")):
        notes.append(
            ContextNote(
                kind="route",
                title_zh="ERT / 外源酶",
                title_en="ERT / exogenous enzyme",
                empiric_delta="基线已高",
                body_zh="即使 IV，ADA/NAb 仍常高；基因型与既往暴露更关键。",
                body_en="Even on IV, ADA/NAb often high; genotype and prior exposure dominate.",
            )
        )

    if any(k in blob for k in ("interferon", "cytokine", "gcsf", "epo")):
        notes.append(
            ContextNote(
                kind="route",
                title_zh="细胞因子长期皮下",
                title_en="Long-term SC cytokine",
                empiric_delta="+5~+25",
                body_zh="IFN 等长期 SC 的 ADA/NAb 差异大；制剂与频率可调制。",
                body_en="Long-term SC IFN etc. vary widely; formulation and frequency modulate rates.",
            )
        )

    # Population — only when neighbor class matches
    if any(k in blob for k in ("tnf", "anti-tnf", "ra", "ibd", "autoimmun")) or any(
        k in inns for k in ("adalimumab", "infliximab", "golimumab", "certolizumab")
    ):
        notes.append(
            ContextNote(
                kind="population",
                title_zh="自身免疫人群 / 联用 MTX",
                title_en="Autoimmune setting / MTX combo",
                empiric_delta="+5~+25；联用 MTX −",
                body_zh="RA/IBD 背景下抗 TNF 常更高；甲氨蝶呤等可下调。",
                body_en="Anti-TNF often higher in RA/IBD; MTX etc. can lower ADA.",
            )
        )

    if any(k in blob for k in ("enzyme", "ert", "alglucosidase", "agalsidase", "laronidase", "idursulfase")):
        notes.append(
            ContextNote(
                kind="population",
                title_zh="CRIM 阴性 / 酶缺失",
                title_en="CRIM-negative / enzyme-null",
                empiric_delta="+20~+40",
                body_zh="无内源交叉耐受时高滴度 ADA/NAb 更常见。",
                body_en="Without endogenous cross-tolerance, high-titer ADA/NAb are more common.",
            )
        )

    if any(k in blob for k in ("factor", "coagulation", "octocog", "efmoroctocog", "nonacog", "porcine")):
        notes.append(
            ContextNote(
                kind="population",
                title_zh="凝血因子：抑制物另计",
                title_en="Clotting factors: inhibitors separate",
                empiric_delta="序列近内源",
                body_zh="人源 FVIII/FIX 结构分应低；临床抑制物多因缺失/频繁暴露/猪源。",
                body_en="Human FVIII/FIX stay structurally low; inhibitors often from deficiency/exposure/porcine.",
            )
        )

    if any(k in blob for k in ("hormone-gh", "somatropin", "insulin", "glp1", "peptide")):
        notes.append(
            ContextNote(
                kind="population",
                title_zh="儿科长期替代",
                title_en="Pediatric long-term use",
                empiric_delta="+0~+10",
                body_zh="结合抗体不少见，近内源激素多数无临床意义。",
                body_en="Binding Abs not rare; near-endogenous hormones often clinically silent.",
            )
        )

    if any(k in blob for k in ("asparaginase", "rasburicase", "pegloticase")):
        notes.append(
            ContextNote(
                kind="population",
                title_zh="再次暴露",
                title_en="Re-exposure",
                empiric_delta="+15~+35",
                body_zh="高免疫原性酶再次暴露时超敏与中和风险明显上升。",
                body_en="Re-exposure markedly raises hypersensitivity and neutralizing risk.",
            )
        )

    if any(
        k in blob or k in inns
        for k in (
            "imlifidase", "ides", "endos", "ricefidase", "bacterial", "streptococ",
            "streptokinase", "pegvaliase", "glucarpidase", "asparaginase",
        )
    ) or any(k in blob for k in ("enzyme-bacterial", "enzyme-xeno", "pathogen")):
        notes.append(
            ContextNote(
                kind="population",
                title_zh="病原酶 / 预存 ADA",
                title_en="Pathogen enzyme / pre-existing ADA",
                empiric_delta="预存↑↑",
                body_zh="常见病原体来源蛋白常有人群预存抗体；风险由暴露史主导，非单纯表位密度。",
                body_en="Common-pathogen proteins often have pre-existing Abs; exposure history dominates over epitope density.",
                examples_zh="Imlifidase / IdeS",
                examples_en="Imlifidase / IdeS",
            )
        )

    if features.neighbor_similarity < 0.5 or not features.neighbors:
        notes.append(
            ContextNote(
                kind="general",
                title_zh="无相近上市邻居",
                title_en="No close marketed neighbor",
                empiric_delta="不确定↑",
                body_zh="途径/人群注释外推价值下降，更依赖试验设计。",
                body_en="Route/population notes are weaker; lean on trial design.",
            )
        )

    cal = getattr(features, "ab_epitope_calibration", None) or {}
    flags = set(cal.get("risk_flags") or [])
    depleting_flag = any(
        k in inns
        for k in ("rituximab", "ofatumumab", "obinutuzumab", "ocrelizumab", "alemtuzumab", "ublituximab")
    )
    if (
        flags & {"hetero_cdr", "germline_hetero", "t_on_hetero", "b_on_hetero", "agg", "t_cdr"}
        and not depleting_flag
        and float(cal.get("risk_uplift") or 0) > 0
    ):
        notes.append(
            ContextNote(
                kind="general",
                title_zh="多信号上浮",
                title_en="Multi-signal uplift",
                empiric_delta=f"+{float(cal.get('risk_uplift') or 0):.0f}",
                body_zh="种系异源 / MHC·B 错配覆盖 / 聚集等触发时，表位可有限上浮。",
                body_en="Germline hetero / MHC·B-on-mismatch / agg can give a capped epitope uplift.",
            )
        )

    return notes
