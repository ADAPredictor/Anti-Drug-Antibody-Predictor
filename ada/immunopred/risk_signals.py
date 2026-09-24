"""
Multi-signal risk flags and capped uplift on calibrated Ab epitope scores.

MHC / CDR density alone stays dampened unless several independent risk signals
fire. Depleting mechanisms (rituximab-class) get flags/notes only — no uplift.

When a GermlineMismatchReport is available, prefer germline_hetero / t_on_hetero /
b_on_hetero over the k-mer hetero_cdr proxy.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .germline_align import (
    GermlineMismatchReport,
    epitope_overlap_mismatches,
    span_overlap_mismatches,
)
from .germline_proxy import cdr_foreignness as compute_cdr_foreignness
from .germline_proxy import is_b_cell_depleting_cd19_20
from .t_cell import PEPTIDE_LENGTH


_FLAG_LABELS: dict[str, tuple[str, str]] = {
    "t_cdr": ("CDR 非人源种系表位偏高；", "Elevated CDR non-human-germline epitope burden;"),
    "b_epi": ("B 细胞线性表位覆盖偏高；", "Elevated B-cell linear epitope coverage;"),
    "hetero_cdr": (
        "CDR 异源性代理偏高（非 IMGT 种系）；",
        "Elevated CDR foreignness proxy (not IMGT germline);",
    ),
    "germline_hetero": (
        "种系 V 基因异源残基偏高；",
        "Elevated heterologous residues vs best-hit human V germline;",
    ),
    "t_on_hetero": (
        "MHC-II 肽覆盖种系异源 CDR 错配；",
        "MHC-II peptides overlap germline CDR mismatch residues;",
    ),
    "b_on_hetero": (
        "B 表位覆盖种系 CDR 错配；",
        "B-cell epitopes overlap germline CDR mismatch residues;",
    ),
    "agg": ("聚集相关风险提示偏高；", "Elevated aggregation-related risk hint;"),
    "fr_distant": ("框架远离低 ADA 人源库；", "Framework distant from low-ADA human pool;"),
    "charge": ("电荷密度偏高；", "Elevated absolute charge density;"),
}

# Prefer FR scaffold foreignness; CDR somatic mutations alone should not flag high risk
_GERMLINE_HETERO_FRAC = 0.35
_GERMLINE_FR_HETERO_FRAC = 0.22
_GERMLINE_CDR_HETERO_FRAC = 0.55
_T_ON_HETERO_MIN = 3
_B_ON_HETERO_MIN = 1


@dataclass
class RiskSignalReport:
    uplift: float  # 0..15
    flags: list[str]  # machine ids
    flag_labels_zh: list[str]
    flag_labels_en: list[str]
    cdr_foreignness: float
    depleting: bool
    note_zh: str
    note_en: str
    t_on_hetero_count: int = 0
    b_on_hetero_count: int = 0
    germline: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "uplift": round(float(self.uplift), 3),
            "flags": list(self.flags),
            "flag_labels_zh": list(self.flag_labels_zh),
            "flag_labels_en": list(self.flag_labels_en),
            "cdr_foreignness": round(float(self.cdr_foreignness), 4),
            "depleting": bool(self.depleting),
            "note_zh": self.note_zh,
            "note_en": self.note_en,
            "t_on_hetero_count": int(self.t_on_hetero_count),
            "b_on_hetero_count": int(self.b_on_hetero_count),
            "germline": self.germline,
        }


def _uplift_for_flags(flags: list[str], depleting: bool) -> float:
    if depleting:
        return 0.0
    n = len(flags)
    if n == 0:
        return 0.0
    if n == 1:
        only = flags[0]
        # Single MHC/B epitope density alone does not lift; need heterologous /
        # framework / aggregation corroboration.
        if only in {"t_cdr", "b_epi"}:
            return 0.0
        # Germline-overlap epitopes are stronger biology than lone t_cdr
        if only in {"t_on_hetero", "b_on_hetero"}:
            return 4.0
        if only in {"germline_hetero", "hetero_cdr", "fr_distant", "agg"}:
            return 3.0
        # lone charge: mild; keep at 0 to avoid noise
        return 0.0
    if n == 2:
        return 7.0
    return 12.0  # 3+


def evaluate_ab_risk_signals(
    *,
    sequence: str,
    regions,
    inn: str,
    cdr_epitope_density: float,
    non_germline_burden: float,
    b_coverage: float,
    agg_risk_hint: float,
    fr_mismatch: float,
    charge_density: float,
    cdr_foreignness: float | None = None,
    germline_report: GermlineMismatchReport | None = None,
    t_strong_starts: list[int] | None = None,
    b_epitope_spans: list[tuple[int, int]] | None = None,
    peptide_len: int = PEPTIDE_LENGTH,
) -> RiskSignalReport:
    """Evaluate multi-signal flags and capped uplift (0–15)."""
    foreign = (
        float(cdr_foreignness)
        if cdr_foreignness is not None
        else float(compute_cdr_foreignness(sequence, regions))
    )
    depleting = is_b_cell_depleting_cd19_20(inn)  # CD19/CD20 B-cell depleting only

    t_on_count = 0
    b_on_count = 0
    gl_dict: dict[str, Any] | None = None

    flags: list[str] = []
    if float(cdr_epitope_density) >= 35.0 or float(non_germline_burden) >= 18.0:
        flags.append("t_cdr")
    if float(b_coverage) >= 0.12:
        flags.append("b_epi")

    if germline_report is not None:
        gl_dict = germline_report.to_dict()
        # Prefer CDR-mismatch overlaps; fall back to any mismatch if CDR set empty
        mm_for_t = germline_report.cdr_mismatch_abs or germline_report.mismatch_abs
        if t_strong_starts:
            t_on_count = epitope_overlap_mismatches(
                list(t_strong_starts),
                int(peptide_len),
                mm_for_t,
                min_hit=3,
            )
        if b_epitope_spans:
            b_on_count = span_overlap_mismatches(
                list(b_epitope_spans),
                germline_report.cdr_mismatch_abs or germline_report.mismatch_abs,
                min_hit=3,
                one_based=True,
            )
        if (
            float(germline_report.fr_hetero_frac) >= _GERMLINE_FR_HETERO_FRAC
            or (
                float(germline_report.hetero_frac) >= _GERMLINE_HETERO_FRAC
                and float(fr_mismatch) >= 0.50
            )
            or (
                float(germline_report.cdr_hetero_frac) >= _GERMLINE_CDR_HETERO_FRAC
                and float(fr_mismatch) >= 0.50
            )
        ):
            flags.append("germline_hetero")
        # Human-like FR: count T∩hetero on FR mismatches only (CDR somatic hits are expected)
        if float(fr_mismatch) < 0.48 and t_strong_starts:
            t_on_count = epitope_overlap_mismatches(
                list(t_strong_starts),
                int(peptide_len),
                germline_report.fr_mismatch_abs or germline_report.mismatch_abs,
                min_hit=3,
            )
        if t_on_count >= _T_ON_HETERO_MIN:
            flags.append("t_on_hetero")
        # B∩hetero: require somewhat non-human FR so CDR-only somatic hits do not lift
        if b_on_count >= _B_ON_HETERO_MIN and float(fr_mismatch) >= 0.42:
            flags.append("b_on_hetero")
    else:
        # Fallback k-mer CDR foreignness proxy when germline alignment unavailable
        # CDR k-mer distances are naturally high vs FR; pool p50≈0.77 for human/humanized.
        # Threshold ~0.85 marks clearly atypical CDRs (not the literal 0.45 FR-scale cut).
        if foreign >= 0.85:
            flags.append("hetero_cdr")

    if float(agg_risk_hint) >= 45.0:
        flags.append("agg")
    if float(fr_mismatch) >= 0.55:
        flags.append("fr_distant")
    if abs(float(charge_density)) >= 0.025:
        flags.append("charge")

    uplift = min(_uplift_for_flags(flags, depleting), 15.0)
    labels_zh = [_FLAG_LABELS[f][0] for f in flags if f in _FLAG_LABELS]
    labels_en = [_FLAG_LABELS[f][1] for f in flags if f in _FLAG_LABELS]

    gl_bit_zh = ""
    gl_bit_en = ""
    if germline_report is not None:
        hits = ", ".join(f"{h.chain}:{h.v_gene}({h.identity:.0%})" for h in germline_report.hits)
        gl_bit_zh = (
            f"；种系最佳命中[{hits}]，FR错配分率={germline_report.fr_hetero_frac:.2f}，"
            f"t_on_hetero={t_on_count}"
        )
        gl_bit_en = (
            f"; germline best-hit[{hits}], FR mismatch frac={germline_report.fr_hetero_frac:.2f}, "
            f"t_on_hetero={t_on_count}"
        )

    if depleting:
        note_zh = (
            "CD19/CD20 B细胞耗竭靶点：仅标注风险信号，"
            f"不加表位上浮（flags={','.join(flags) or '无'}；"
            f"CDR异源性代理={foreign:.2f}{gl_bit_zh}）。"
        )
        note_en = (
            "CD19/CD20 B-cell depleting target: flags only, "
            f"no epitope uplift (flags={','.join(flags) or 'none'}; "
            f"CDR foreignness proxy={foreign:.2f}{gl_bit_en})."
        )
    elif uplift > 0:
        note_zh = (
            f"多信号风险：{len(flags)} 项触发（{', '.join(labels_zh)}），"
            f"表位上浮 +{uplift:.0f}（封顶 15）。"
            "t_on_hetero/b_on_hetero 单独即可 +4；单靠 MHC 密度通常不上浮。"
            f"{gl_bit_zh}"
        )
        note_en = (
            f"Multi-signal risk: {len(flags)} flags ({', '.join(labels_en)}), "
            f"epitope uplift +{uplift:.0f} (cap 15). "
            "t_on_hetero/b_on_hetero alone can give +4; lone MHC density usually does not."
            f"{gl_bit_en}"
        )
    elif flags:
        note_zh = (
            f"已记录风险信号（{', '.join(labels_zh)}），但未达多信号上浮门槛"
            "（单独 t_cdr/b_epi 或弱信号不加分）。"
            f"{gl_bit_zh}"
        )
        note_en = (
            f"Risk flags recorded ({', '.join(labels_en)}), but below multi-signal "
            "uplift threshold (lone t_cdr/b_epi or weak signals do not add)."
            f"{gl_bit_en}"
        )
    else:
        note_zh = "未触发多信号风险标志；标定表位保持 kNN 基线。" + gl_bit_zh
        note_en = "No multi-signal risk flags; calibrated epitope stays on kNN baseline." + gl_bit_en

    return RiskSignalReport(
        uplift=float(uplift),
        flags=flags,
        flag_labels_zh=labels_zh,
        flag_labels_en=labels_en,
        cdr_foreignness=foreign,
        depleting=depleting,
        note_zh=note_zh,
        note_en=note_en,
        t_on_hetero_count=int(t_on_count),
        b_on_hetero_count=int(b_on_count),
        germline=gl_dict,
    )
