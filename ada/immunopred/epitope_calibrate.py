"""
Predictive antibody epitope score for early screening.

Clinical ADA is an outcome reference — not used to correct the prediction.
Pipeline:
  1) Foreignness (germline hetero + FR distance + CDR proxy + format exposure)
  2) T/B epitope burden on CDR / non-germline
  3) Lineage tolerance (human > humanized > chimeric)
  4) Optional multi-signal uplift (hetero∩T/B flags; not clinical kNN)
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import numpy as np

from .cdr import AntibodyRegionMap, annotate_antibody_regions, region_epitope_burden
from .t_cell import PEPTIDE_LENGTH, predict_t_cell_epitopes
from .sequence import from_string


_LINEAGE_CODE = {"human": 0.0, "humanized": 0.5, "chimeric": 1.0, "unknown": 0.45}
_KMER = 4


@dataclass
class AbEpitopeBreakdown:
    score: float
    lineage_inferred: str
    fr_mismatch: float
    ng_burden: float
    charge_term: float
    humanness_term: float
    burden_term: float
    knn_clinical: float
    knn_n: int
    note_zh: str
    note_en: str
    cdr_foreignness: float = 0.0
    risk_uplift: float = 0.0
    risk_flags: list[str] | None = None
    risk_flag_labels_zh: list[str] | None = None
    risk_flag_labels_en: list[str] | None = None
    risk_note_zh: str = ""
    risk_note_en: str = ""
    base_score: float | None = None
    germline: dict[str, Any] | None = None
    t_on_hetero_count: int = 0
    b_on_hetero_count: int = 0
    ab_format: dict[str, Any] | None = None
    format_exposure: float = 0.0
    pre_context_score: float = 0.0
    context_tags: list[str] | None = None
    context_note_zh: str = ""
    context_note_en: str = ""
    context_mult: float = 1.0
    context_cap: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "score": round(self.score, 3),
            "lineage_inferred": self.lineage_inferred,
            "fr_mismatch": round(self.fr_mismatch, 4),
            "ng_burden": round(self.ng_burden, 4),
            "charge_term": round(self.charge_term, 3),
            "humanness_term": round(self.humanness_term, 3),
            "burden_term": round(self.burden_term, 3),
            "knn_clinical": round(self.knn_clinical, 3),
            "knn_n": self.knn_n,
            "note_zh": self.note_zh,
            "note_en": self.note_en,
            "cdr_foreignness": round(float(self.cdr_foreignness), 4),
            "risk_uplift": round(float(self.risk_uplift), 3),
            "risk_flags": list(self.risk_flags or []),
            "flag_labels_zh": list(self.risk_flag_labels_zh or []),
            "flag_labels_en": list(self.risk_flag_labels_en or []),
            "risk_note_zh": self.risk_note_zh,
            "risk_note_en": self.risk_note_en,
            "base_score": round(float(self.base_score if self.base_score is not None else self.score), 3),
            "germline": self.germline,
            "t_on_hetero_count": int(self.t_on_hetero_count),
            "b_on_hetero_count": int(self.b_on_hetero_count),
            "ab_format": self.ab_format,
            "format_exposure": round(float(self.format_exposure), 4),
            "pre_context_score": round(float(self.pre_context_score), 3),
            "context_tags": list(self.context_tags or []),
            "context_note_zh": self.context_note_zh,
            "context_note_en": self.context_note_en,
            "context_mult": round(float(self.context_mult), 3),
            "context_cap": self.context_cap,
        }


def _fr_only(sequence: str, regions: AntibodyRegionMap) -> str:
    return "".join(
        aa
        for i, aa in enumerate(sequence)
        if i < len(regions.cdr_mask) and not regions.cdr_mask[i]
    )


def _kmers(seq: str, k: int = _KMER) -> set[str]:
    s = "".join(c for c in (seq or "").upper() if c.isalpha())
    if len(s) < k:
        return {s} if s else set()
    return {s[i : i + k] for i in range(len(s) - k + 1)}


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def _inn_stem(inn: str) -> str:
    """Strip conjugate / combo suffixes so related products leave out together."""
    t = (inn or "").strip().lower().replace("-", " ")
    # keep first token; drop common ADC/payload words from multiword names
    parts = t.split()
    if not parts:
        return ""
    stem = parts[0]
    for suffix in (
        "pasudotox",
        "vedotin",
        "emtansine",
        "ozogamicin",
        "tesirine",
        "deruxtecan",
    ):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)].rstrip()
    return stem


@lru_cache(maxsize=1)
def _human_fr_pool() -> tuple[str, ...]:
    """Framework-only sequences from low-ADA human / humanized references."""
    from .ada_db import ADAReferenceDB
    from .model import _inn_antibody_lineage

    try:
        db = ADAReferenceDB.from_csv()
    except FileNotFoundError:
        return tuple()

    frs: list[str] = []
    for _, r in db.df.iterrows():
        if float(r.get("ada_rate_pct") or 99) >= 10:
            continue
        lin = _inn_antibody_lineage(str(r.get("inn") or ""))
        if lin not in {"human", "humanized"}:
            continue
        heavy = str(r.get("heavy") or "")
        light = str(r.get("light") or "").rstrip("R")
        if len(heavy) < 80:
            continue
        seq = heavy + (light if len(light) >= 80 else "")
        reg = annotate_antibody_regions(seq)
        if not reg.is_antibody_like:
            continue
        fr = _fr_only(seq, reg)
        if len(fr) >= 40:
            frs.append(fr)
    return tuple(frs)


def framework_mismatch(sequence: str, regions: AntibodyRegionMap) -> float:
    """
    Distance to low-ADA human frameworks in [0, 1].

    Uses 4-mer Jaccard against the nearest human-like FRs (not columnar consensus),
    so CDR-length differences do not inflate mismatch.
    """
    pool = _human_fr_pool()
    fr = _fr_only(sequence, regions)
    if len(fr) < 40 or len(pool) < 5:
        return 0.5
    q = _kmers(fr)
    sims = sorted((_jaccard(q, _kmers(p)) for p in pool), reverse=True)
    top = sims[: min(5, len(sims))]
    sim = float(np.mean(top)) if top else 0.0
    # Map similarity→distance; typical human Fv ~0.55–0.85 k-mer Jaccard
    return float(np.clip(1.0 - sim, 0.0, 1.0))


def infer_lineage(
    fr_mismatch: float,
    inn: str = "",
    germline_report=None,
) -> str:
    """
    Structure-only lineage — ignore WHO INN suffix (-umab/-zumab/-ximab).

    Prefer human V-gene germline best-hit (standard humanness proxy):
      chimeric  ≈ murine V on human C → high FR hetero / low V identity
      humanized ≈ CDR graft on human FR → moderate FR hetero
      human     ≈ human repertoire / transgenic → low FR hetero

    FR k-mer distance to a marketed pool is fallback only (can false-call
    distant-but-human V genes as chimeric).
    """
    del inn
    fr_h = None
    mean_id = None
    if germline_report is not None:
        gl = germline_report.to_dict() if hasattr(germline_report, "to_dict") else germline_report
        if isinstance(gl, dict):
            fr_h = float(gl.get("fr_hetero_frac") or 0.0) or None
            hits = gl.get("hits") or []
            idents = [float(h.get("identity") or 0.0) for h in hits if float(h.get("identity") or 0) > 0]
            if idents:
                mean_id = float(sum(idents) / len(idents))
        else:
            fr_h = float(getattr(germline_report, "fr_hetero_frac", 0) or 0) or None
            hits = getattr(germline_report, "hits", None) or []
            idents = [float(h.identity) for h in hits if float(getattr(h, "identity", 0) or 0) > 0]
            if idents:
                mean_id = float(sum(idents) / len(idents))

    if fr_h is not None and fr_h > 0:
        # Thresholds from IMGT-style V identity / FR mismatch (Infliximab-class vs Evolocumab-class)
        if fr_h >= 0.28 or (mean_id is not None and mean_id < 0.68):
            return "chimeric"
        if fr_h >= 0.14 or (mean_id is not None and mean_id < 0.82):
            return "humanized"
        return "human"

    # Fallback: marketed-pool FR k-mer distance
    if fr_mismatch >= 0.55:
        return "chimeric"
    if fr_mismatch >= 0.38:
        return "humanized"
    return "human"


def scaffold_foreignness(
    fr_mismatch: float,
    germline_report=None,
) -> float:
    """
    Scaffold non-humanness in [0, 1].

    Prefer germline FR heterologous fraction (best-hit human V); fall back to
    FR k-mer distance. Maps fr_hetero onto a comparable 0–1 scale where
    chimeric FR hetero ~0.33 → ~0.70 distance-equivalent for scoring weights.
    """
    fr_h = 0.0
    if germline_report is not None:
        gl = germline_report.to_dict() if hasattr(germline_report, "to_dict") else germline_report
        if isinstance(gl, dict):
            fr_h = float(gl.get("fr_hetero_frac") or 0.0)
        else:
            fr_h = float(getattr(germline_report, "fr_hetero_frac", 0) or 0.0)
    if fr_h > 0:
        # fr_hetero 0.08 (Ada) → ~0.20; 0.15 (Evo) → ~0.35; 0.33 (Inflix) → ~0.75
        return float(np.clip(fr_h / 0.44, 0.0, 1.0))
    return float(np.clip(fr_mismatch, 0.0, 1.0))


def _feature_vector(
    fr_mism: float,
    ng: float,
    charge: float,
    lineage: str,
    agg: float = 0.0,
) -> np.ndarray:
    # physchem.agg_risk_hint is typically ~20–50 for Fv; scale into ~[0, 1.2]
    agg_n = min(max(float(agg), 0.0) / 45.0, 1.2)
    return np.array(
        [
            fr_mism,
            ng / 25.0,
            min(abs(charge) * 10.0, 1.5),
            _LINEAGE_CODE.get(lineage, 0.45),
            agg_n,
        ],
        dtype=float,
    )


@lru_cache(maxsize=1)
def _reference_ab_table() -> tuple[np.ndarray, np.ndarray, tuple[str, ...], tuple[str, ...]]:
    """Return (feature_matrix, clinical_ada, inns, stems) for antibody reference rows."""
    from .ada_db import ADAReferenceDB
    from .physchem import profile_physchem

    db = ADAReferenceDB.from_csv()
    feats: list[np.ndarray] = []
    y: list[float] = []
    inns: list[str] = []
    stems: list[str] = []
    for _, r in db.df.iterrows():
        heavy = str(r.get("heavy") or "")
        if len(heavy) < 80:
            continue
        inn = str(r.get("inn") or "")
        light = str(r.get("light") or "").rstrip("R")
        seq = heavy + (light if len(light) >= 80 else "")
        reg = annotate_antibody_regions(seq)
        if not reg.is_antibody_like:
            continue
        prot = from_string(seq, inn)
        t = predict_t_cell_epitopes(prot)
        burden = region_epitope_burden(list(t.strong_starts or []), PEPTIDE_LENGTH, reg.cdr_mask)
        ng = float(burden["non_germline_burden"])
        fr_mism = framework_mismatch(seq, reg)
        lineage = infer_lineage(fr_mism, inn=inn)
        phys = profile_physchem(seq)
        feats.append(
            _feature_vector(
                fr_mism,
                ng,
                phys.charge_density,
                lineage,
                float(phys.agg_risk_hint or 0.0),
            )
        )
        y.append(float(r["ada_rate_pct"]))
        inns.append(inn.lower())
        stems.append(_inn_stem(inn))
    if not feats:
        return np.zeros((0, 5)), np.zeros(0), tuple(), tuple()
    return np.vstack(feats), np.asarray(y, dtype=float), tuple(inns), tuple(stems)


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cdf = np.cumsum(w)
    cutoff = 0.5 * float(w.sum())
    return float(v[int(np.searchsorted(cdf, cutoff, side="left"))])


@lru_cache(maxsize=1)
def _fr_stratum_prior_table() -> tuple[float, float, float]:
    """Mean clinical ADA% for low / mid / high FR-distance marketed Abs."""
    X, y, _, _ = _reference_ab_table()
    if len(y) == 0:
        return 5.0, 8.0, 16.0
    fr = X[:, 0]
    low = float(np.mean(y[fr < 0.35])) if np.any(fr < 0.35) else 5.0
    mid = float(np.mean(y[(fr >= 0.35) & (fr < 0.55)])) if np.any((fr >= 0.35) & (fr < 0.55)) else 8.0
    high = float(np.mean(y[fr >= 0.55])) if np.any(fr >= 0.55) else 16.0
    return low, mid, high


def _fr_stratum_prior(fr_mism: float) -> float:
    low, mid, high = _fr_stratum_prior_table()
    if fr_mism < 0.35:
        return low
    if fr_mism < 0.55:
        return mid
    return high


def _knn_clinical(vec: np.ndarray, inn: str = "", k: int = 11) -> tuple[float, int]:
    X, y, inns, stems = _reference_ab_table()
    if len(y) == 0:
        return 12.0, 0

    stem = _inn_stem(inn) if inn else ""
    mask = np.ones(len(y), dtype=bool)
    if stem:
        mask = np.array([s != stem for s in stems], dtype=bool)
        if int(mask.sum()) < max(8, k):
            key = inn.strip().lower()
            mask = np.array([i != key for i in inns], dtype=bool)
        if int(mask.sum()) < max(5, k):
            mask = np.ones(len(y), dtype=bool)

    X2, y2 = X[mask], y[mask]
    mu = X2.mean(axis=0)
    sd = X2.std(axis=0)
    sd = np.where(sd < 1e-6, 1.0, sd)
    # Emphasize FR distance & non-germline burden over charge/agg (noisier)
    feat_w = np.array([1.8, 1.3, 0.55, 0.9, 0.7], dtype=float)
    Z = (X2 - mu) / sd * feat_w
    q = (vec - mu) / sd * feat_w
    d = np.sqrt(((Z - q) ** 2).sum(axis=1))
    kk = min(k, len(y2))
    idx = np.argpartition(d, kk - 1)[:kk]
    w = 1.0 / (d[idx] + 0.06)

    # Robust center; mild high-tail for non-humanlike FR (humanized/chimeric risk)
    knn_med = _weighted_median(y2[idx], w)
    knn_mean = float(np.average(y2[idx], weights=w))
    knn_core = 0.70 * knn_med + 0.30 * knn_mean
    fr = float(vec[0])
    order = np.argsort(y2[idx])
    top = y2[idx][order[max(0, (2 * len(order)) // 3) :]]
    knn_hi = float(np.average(top)) if len(top) else knn_core
    hi_w = float(np.clip((fr - 0.48) / 0.32, 0.0, 0.40)) if fr >= 0.48 else 0.0
    knn = (1.0 - hi_w) * knn_core + hi_w * knn_hi
    return knn, int(kk)


_LINEAGE_TOLERANCE = {
    "human": 0.55,
    "humanized": 0.75,
    "chimeric": 1.15,
    "unknown": 0.70,
}


def calibrated_ab_epitope_score(
    sequence: str,
    regions: AntibodyRegionMap,
    non_germline_burden: float,
    charge_density: float,
    agg_risk_hint: float,
    inn: str = "",
    drug_class: str = "",
    b_coverage: float = 0.0,
    cdr_epitope_density: float | None = None,
    germline_report=None,
    t_strong_starts: list[int] | None = None,
    b_epitope_spans: list[tuple[int, int]] | None = None,
) -> AbEpitopeBreakdown:
    """Predictive Ab score: FR/germline foreignness × T/B burden × lineage weight.

    Clinical ADA neighbors are not used to correct the score.
    CDR k-mer foreignness is only a small excess term — CDRs are meant to vary.
    """
    from .ab_format import detect_ab_format
    from .germline_proxy import (
        cdr_foreignness,
        infer_ab_use_context,
        is_b_cell_depleting_cd19_20,
    )
    from .risk_signals import evaluate_ab_risk_signals

    fr_mism = framework_mismatch(sequence, regions)
    ng = float(non_germline_burden)
    agg = float(agg_risk_hint or 0.0)
    foreign = cdr_foreignness(sequence, regions)
    cdr_dens = float(cdr_epitope_density) if cdr_epitope_density is not None else 0.0
    fmt = detect_ab_format(sequence, regions, inn=inn, drug_class=drug_class)
    fmt_exp = float(fmt.exposure)
    ctx = infer_ab_use_context(inn, drug_class)
    # Ensure CD19/CD20 depleting tag even if only inn stem matched
    if is_b_cell_depleting_cd19_20(inn, drug_class) and "b_cell_depleting_cd19_20" not in ctx:
        ctx.append("b_cell_depleting_cd19_20")

    gl = germline_report
    if gl is not None and hasattr(gl, "to_dict"):
        gl_d = gl.to_dict()
    elif isinstance(gl, dict):
        gl_d = gl
    else:
        gl_d = {}
    lineage = infer_lineage(fr_mism, inn=inn, germline_report=germline_report or gl_d)
    scaffold = scaffold_foreignness(fr_mism, germline_report=germline_report or gl_d)
    hetero_frac = float(gl_d.get("hetero_frac") or 0.0)
    fr_hetero = float(gl_d.get("fr_hetero_frac") or 0.0)
    cdr_hetero = float(gl_d.get("cdr_hetero_frac") or 0.0)
    if hetero_frac <= 0:
        hetero_frac = float(scaffold)
    # Germline FR hetero for scaffold; overall hetero as mild mix
    hetero_for_f = (
        0.75 * fr_hetero + 0.25 * hetero_frac if fr_hetero > 0 else hetero_frac
    )

    # 1) Foreignness 0–100 — scaffold from germline FR (not marketed FR k-mer)
    cdr_excess = float(np.clip((float(foreign) - 0.78) / 0.22, 0.0, 1.0))
    cdr_gl = float(np.clip(cdr_hetero, 0.0, 1.0))
    foreignness = float(
        np.clip(
            100.0
            * (
                0.35 * float(scaffold)
                + 0.20 * float(hetero_for_f)
                + 0.20 * cdr_gl
                + 0.05 * cdr_excess
                + 0.20 * fmt_exp
            ),
            0.0,
            100.0,
        )
    )
    relevance = 0.20 + 0.80 * (foreignness / 100.0)

    # 2) Epitope burden — CDR-local T/B; scaffold term from germline FR
    t_term = float(np.clip(cdr_dens / 40.0 * 100.0, 0, 100))
    ng_term = float(np.clip(ng / 35.0 * 100.0, 0, 100))
    b_term = float(np.clip(float(b_coverage or 0.0) * 100.0, 0, 100))
    fr_term = float(np.clip(max(0.0, scaffold - 0.25) / 0.55 * 100.0, 0, 100))
    epitope = 0.32 * t_term + 0.33 * ng_term + 0.20 * b_term + 0.15 * fr_term

    humanness_term = 24.0 * max(0.0, scaffold - 0.28)
    burden_term = 0.40 * (ng - 14.0)
    charge_term = 40.0 * abs(float(charge_density))

    # 3) Combine: foreignness floor + epitope gated by foreignness, then lineage weight
    lw = float(_LINEAGE_TOLERANCE.get(lineage, 0.70))
    base_score = float(
        np.clip(
            (0.45 * foreignness * lw + 0.55 * epitope * relevance),
            1.0,
            80.0,
        )
    )

    # Soft ceiling for highly human scaffolds (germline FR / lineage)
    if fmt_exp < 0.50:
        if lineage == "human" and scaffold < 0.40:
            base_score = min(base_score, 22.0)
        elif lineage == "humanized" and scaffold < 0.45:
            base_score = min(base_score, 16.0)
        elif lineage == "humanized" and scaffold < 0.55:
            base_score = min(base_score, 19.0)
    if lineage == "chimeric":
        base_score = min(base_score, 58.0)

    risk = evaluate_ab_risk_signals(
        sequence=sequence,
        regions=regions,
        inn=inn,
        cdr_epitope_density=cdr_dens,
        non_germline_burden=ng,
        b_coverage=float(b_coverage or 0.0),
        agg_risk_hint=agg,
        fr_mismatch=scaffold,  # risk flags use germline scaffold distance
        charge_density=float(charge_density),
        cdr_foreignness=foreign,
        germline_report=germline_report,
        t_strong_starts=t_strong_starts,
        b_epitope_spans=b_epitope_spans,
    )
    # Uplift only when framework is meaningfully non-human
    uplift = float(risk.uplift)
    if lineage == "human":
        uplift *= float(np.clip((scaffold - 0.28) / 0.35, 0.0, 1.0))
    elif lineage == "humanized":
        uplift *= float(np.clip((scaffold - 0.32) / 0.40, 0.0, 1.0))
    if lineage == "humanized" and cdr_hetero >= 0.58 and scaffold >= 0.45:
        uplift = max(
            uplift,
            10.0 + 5.0 * float(np.clip((cdr_hetero - 0.58) / 0.25, 0, 1)),
        )
    # Format-driven uplift: scFv / VHH expose non-natural surfaces → pre-existing ADA
    if fmt.format in {"scfv", "fv", "vhh"} and fmt_exp >= 0.55:
        uplift = max(
            uplift,
            14.0 + 10.0 * float(np.clip((fmt_exp - 0.55) / 0.30, 0, 1)),
        )

    score = float(np.clip(base_score + uplift, 1.0, 85.0))
    pre_context_score = score

    # Use-context tolerance (NOT clinical ADA). Chimeric stays high unless tagged.
    ctx_note = ""
    ctx_note_en = ""
    ctx_mult = 1.0
    ctx_cap: float | None = None
    if "b_cell_depleting_cd19_20" in ctx:
        # Target-mediated B-cell depletion (CD19/CD20): suppress measured-ADA expectation
        # Clinical rates for this class are typically ~0–5%.
        score = float(np.clip(score * 0.12, 1.0, 6.0))
        ctx_mult, ctx_cap = 0.12, 6.0
        ctx_note = "CD19/CD20 B细胞耗竭靶点：测得ADA预期×0.12、封顶6；"
        ctx_note_en = "CD19/CD20 B-cell depleting target: measured-ADA expectation ×0.12, cap 6; "
    elif "engager_adc" in ctx:
        score = float(np.clip(score * 0.48, 1.0, 28.0))
        ctx_mult, ctx_cap = 0.48, 28.0
        ctx_note = "双抗/ADC/衔接器情境耐受↑；"
        ctx_note_en = "BsAb/ADC/engager context tolerance ↑; "
    elif "oncology_mab" in ctx:
        score = float(np.clip(score * 0.40, 1.0, 24.0))
        ctx_mult, ctx_cap = 0.40, 24.0
        ctx_note = "肿瘤向全抗情境耐受↑；"
        ctx_note_en = "Oncology mAb context tolerance ↑; "
    elif "fragment_acute" in ctx:
        score = float(np.clip(score * 0.40, 1.0, 18.0))
        ctx_mult, ctx_cap = 0.40, 18.0
        ctx_note = "急救/解毒片段单次情境耐受↑；"
        ctx_note_en = "Rescue/antidote fragment single-dose context tolerance ↑; "
    # else: chimeric / foreign scaffolds keep full structural score

    if (
        lineage == "humanized"
        and cdr_hetero >= 0.62
        and scaffold >= 0.45
        and "b_cell_depleting_cd19_20" not in ctx
    ):
        score = max(score, 42.0)
        ctx_note += "高CDR种系异源人源化下限；"
        ctx_note_en += "High CDR germline hetero humanized floor; "
    # scFv floor only when not B-cell depleting (depletion overrides measured-ADA expectation)
    if (
        fmt.format in {"scfv", "fv"}
        and fmt_exp >= 0.60
        and "b_cell_depleting_cd19_20" not in ctx
    ):
        score = max(score, 52.0)
        ctx_note += f"格式暴露面抬升（{fmt.format}）；"
        ctx_note_en += f"Format exposure uplift ({fmt.format}); "

    note_zh = (
        f"预测链路：异源度={foreignness:.0f}（骨架异源={scaffold:.2f}，"
        f"FR种系异源={fr_hetero:.2f}，CDR种系异源={cdr_hetero:.2f}，"
        f"FR k-mer距离={fr_mism:.2f}，格式暴露={fmt_exp:.2f}/{fmt.format}；"
        f"lineage={lineage}按种系FR）→ "
        f"表位相关度={relevance:.2f}；"
        f"CDR加权表位负担={epitope:.1f}；人源化权重={lineage}×{lw:.2f} → "
        f"基线={base_score:.1f}，FR/格式上浮=+{uplift:.1f}。"
        f"{ctx_note}"
        f"情境标签={{{','.join(ctx) or '无'}}}。"
        "临床 ADA 为结果对照，不参与矫正。"
    )
    note_en = (
        f"Predictive path: foreignness={foreignness:.0f} (scaffold={scaffold:.2f}, "
        f"FR germline hetero={fr_hetero:.2f}, CDR germline hetero={cdr_hetero:.2f}, "
        f"FR k-mer distance={fr_mism:.2f}, format exposure={fmt_exp:.2f}/{fmt.format}; "
        f"lineage={lineage} from germline FR) → "
        f"epitope relevance={relevance:.2f}; "
        f"CDR-weighted epitope burden={epitope:.1f}; lineage weight={lineage}×{lw:.2f} → "
        f"base={base_score:.1f}, FR/format uplift=+{uplift:.1f}. "
        f"Context tags={{{','.join(ctx) or 'none'}}}. "
        "Clinical ADA is an outcome reference, not a correction."
    )
    if risk.note_zh:
        note_zh = f"{note_zh} {risk.note_zh}"
    if risk.note_en:
        note_en = f"{note_en} {risk.note_en}"
    note_zh = f"{note_zh} {fmt.method_note_zh}"
    note_en = f"{note_en} {fmt.method_note_en}"

    return AbEpitopeBreakdown(
        score=score,
        lineage_inferred=lineage,
        fr_mismatch=fr_mism,
        ng_burden=ng,
        charge_term=charge_term,
        humanness_term=humanness_term,
        burden_term=burden_term,
        knn_clinical=foreignness,
        knn_n=0,
        note_zh=note_zh,
        note_en=note_en,
        cdr_foreignness=foreign,
        risk_uplift=uplift,
        risk_flags=list(risk.flags or []),
        risk_flag_labels_zh=list(getattr(risk, "flag_labels_zh", None) or []),
        risk_flag_labels_en=list(getattr(risk, "flag_labels_en", None) or []),
        risk_note_zh=risk.note_zh,
        risk_note_en=risk.note_en,
        base_score=base_score,
        germline=gl_d or None,
        t_on_hetero_count=int(getattr(risk, "t_on_hetero_count", 0) or 0),
        b_on_hetero_count=int(getattr(risk, "b_on_hetero_count", 0) or 0),
        ab_format=fmt.to_dict(),
        format_exposure=fmt_exp,
        pre_context_score=pre_context_score,
        context_tags=list(ctx or []),
        context_note_zh=ctx_note,
        context_note_en=ctx_note_en,
        context_mult=ctx_mult,
        context_cap=ctx_cap,
    )
