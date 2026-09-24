"""
CDR foreignness proxy (NOT IMGT germline / ANARCI).

Windows-friendly heterologous signal: 4-mer Jaccard distance of CDR-only sequence
to a pool of low-ADA human / humanized marketed antibody CDRs. Mirrors the
framework_mismatch style in epitope_calibrate.py — a screening proxy only.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from .cdr import AntibodyRegionMap, annotate_antibody_regions

_KMER = 4

# B-cell depleting via CD19 / CD20 (target-mediated). NOT CD52 / alemtuzumab.
_B_CELL_DEPLETING_EXACT = frozenset(
    {
        "rituximab",
        "ofatumumab",
        "obinutuzumab",
        "ocrelizumab",
        "ublituximab",
        "tafasitamab",
        "inebilizumab",
        "blinatumomab",
        "epcoritamab",
        "glofitamab",
        "mosunetuzumab",
        "odronextamab",
        "loncastuximab",
    }
)
_B_CELL_DEPLETING_SUBSTR = tuple(sorted(_B_CELL_DEPLETING_EXACT))


def _cdr_only(sequence: str, regions: AntibodyRegionMap) -> str:
    return "".join(
        aa
        for i, aa in enumerate(sequence)
        if i < len(regions.cdr_mask) and regions.cdr_mask[i]
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


@lru_cache(maxsize=1)
def _low_ada_human_cdr_pool() -> tuple[str, ...]:
    """CDR-only sequences from low-ADA human / humanized marketed Abs (FR-based lineage)."""
    from .ada_db import ADAReferenceDB
    from .epitope_calibrate import framework_mismatch, infer_lineage

    try:
        db = ADAReferenceDB.from_csv()
    except FileNotFoundError:
        return tuple()

    cdrs: list[str] = []
    for _, r in db.df.iterrows():
        if float(r.get("ada_rate_pct") or 99) >= 10:
            continue
        heavy = str(r.get("heavy") or "")
        light = str(r.get("light") or "").rstrip("R")
        if len(heavy) < 80:
            continue
        seq = heavy + (light if len(light) >= 80 else "")
        reg = annotate_antibody_regions(seq)
        if not reg.is_antibody_like:
            continue
        fr = framework_mismatch(seq, reg)
        if infer_lineage(fr) not in {"human", "humanized"}:
            continue
        cdr = _cdr_only(seq, reg)
        if len(cdr) >= 20:
            cdrs.append(cdr)
    return tuple(cdrs)


def cdr_foreignness(sequence: str, regions: AntibodyRegionMap) -> float:
    """
    Heterologous / foreignness proxy in [0, 1].

    1 − mean top-5 4-mer Jaccard of CDR-only sequence vs low-ADA human CDR pool.
    Not IMGT germline identity — screening proxy only (no ANARCI required).
    """
    pool = _low_ada_human_cdr_pool()
    cdr = _cdr_only(sequence, regions)
    if len(cdr) < 20 or len(pool) < 5:
        return 0.5
    q = _kmers(cdr)
    sims = sorted((_jaccard(q, _kmers(p)) for p in pool), reverse=True)
    top = sims[: min(5, len(sims))]
    sim = float(np.mean(top)) if top else 0.0
    return float(np.clip(1.0 - sim, 0.0, 1.0))


def is_b_cell_depleting_cd19_20(inn: str = "", drug_class: str = "") -> bool:
    """
    Target-mediated B-cell depletion via CD19 / CD20.
    Used to dampen measured-ADA expectation; does NOT include CD52/alemtuzumab.
    """
    import re

    blob = f"{inn} {drug_class}".lower().replace("-", " ")
    if not blob.strip():
        return False
    if re.search(r"\bcd19\b|\bcd20\b|anti\s*cd19|anti\s*cd20", blob):
        return True
    token = blob.split()[0] if blob.split() else ""
    if token in _B_CELL_DEPLETING_EXACT:
        return True
    return any(s in token or token.startswith(s) for s in _B_CELL_DEPLETING_SUBSTR)


def is_depleting_mab(inn: str) -> bool:
    """Backward-compatible alias: CD19/CD20 B-cell depleting only."""
    return is_b_cell_depleting_cd19_20(inn, "")


def is_cd20_depleting_mab(inn: str) -> bool:
    """Alias for CD19/CD20 B-cell depleting (name kept for callers)."""
    return is_b_cell_depleting_cd19_20(inn, "")


def infer_ab_use_context(inn: str = "", drug_class: str = "") -> list[str]:
    """
    Product / use-context tags (NOT clinical ADA rates).
    b_cell_depleting_cd19_20: target-mediated B-cell depletion → dampen score.
    Chimeric scaffolds stay high unless this (or similar format) tag fires.
    """
    import re

    blob = f"{inn} {drug_class}".lower()
    tags: list[str] = []
    if is_b_cell_depleting_cd19_20(inn, drug_class):
        tags.append("b_cell_depleting_cd19_20")
    if re.search(
        r"bispecific|multispecific|adc|antibody-drug|conjugate|immunotoxin|"
        r"vedotin|ozogamicin|tesirine|pasudotox",
        blob,
    ) and "b_cell_depleting_cd19_20" not in tags:
        # Non-CD19/20 engagers/ADCs still get format dampening
        tags.append("engager_adc")
    # Acute antidote / single-use fragments only — do NOT tag chronic scFv (e.g. brolucizumab).
    # Format exposure for scFv/Fab is handled by ab_format.detect_ab_format (raises foreignness).
    if re.search(r"antidote|idarucizumab|digoxin|andexanet", blob):
        tags.append("fragment_acute")
    if re.search(
        r"cetuxi|panitumu|necitumu|dinutuxi|naxitam|isatuxi|margetuxi|"
        r"siltuxi|tositumo|ibritumo|obilt|brentux",
        blob,
    ) and "b_cell_depleting_cd19_20" not in tags:
        tags.append("oncology_mab")
    return tags
