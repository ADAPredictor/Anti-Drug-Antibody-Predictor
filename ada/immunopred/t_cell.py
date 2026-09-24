"""
MHC class II (CD4 T-cell) epitope risk heuristics for protein therapeutics.

Clinical ADA responses are typically T-cell dependent. This module scores
15-mer peptides for binding propensity to common HLA-DR alleles using
transparent pocket-preference matrices (research-grade, not a NetMHCIIpan
replacement).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from .sequence import ProteinSequence

# Default population-relevant HLA-DRB1 alleles
DEFAULT_ALLELES = [
    "DRB1*01:01",
    "DRB1*03:01",
    "DRB1*04:01",
    "DRB1*07:01",
    "DRB1*11:01",
    "DRB1*15:01",
]

PEPTIDE_LENGTH = 15
CORE_LENGTH = 9


@dataclass
class TCellEpitope:
    start: int  # 1-based
    peptide: str
    core: str
    allele: str
    score: float
    percentile: float
    strong: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _builtin_matrices() -> dict[str, Any]:
    """Fallback pocket preferences if data file is missing."""
    hydrophobic = list("FILMWYVA")
    aromatic = list("FYW")
    positive = list("KRH")
    negative = list("DE")
    polar = list("STNQ")
    small = list("AGST")

    def base_pssm() -> dict[str, list[float]]:
        aa = "ACDEFGHIKLMNPQRSTVWY"
        # positions P1..P9
        mat = {a: [0.0] * 9 for a in aa}
        for a in hydrophobic:
            mat[a][0] += 1.2  # P1
            mat[a][3] += 0.4
            mat[a][8] += 0.6
        for a in aromatic:
            mat[a][0] += 0.8
        for a in positive:
            mat[a][0] -= 1.0
            mat[a][3] -= 0.3
        for a in negative:
            mat[a][3] -= 0.2
            mat[a][8] -= 0.5
        for a in polar:
            mat[a][5] += 0.3
        for a in small:
            mat[a][5] += 0.4
        return mat

    matrices = {}
    for allele in DEFAULT_ALLELES:
        pssm = base_pssm()
        # Mild allele-specific tweaks
        if "03:01" in allele:
            for a in negative:
                pssm[a][3] += 0.8
            for a in positive:
                pssm[a][8] += 0.5
        elif "04:01" in allele:
            for a in list("DE"):
                pssm[a][3] -= 0.6
            for a in hydrophobic:
                pssm[a][3] += 0.5
        elif "07:01" in allele:
            for a in aromatic + list("ILV"):
                pssm[a][0] += 0.5
                pssm[a][8] += 0.4
        elif "11:01" in allele:
            for a in list("R"):
                pssm[a][5] += 0.6
            for a in list("DE"):
                pssm[a][5] -= 0.4
        elif "15:01" in allele:
            for a in list("ILVFM"):
                pssm[a][0] += 0.4
                pssm[a][8] += 0.5
        matrices[allele] = {"pssm": pssm, "threshold_strong": 1.8}
    return {"alleles": matrices, "core_offset": 3}


@lru_cache(maxsize=1)
def load_hla_matrices(path: str | None = None) -> dict[str, Any]:
    if path:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return data
    # Prefer packaged / repo data file
    candidates = [
        Path(__file__).resolve().parents[2] / "data" / "hla_dr_matrices.json",
        Path.cwd() / "data" / "hla_dr_matrices.json",
    ]
    for cand in candidates:
        if cand.exists():
            return json.loads(cand.read_text(encoding="utf-8"))
    return _builtin_matrices()


def _score_core(core: str, pssm: dict[str, list[float]]) -> float:
    """Sum of pocket PSSM weights over the 9-mer core (not per-residue mean).

    ``threshold_strong`` (~1.8) and downstream ``mean_top / 2.5`` scaling assume
    this sum scale. Averaging by CORE_LENGTH shrinks scores to ~0.3–0.5 so that
    the absolute floor ``threshold_strong * 0.5`` is never reached.
    """
    if len(core) != CORE_LENGTH:
        raise ValueError("core must be 9 residues")
    total = 0.0
    for i, aa in enumerate(core):
        total += pssm.get(aa, [0.0] * 9)[i]
    return total


def _best_core(peptide: str, pssm: dict[str, list[float]], core_offset: int = 3) -> tuple[str, float, int]:
    """
    Evaluate possible 9-mer cores within a 15-mer.
    Default offset=3 mirrors common NetMHCII framing (positions 4-12),
    but we also scan nearby frames and keep the best.
    """
    best_score = -1e9
    best_core = peptide[core_offset : core_offset + CORE_LENGTH]
    best_pos = core_offset
    for start in range(0, len(peptide) - CORE_LENGTH + 1):
        core = peptide[start : start + CORE_LENGTH]
        score = _score_core(core, pssm)
        # Prefer the canonical offset when scores are close
        adjusted = score + (0.05 if start == core_offset else 0.0)
        if adjusted > best_score:
            best_score = adjusted
            best_core = core
            best_pos = start
    return best_core, best_score, best_pos


@dataclass
class TCellResult:
    epitopes: list[TCellEpitope]
    strong_count: int
    epitope_density: float  # strong binders / 100 aa
    mean_top_score: float
    allele_summary: dict[str, dict[str, float]]
    strong_starts: list[int] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "strong_count": self.strong_count,
            "epitope_density": self.epitope_density,
            "mean_top_score": self.mean_top_score,
            "allele_summary": self.allele_summary,
            "epitopes": [e.to_dict() for e in self.epitopes],
            "strong_starts": list(self.strong_starts or []),
        }


def predict_t_cell_epitopes(
    protein: ProteinSequence,
    alleles: list[str] | None = None,
    top_n: int = 30,
    strong_percentile: float = 10.0,
    matrices_path: str | None = None,
) -> TCellResult:
    """
    Predict MHC-II binding propensity along a protein sequence.

    Parameters
    ----------
    strong_percentile:
        Peptides at or above this percentile of scores (per allele) are
        marked as strong binder candidates.
    """
    alleles = alleles or DEFAULT_ALLELES
    data = load_hla_matrices(matrices_path)
    allele_mats = data["alleles"]
    core_offset = int(data.get("core_offset", 3))

    peptides = protein.peptides(PEPTIDE_LENGTH, step=1)
    all_hits: list[TCellEpitope] = []
    allele_summary: dict[str, dict[str, float]] = {}

    for allele in alleles:
        if allele not in allele_mats:
            raise KeyError(f"No matrix for allele {allele}")
        pssm = allele_mats[allele]["pssm"]
        fixed_thr = float(allele_mats[allele].get("threshold_strong", 1.8))

        scored: list[tuple[int, str, str, float]] = []
        for start, pep in peptides:
            core, score, _ = _best_core(pep, pssm, core_offset=core_offset)
            scored.append((start, pep, core, score))

        if not scored:
            allele_summary[allele] = {
                "n_peptides": 0,
                "strong": 0,
                "mean_score": 0.0,
                "max_score": 0.0,
            }
            continue

        scores = np.array([s[3] for s in scored], dtype=float)
        # Percentile rank within this protein/allele (higher = stronger)
        order = scores.argsort()
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.linspace(0, 100, len(scores))

        strong = 0
        for (start, pep, core, score), pct in zip(scored, ranks):
            is_strong = pct >= (100.0 - strong_percentile) and score >= fixed_thr * 0.5
            if is_strong:
                strong += 1
            all_hits.append(
                TCellEpitope(
                    start=start,
                    peptide=pep,
                    core=core,
                    allele=allele,
                    score=float(score),
                    percentile=float(pct),
                    strong=bool(is_strong),
                )
            )

        allele_summary[allele] = {
            "n_peptides": float(len(scored)),
            "strong": float(strong),
            "mean_score": float(scores.mean()),
            "max_score": float(scores.max()),
        }

    strong_epitopes = [e for e in all_hits if e.strong]
    strong_epitopes.sort(key=lambda e: e.score, reverse=True)
    # Deduplicate overlapping reports for density: unique start positions across alleles
    unique_strong_starts = {e.start for e in strong_epitopes}
    density = (len(unique_strong_starts) / max(protein.length, 1)) * 100.0

    top = strong_epitopes[:top_n]
    mean_top = float(np.mean([e.score for e in top])) if top else 0.0

    return TCellResult(
        epitopes=top,
        strong_count=len(unique_strong_starts),
        epitope_density=density,
        mean_top_score=mean_top,
        allele_summary=allele_summary,
        strong_starts=sorted(unique_strong_starts),
    )
