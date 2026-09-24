"""
MHC class I (CD8 T-cell) epitope heuristics for protein therapeutics.

ADA formation is primarily CD4 / MHC-II driven; MHC-I binders are reported for
completeness (cytotoxic T-cell risk, vaccine-style epitope maps). Matrices are
research heuristics, not NetMHCpan.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from .sequence import ProteinSequence

DEFAULT_MHC1_ALLELES = [
    "HLA-A*02:01",
    "HLA-A*01:01",
    "HLA-A*03:01",
    "HLA-A*24:02",
    "HLA-B*07:02",
    "HLA-B*08:01",
    "HLA-B*44:02",
]

PEPTIDE_LENGTH = 9


@dataclass
class MHC1Epitope:
    start: int
    peptide: str
    allele: str
    score: float
    percentile: float
    strong: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class MHC1Result:
    epitopes: list[MHC1Epitope]
    strong_count: int
    epitope_density: float
    mean_top_score: float
    allele_summary: dict[str, dict[str, float]]

    def to_dict(self) -> dict[str, Any]:
        return {
            "strong_count": self.strong_count,
            "epitope_density": self.epitope_density,
            "mean_top_score": self.mean_top_score,
            "allele_summary": self.allele_summary,
            "epitopes": [e.to_dict() for e in self.epitopes],
        }


def _empty() -> dict[str, list[float]]:
    return {a: [0.0] * 9 for a in "ACDEFGHIKLMNPQRSTVWY"}


def _mhc1_matrices() -> dict[str, dict[str, list[float]]]:
    """Allele-specific 9-mer anchor preferences (positions 1..9 → index 0..8)."""
    mats: dict[str, dict[str, list[float]]] = {}

    def make(allele: str, anchors: dict[int, set[str]], anti: dict[int, set[str]] | None = None):
        mat = _empty()
        for pos, residues in anchors.items():
            for aa in residues:
                mat[aa][pos] += 1.4
        if anti:
            for pos, residues in anti.items():
                for aa in residues:
                    mat[aa][pos] -= 0.8
        # Mild preference for hydrophobic interiors
        for aa in "FILMWYV":
            mat[aa][2] += 0.2
            mat[aa][5] += 0.2
        mats[allele] = mat

    make("HLA-A*02:01", {1: set("LM"), 8: set("VLIM")}, {1: set("DEKR"), 8: set("DEKR")})
    make("HLA-A*01:01", {1: set("YTS"), 8: set("Y")}, {8: set("DEKR")})
    make("HLA-A*03:01", {1: set("AILMV"), 8: set("KR")}, {8: set("DE")})
    make("HLA-A*24:02", {1: set("YF"), 8: set("FILMWY")}, {1: set("DEKR")})
    make("HLA-B*07:02", {1: set("P"), 8: set("LMF")}, {1: set("DEKR")})
    make("HLA-B*08:01", {2: set("K"), 8: set("L")}, {2: set("DE")})
    make("HLA-B*44:02", {1: set("E"), 8: set("FY")}, {1: set("KR")})
    return mats


_MATRICES = _mhc1_matrices()


def _score_peptide(pep: str, pssm: dict[str, list[float]]) -> float:
    return float(sum(pssm.get(aa, [0.0] * 9)[i] for i, aa in enumerate(pep)) / len(pep))


def predict_mhc1_epitopes(
    protein: ProteinSequence,
    alleles: list[str] | None = None,
    top_n: int = 30,
    strong_percentile: float = 10.0,
) -> MHC1Result:
    alleles = alleles or DEFAULT_MHC1_ALLELES
    peptides = protein.peptides(PEPTIDE_LENGTH, step=1)
    all_hits: list[MHC1Epitope] = []
    allele_summary: dict[str, dict[str, float]] = {}

    for allele in alleles:
        if allele not in _MATRICES:
            raise KeyError(f"No MHC-I matrix for {allele}")
        pssm = _MATRICES[allele]
        scored = [(start, pep, _score_peptide(pep, pssm)) for start, pep in peptides]
        if not scored:
            allele_summary[allele] = {
                "n_peptides": 0,
                "strong": 0,
                "mean_score": 0.0,
                "max_score": 0.0,
            }
            continue
        scores = np.array([s[2] for s in scored], dtype=float)
        order = scores.argsort()
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.linspace(0, 100, len(scores))
        strong = 0
        for (start, pep, score), pct in zip(scored, ranks):
            is_strong = pct >= (100.0 - strong_percentile) and score >= 0.15
            if is_strong:
                strong += 1
            all_hits.append(
                MHC1Epitope(
                    start=start,
                    peptide=pep,
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
    unique_starts = {e.start for e in strong_epitopes}
    density = (len(unique_starts) / max(protein.length, 1)) * 100.0
    top = strong_epitopes[:top_n]
    mean_top = float(np.mean([e.score for e in top])) if top else 0.0
    return MHC1Result(
        epitopes=top,
        strong_count=len(unique_starts),
        epitope_density=density,
        mean_top_score=mean_top,
        allele_summary=allele_summary,
    )
