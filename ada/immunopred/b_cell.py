"""
Linear B-cell epitope prediction using classical propensity scales.

Methods combined (research heuristics, not crystal-structure discontinuous
epitope prediction):
  - Kolaskar-Tongaonkar antigenicity
  - Parker hydrophilicity
  - Emini surface accessibility
  - Chou-Fasman beta-turn propensity
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from .scales import (
    CHOU_FASMAN_TURN,
    EMINI,
    KOLASKAR,
    PARKER,
    window_scores,
)
from .sequence import ProteinSequence


@dataclass
class BCellEpitope:
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    sequence: str
    score: float
    method: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class BCellResult:
    epitopes: list[BCellEpitope]
    coverage: float  # fraction of residues in any predicted epitope
    mean_antigenicity: float
    residue_scores: list[float]

    def to_dict(self) -> dict[str, Any]:
        return {
            "coverage": self.coverage,
            "mean_antigenicity": self.mean_antigenicity,
            "epitopes": [e.to_dict() for e in self.epitopes],
        }


def _zscore(xs: list[float]) -> np.ndarray:
    arr = np.asarray(xs, dtype=float)
    if arr.size == 0:
        return arr
    mu = arr.mean()
    sigma = arr.std()
    if sigma < 1e-9:
        return np.zeros_like(arr)
    return (arr - mu) / sigma


def _merge_runs(mask: np.ndarray, min_len: int = 6) -> list[tuple[int, int]]:
    """Return 0-based inclusive (start, end) runs where mask is True."""
    runs: list[tuple[int, int]] = []
    n = len(mask)
    i = 0
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and mask[j + 1]:
            j += 1
        if j - i + 1 >= min_len:
            runs.append((i, j))
        i = j + 1
    return runs


def predict_b_cell_epitopes(
    protein: ProteinSequence,
    window: int = 7,
    min_epitope_len: int = 6,
    score_threshold: float = 0.5,
    top_n: int = 20,
) -> BCellResult:
    """
    Predict linear B-cell epitope stretches.

    Residue-level consensus score is the mean of z-scored Kolaskar, Parker,
    Emini, and Chou-Fasman turn scales. Contiguous stretches above
    ``score_threshold`` (in z-score units) are reported as epitopes.
    """
    seq = protein.sequence
    if not seq:
        return BCellResult([], 0.0, 0.0, [])

    kolaskar = window_scores(seq, KOLASKAR, window)
    parker = window_scores(seq, PARKER, window)
    emini = window_scores(seq, EMINI, window)
    turns = window_scores(seq, CHOU_FASMAN_TURN, window)

    consensus = (
        _zscore(kolaskar)
        + _zscore(parker)
        + _zscore(emini)
        + _zscore(turns)
    ) / 4.0

    mask = consensus >= score_threshold
    runs = _merge_runs(mask, min_len=min_epitope_len)

    epitopes: list[BCellEpitope] = []
    covered = np.zeros(len(seq), dtype=bool)
    for start0, end0 in runs:
        covered[start0 : end0 + 1] = True
        sub = seq[start0 : end0 + 1]
        score = float(consensus[start0 : end0 + 1].mean())
        epitopes.append(
            BCellEpitope(
                start=start0 + 1,
                end=end0 + 1,
                sequence=sub,
                score=score,
                method="consensus_KT_Parker_Emini_CF",
            )
        )

    epitopes.sort(key=lambda e: e.score, reverse=True)
    epitopes = epitopes[:top_n]
    coverage = float(covered.mean()) if len(seq) else 0.0

    return BCellResult(
        epitopes=epitopes,
        coverage=coverage,
        mean_antigenicity=float(np.mean(kolaskar)) if kolaskar else 0.0,
        residue_scores=consensus.tolist(),
    )
