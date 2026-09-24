"""Amino-acid physicochemical scales used by classical epitope predictors."""

from __future__ import annotations

# Kolaskar & Tongaonkar (1990) antigenicity
KOLASKAR = {
    "A": 1.064,
    "R": 0.873,
    "N": 0.776,
    "D": 0.866,
    "C": 1.412,
    "Q": 1.015,
    "E": 0.851,
    "G": 0.874,
    "H": 1.105,
    "I": 1.152,
    "L": 1.250,
    "K": 0.930,
    "M": 0.826,
    "F": 1.091,
    "P": 1.055,
    "S": 1.012,
    "T": 0.909,
    "W": 0.893,
    "Y": 1.161,
    "V": 1.383,
}

# Parker hydrophilicity (Hopp/Woods-like Parker 1986)
PARKER = {
    "A": 2.1,
    "R": 4.2,
    "N": 7.0,
    "D": 10.0,
    "C": -1.4,
    "Q": 6.0,
    "E": 7.8,
    "G": 5.7,
    "H": 2.1,
    "I": -8.0,
    "L": -9.2,
    "K": 5.7,
    "M": -4.2,
    "F": -9.2,
    "P": 2.1,
    "S": 6.5,
    "T": 5.2,
    "W": -10.0,
    "Y": -1.9,
    "V": -3.7,
}

# Emini surface accessibility / fractional surface probability factors
EMINI = {
    "A": 0.49,
    "R": 0.95,
    "N": 0.78,
    "D": 0.78,
    "C": 0.26,
    "Q": 0.84,
    "E": 0.84,
    "G": 0.48,
    "H": 0.66,
    "I": 0.34,
    "L": 0.40,
    "K": 0.97,
    "M": 0.48,
    "F": 0.42,
    "P": 0.75,
    "S": 0.65,
    "T": 0.70,
    "W": 0.51,
    "Y": 0.76,
    "V": 0.36,
}

# Chou-Fasman beta-turn propensity (approximate)
CHOU_FASMAN_TURN = {
    "A": 0.66,
    "R": 0.95,
    "N": 1.56,
    "D": 1.46,
    "C": 1.19,
    "Q": 0.98,
    "E": 0.74,
    "G": 1.56,
    "H": 0.95,
    "I": 0.47,
    "L": 0.59,
    "K": 1.01,
    "M": 0.60,
    "F": 0.60,
    "P": 1.52,
    "S": 1.43,
    "T": 0.96,
    "W": 0.96,
    "Y": 1.14,
    "V": 0.50,
}

# Kyte-Doolittle hydrophobicity (for MHC pocket heuristics)
KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}


def scale_value(aa: str, scale: dict[str, float], default: float = 0.0) -> float:
    return scale.get(aa.upper(), default)


def window_scores(sequence: str, scale: dict[str, float], window: int) -> list[float]:
    """Centered sliding-window mean of a physicochemical scale."""
    n = len(sequence)
    if n == 0:
        return []
    if window > n:
        window = n
    half = window // 2
    vals = [scale_value(aa, scale) for aa in sequence]
    out: list[float] = []
    for i in range(n):
        start = max(0, i - half)
        end = min(n, start + window)
        start = max(0, end - window)
        out.append(sum(vals[start:end]) / (end - start))
    return out
