"""Simple physicochemical sequence descriptors (research heuristics)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

HYDROPHOBIC = set("AILMFWYV")
POSITIVE = set("KR")
NEGATIVE = set("DE")
AROMATIC = set("FYW")


@dataclass
class PhysChemProfile:
    length: int
    hydrophobic_frac: float
    charge_density: float  # (K+R − D−E) / length
    aromatic_frac: float
    cysteine_count: int
    agg_risk_hint: float  # 0–100 heuristic, NOT clinical ADA

    def to_dict(self) -> dict[str, Any]:
        return {
            "length": self.length,
            "hydrophobic_frac": round(self.hydrophobic_frac, 4),
            "charge_density": round(self.charge_density, 4),
            "aromatic_frac": round(self.aromatic_frac, 4),
            "cysteine_count": self.cysteine_count,
            "agg_risk_hint": round(self.agg_risk_hint, 2),
        }


def profile_physchem(sequence: str) -> PhysChemProfile:
    seq = "".join(c for c in sequence.upper() if c.isalpha())
    n = max(len(seq), 1)
    hydro = sum(1 for c in seq if c in HYDROPHOBIC) / n
    charge = (sum(1 for c in seq if c in POSITIVE) - sum(1 for c in seq if c in NEGATIVE)) / n
    aromatic = sum(1 for c in seq if c in AROMATIC) / n
    cys = seq.count("C")
    # Crude aggregation-related hint: high hydrophobicity + low |charge|
    agg = 100.0 * (0.65 * max(0.0, hydro - 0.35) / 0.25 + 0.35 * max(0.0, 0.08 - abs(charge)) / 0.08)
    agg = max(0.0, min(100.0, agg))
    return PhysChemProfile(
        length=len(seq),
        hydrophobic_frac=hydro,
        charge_density=charge,
        aromatic_frac=aromatic,
        cysteine_count=cys,
        agg_risk_hint=agg,
    )
