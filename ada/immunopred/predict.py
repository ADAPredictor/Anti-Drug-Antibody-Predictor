"""End-to-end immunogenicity / ADA-level / NAb-risk prediction."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .ada_db import ADAReferenceDB
from .b_cell import BCellResult
from .context_notes import build_context_notes
from .features import extract_features
from .mhc_backend import meta_for_result
from .mhc_i import MHC1Result
from .model import RiskScore, score_risk
from .nab import NabScore, score_nab
from .sequence import ProteinSequence
from .t_cell import TCellResult


@dataclass
class ImmunogenicityResult:
    sequence_id: str
    length: int
    risk: RiskScore
    nab: NabScore
    t_cell: TCellResult  # MHC-II
    mhc1: MHC1Result
    b_cell: BCellResult
    features: dict[str, Any]
    context_notes: list[dict[str, Any]]
    mhc_backend_meta: dict[str, str] = field(default_factory=dict)
    disclaimer_zh: str = (
        "仅供研究。MHC 结合默认为 Heuristic；可选 IEDB / 本地 NetMHC，失败时回退 Heuristic。"
        "不能代替体外 MHC 实验、ADA/NAb 生物分析或监管免疫原性评估。"
    )
    disclaimer_en: str = (
        "Research use only. MHC binders default to heuristic; optional IEDB / local NetMHC "
        "with auto-fallback. Not a substitute for in vitro MHC assays, ADA/NAb bioanalysis, "
        "or regulatory immunogenicity assessment."
    )

    def to_dict(self) -> dict[str, Any]:
        return {
            "sequence_id": self.sequence_id,
            "length": self.length,
            "risk": self.risk.to_dict(),
            "ada": self.risk.to_dict(),
            "nab": self.nab.to_dict(),
            "mhc2": self.t_cell.to_dict(),
            "mhc1": self.mhc1.to_dict(),
            "t_cell": self.t_cell.to_dict(),
            "b_cell": self.b_cell.to_dict(),
            "features": self.features,
            "context_notes": self.context_notes,
            "mhc_backend": meta_for_result(self.mhc_backend_meta)
            if self.mhc_backend_meta
            else {
                "requested": "heuristic",
                "used": "heuristic",
                "fallback_reason": "",
                "note_zh": "",
                "note_en": "",
            },
            "disclaimer": self.disclaimer_en,
            "disclaimer_zh": self.disclaimer_zh,
            "disclaimer_en": self.disclaimer_en,
        }


def predict_immunogenicity(
    protein: ProteinSequence,
    ada_db: ADAReferenceDB | None = None,
    use_builtin_db: bool = True,
    mhc_backend: str | None = None,
    drug_class: str = "",
) -> ImmunogenicityResult:
    """Predict ADA level and neutralizing-ADA risk for a protein therapeutic."""
    if ada_db is None and use_builtin_db:
        try:
            ada_db = ADAReferenceDB.from_csv()
        except FileNotFoundError:
            ada_db = None

    features, t_result, mhc1_result, b_result, mhc_meta = extract_features(
        protein, ada_db=ada_db, mhc_backend=mhc_backend, drug_class=drug_class or ""
    )
    risk = score_risk(features)
    nab = score_nab(features, ada_score=risk.score)
    context = [n.to_dict() for n in build_context_notes(features)]
    return ImmunogenicityResult(
        sequence_id=protein.id,
        length=protein.length,
        risk=risk,
        nab=nab,
        t_cell=t_result,
        mhc1=mhc1_result,
        b_cell=b_result,
        features=features.to_dict(),
        context_notes=context,
        mhc_backend_meta=mhc_meta,
    )
