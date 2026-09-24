"""Immunopred: immunogenicity risk assessment for protein therapeutics."""

from .predict import ImmunogenicityResult, predict_immunogenicity

__version__ = "0.1.0"

__all__ = [
    "ImmunogenicityResult",
    "predict_immunogenicity",
    "__version__",
]
