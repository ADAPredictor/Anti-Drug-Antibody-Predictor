"""
Predictive immunogenicity score for non-antibody protein drugs.

Clinical ADA is an *outcome* — not used to correct the prediction.
Pipeline:
  1) Foreignness (sequence vs human self pool + soft class prior)
  2) T/B epitope burden (sequence)
  3) Product-type tolerance multiplier (host central/peripheral tolerance)
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

import numpy as np


# Foreignness 0–100: how non-self the product is expected to look to the immune system
_FOREIGNNESS = {
    "insulin": 8.0,
    "plasma": 12.0,
    "hormone": 18.0,
    "peptide": 35.0,
    "cytokine": 40.0,
    "cytokine_ifn": 48.0,
    "coagulation": 18.0,  # same near-self prior as hormones
    "coagulation_xeno": 78.0,  # non-human species (e.g. porcine FVIII)
    "fusion": 42.0,
    "enzyme_other": 55.0,
    "ert_lysosomal": 88.0,  # human-like seq, but missing-self foreignness tag
    "enzyme_xeno": 90.0,
    "other": 50.0,
}

# Tolerance multiplier: <1 = more host tolerance (dampen); >1 = less tolerance (amplify)
_TOLERANCE = {
    "insulin": 0.22,
    "plasma": 0.28,
    "hormone": 0.35,
    "peptide": 0.55,
    "cytokine": 0.85,
    "cytokine_ifn": 0.95,
    "coagulation": 0.35,  # same host tolerance as hormones
    "coagulation_xeno": 1.15,
    "fusion": 0.85,
    "enzyme_other": 0.95,
    "ert_lysosomal": 1.40,
    "enzyme_xeno": 1.25,
    "other": 1.0,
}


@dataclass
class ProteinEpitopeBreakdown:
    score: float
    protein_class: str
    foreignness: float
    epitope_raw: float
    tolerance_mult: float
    gate: float
    # Legacy aliases kept for UI/compat (not clinical corrections)
    class_prior: float
    knn_clinical: float
    knn_n: int
    raw_epitope: float
    note_zh: str
    note_en: str
    homology: dict[str, Any] | None = None
    class_foreignness: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "score": round(self.score, 3),
            "protein_class": self.protein_class,
            "foreignness": round(self.foreignness, 3),
            "epitope_raw": round(self.epitope_raw, 3),
            "tolerance_mult": round(self.tolerance_mult, 3),
            "gate": round(self.gate, 3),
            "class_prior": round(self.class_prior, 3),
            "knn_clinical": round(self.knn_clinical, 3),
            "knn_n": self.knn_n,
            "raw_epitope": round(self.raw_epitope, 3),
            "note_zh": self.note_zh,
            "note_en": self.note_en,
            "homology": self.homology,
            "class_foreignness": (
                round(float(self.class_foreignness), 3)
                if self.class_foreignness is not None
                else None
            ),
        }


def _stem(inn: str) -> str:
    t = (inn or "").strip().lower().replace("-", " ")
    if not t:
        return ""
    return re.sub(r"[^a-z0-9]+", "", t.split()[0])


def infer_protein_class(inn: str = "", drug_class: str = "", length: int = 0) -> str:
    """Heuristic class label for foreignness / tolerance tags."""
    s = f"{inn} {drug_class}".lower()
    L = int(length or 0)

    if re.search(r"\binsulin\b", s):
        return "insulin"
    if re.search(
        r"somatropin|pegvisomant|somatrogon|lonapeg|growth.?hormone|hormone-gh|"
        r"metreleptin|tesamorelin|glucagon|igf|pth|teriparatide",
        s,
    ):
        return "hormone"
    if re.search(
        r"liraglutide|exenatide|semaglutide|dulaglutide|glp.?1|pramlintide",
        s,
    ):
        return "peptide"
    if re.search(
        r"interferon|aldesleukin|oprelvekin|interleukin-2|\bil-?2\b|tasonermin",
        s,
    ):
        return "cytokine_ifn"
    if re.search(
        r"epoetin|filgrastim|pegfilgrastim|cytokine|"
        r"g-?csf|interleukin|darbepoetin|molgramostim|sargramostim|"
        r"anakinra|becaplermin|palifermin|balugrastim|eflapegrastim",
        s,
    ):
        return "cytokine"
    if re.search(
        r"factor\s*viii|octocog|moroctocog|efmoroctocog|antihemophilic|"
        r"factor\s*ix|nonacog|eptacog|factor\s*vii|coagulation|susoctocog|"
        r"turoctocog|lonoctocog|damasectocog|von willebrand|catridecacog|"
        r"albutrepenonacog|eftrenonacog|efanesoctocog|rurioctocog|simoctocog|"
        r"damoctocog",
        s,
    ):
        # Non-human species clotting factor (structure/origin — not clinical rate)
        if re.search(r"porcine|pig\b|bovine|murine|canine", s):
            return "coagulation_xeno"
        return "coagulation"
    if re.search(r"albumin|antithrombin|conestat|c1 esterase|plasma-protein", s):
        return "plasma"
    if re.search(
        r"glucosidase|galsulfase|laronidase|idursulfase|agalsidase|vestronidase|"
        r"elosulfase|asfotase|sebelipase|imiglucerase|taliglucerase|velaglucerase|"
        r"cerliponase|avalglucosidase|enzyme-ert|\bert\b",
        s,
    ):
        return "ert_lysosomal"
    if re.search(
        r"rasburicase|pegloticase|asparaginase|crisantaspase|streptokinase|"
        r"botulinum|urate-oxidase|imlifidase|\bides\b|ide.?s\b|ricefidase|"
        r"kj.?103|pegvaliase|glucarpidase|bacterial|streptococ|pathogen.?enzyme|"
        r"endos|endo.?s2?",
        s,
    ):
        return "enzyme_xeno"
    if re.search(
        r"alteplase|reteplase|tenecteplase|urokinase|hyaluronidase|pancrelipase|"
        r"collagenase|tilactase|enzyme-tpa",
        s,
    ):
        return "enzyme_other"
    if re.search(r"fusion|cept\b|etanercept|rilonacept|belatacept|luspatercept", s):
        return "fusion"
    if re.search(r"enzyme", s):
        return "enzyme_other"
    if L and L < 80:
        return "peptide"
    return "other"


def _raw_epitope_score(
    t_density: float,
    t_mean_top: float,
    b_coverage: float,
    mhc1_density: float,
    agg_hint: float,
) -> float:
    t_density_term = float(np.clip(t_density / 28.0 * 100.0, 0, 100))
    t_score_term = float(np.clip(t_mean_top / 2.5 * 100.0, 0, 100))
    b_cov_term = float(np.clip(b_coverage * 100.0, 0, 100))
    mhc1_term = float(np.clip(mhc1_density / 40.0 * 100.0, 0, 100))
    agg_term = float(np.clip(agg_hint, 0, 100))
    return float(
        0.45 * t_density_term
        + 0.20 * t_score_term
        + 0.20 * b_cov_term
        + 0.08 * mhc1_term
        + 0.07 * agg_term
    )


def _is_mab_row(inn: str, drug_class: str) -> bool:
    cls = (drug_class or "").lower()
    token = (inn or "").lower().replace("-", " ").split()[0] if inn else ""
    if re.search(r"antibody|monoclonal|\bmab\b|igg|chimeric|humanized|adc|bispecific", cls):
        return True
    return bool(re.search(r"(mab|zumab|ximab|umab)$", token))


def _foreignness_score(pclass: str, inn: str = "", drug_class: str = "") -> float:
    base = float(_FOREIGNNESS.get(pclass, _FOREIGNNESS["other"]))
    blob = f"{inn} {drug_class}".lower()
    if re.search(
        r"bacterial|streptococ|imlifidase|\bides\b|endos|pathogen|asparaginase|"
        r"pegvaliase|glucarpidase|rasburicase|pegloticase",
        blob,
    ):
        base = max(base, 92.0)
    if re.search(r"exenatide|exendin", blob):
        base = max(base, 70.0)
    if re.search(r"metreleptin|\bleptin\b", blob):
        base = max(base, 55.0)
    return float(np.clip(base, 0.0, 100.0))


def _tolerance_mult(pclass: str, inn: str = "", drug_class: str = "") -> float:
    mult = float(_TOLERANCE.get(pclass, _TOLERANCE["other"]))
    blob = f"{inn} {drug_class}".lower()
    # Missing endogenous protein / pathogen enzymes: lower host tolerance
    if re.search(r"enzyme-ert|\bert\b|bacterial|streptococ|imlifidase|endos", blob):
        mult = max(mult, 1.20)
    if re.search(r"exenatide|exendin", blob):
        mult = max(mult, 1.05)
    return float(np.clip(mult, 0.15, 1.50))


def calibrated_protein_epitope_score(
    sequence: str,
    *,
    inn: str = "",
    drug_class: str = "",
    t_epitope_density: float = 0.0,
    t_mean_top_score: float = 0.0,
    b_coverage: float = 0.0,
    mhc1_epitope_density: float = 0.0,
    charge_density: float = 0.0,
    agg_risk_hint: float = 0.0,
) -> ProteinEpitopeBreakdown:
    """Predictive protein score: homology foreignness × epitope × tolerance."""
    from .human_homology import best_human_homology

    del charge_density  # reserved for future physchem term
    L = len(sequence or "")
    pclass = infer_protein_class(inn, drug_class, L)
    epitope = _raw_epitope_score(
        t_epitope_density,
        t_mean_top_score,
        b_coverage,
        mhc1_epitope_density,
        agg_risk_hint,
    )
    class_f = _foreignness_score(pclass, inn, drug_class)
    hom = best_human_homology(sequence)
    # Sequence vs human self pool is primary; class label is a soft prior / pathogen floor
    if hom.n_refs >= 8 and hom.best_identity > 0:
        foreignness = float(
            np.clip(0.70 * float(hom.foreignness) + 0.30 * float(class_f), 0.0, 100.0)
        )
        # ERT missing-self tag: human homology must not cancel foreignness
        if pclass == "ert_lysosomal":
            foreignness = max(foreignness, float(class_f))
        # Pathogen / bacterial enzymes: never let homology accidental hits collapse risk
        elif pclass == "enzyme_xeno" or class_f >= 85.0:
            foreignness = max(foreignness, float(class_f) * 0.92, float(hom.foreignness))
        # Non-human species clotting factor: keep species-origin floor
        # (DB sequences may be human placeholders — name/species is the structure cue)
        elif pclass == "coagulation_xeno":
            foreignness = max(foreignness, float(class_f) * 0.90)
        # Human coagulation: trust homology (near-self); do not inflate from clinical inhibitors
        elif pclass == "coagulation" and hom.best_identity >= 0.85:
            foreignness = min(foreignness, max(float(hom.foreignness), 18.0))
    else:
        foreignness = float(class_f)

    tolerance = _tolerance_mult(pclass, inn, drug_class)
    # Gate: self-like sequences contribute little even if MHC density is high
    gate = 0.20 + 0.80 * (foreignness / 100.0)
    score = float(np.clip(epitope * gate * tolerance, 1.0, 95.0))

    return ProteinEpitopeBreakdown(
        score=score,
        protein_class=pclass,
        foreignness=foreignness,
        epitope_raw=epitope,
        tolerance_mult=tolerance,
        gate=gate,
        class_prior=foreignness,  # UI legacy: show foreignness where "class prior" was
        knn_clinical=0.0,
        knn_n=0,
        raw_epitope=epitope,
        homology=hom.to_dict(),
        class_foreignness=class_f,
        note_zh=(
            f"预测链路：异源度={foreignness:.0f}"
            f"（序列对人源库={hom.foreignness:.0f}/相似{hom.best_identity:.2f}，"
            f"类别先验={class_f:.0f}/{pclass}"
            f"{'，ERT内源缺失异源标签' if pclass == 'ert_lysosomal' else ''}）→ "
            f"门控={gate:.2f}；T/B 表位负担={epitope:.1f}；"
            f"产品耐受系数={tolerance:.2f} → 分={score:.1f}。"
            f"{hom.method_note_zh}"
            "临床 ADA 为结果对照，不参与矫正。"
        ),
        note_en=(
            f"Predictive path: foreignness={foreignness:.0f}"
            f" (seq-vs-human={hom.foreignness:.0f}/sim{hom.best_identity:.2f}, "
            f"class prior={class_f:.0f}/{pclass}) → "
            f"gate={gate:.2f}; T/B epitope burden={epitope:.1f}; "
            f"tolerance mult={tolerance:.2f} → score={score:.1f}. "
            f"{hom.method_note_en} "
            "Clinical ADA is an outcome reference, not a correction."
        ),
    )
