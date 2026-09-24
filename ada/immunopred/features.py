"""Feature engineering for immunogenicity / ADA-level scoring."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .ada_db import ADAReferenceDB, NeighborHit, _modality_weight, _trust
from .b_cell import BCellResult, predict_b_cell_epitopes
from .cdr import AntibodyRegionMap, annotate_antibody_regions, region_epitope_burden
from .epitope_calibrate import calibrated_ab_epitope_score
from .germline_align import GermlineMismatchReport, analyze_germline_mismatches
from .protein_calibrate import (
    _is_mab_row,
    calibrated_protein_epitope_score,
    infer_protein_class,
)
from .mhc_backend import predict_mhc1 as predict_mhc1_backend
from .mhc_backend import predict_mhc2 as predict_mhc2_backend
from .mhc_i import MHC1Result
from .physchem import PhysChemProfile, profile_physchem
from .sequence import ProteinSequence
from .t_cell import PEPTIDE_LENGTH, TCellResult


@dataclass
class FeatureSet:
    length: int
    t_epitope_density: float
    t_strong_count: int
    t_mean_top_score: float
    mhc1_epitope_density: float
    mhc1_strong_count: int
    mhc1_mean_top_score: float
    b_coverage: float
    b_epitope_count: int
    b_mean_antigenicity: float
    neighbor_ada_pct: float | None
    neighbor_similarity: float
    neighbor_trust: float
    neighbor_nab_pct: float | None
    neighbor_nab_note: str
    neighbor_nab_note_zh: str
    neighbors: list[NeighborHit]
    physchem: PhysChemProfile | None = None
    antibody_regions: AntibodyRegionMap | None = None
    cdr_epitope_density: float | None = None
    framework_epitope_density: float | None = None
    non_germline_burden: float | None = None
    cdr_strong_count: float | None = None
    framework_strong_count: float | None = None
    early_screen_epitope: float | None = None
    ab_epitope_calibration: dict[str, Any] | None = None
    protein_epitope_calibration: dict[str, Any] | None = None
    germline_report: dict[str, Any] | None = None

    def to_vector(self) -> dict[str, float]:
        d = {
            "length": float(self.length),
            "t_epitope_density": self.t_epitope_density,
            "t_strong_count": float(self.t_strong_count),
            "t_mean_top_score": self.t_mean_top_score,
            "mhc1_epitope_density": self.mhc1_epitope_density,
            "mhc1_strong_count": float(self.mhc1_strong_count),
            "mhc1_mean_top_score": self.mhc1_mean_top_score,
            "b_coverage": self.b_coverage,
            "b_epitope_count": float(self.b_epitope_count),
            "b_mean_antigenicity": self.b_mean_antigenicity,
            "neighbor_ada_pct": float(self.neighbor_ada_pct)
            if self.neighbor_ada_pct is not None
            else -1.0,
            "neighbor_similarity": self.neighbor_similarity,
            "neighbor_trust": self.neighbor_trust,
            "neighbor_nab_pct": float(self.neighbor_nab_pct)
            if self.neighbor_nab_pct is not None
            else -1.0,
        }
        if self.physchem is not None:
            d["hydrophobic_frac"] = self.physchem.hydrophobic_frac
            d["charge_density"] = self.physchem.charge_density
            d["agg_risk_hint"] = self.physchem.agg_risk_hint
        if self.non_germline_burden is not None:
            d["non_germline_burden"] = self.non_germline_burden
            d["cdr_epitope_density"] = float(self.cdr_epitope_density or 0.0)
            d["framework_epitope_density"] = float(self.framework_epitope_density or 0.0)
            d["early_screen_epitope"] = float(self.early_screen_epitope or 0.0)
        return d

    def to_dict(self) -> dict[str, Any]:
        d = self.to_vector()
        d["neighbor_nab_note"] = self.neighbor_nab_note
        d["neighbor_nab_note_zh"] = self.neighbor_nab_note_zh
        d["neighbors"] = [n.to_dict() for n in self.neighbors]
        d["physchem"] = None if self.physchem is None else self.physchem.to_dict()
        d["antibody_regions"] = (
            None if self.antibody_regions is None else self.antibody_regions.to_dict()
        )
        d["cdr_strong_count"] = self.cdr_strong_count
        d["framework_strong_count"] = self.framework_strong_count
        d["ab_epitope_calibration"] = self.ab_epitope_calibration
        d["protein_epitope_calibration"] = self.protein_epitope_calibration
        d["germline_report"] = self.germline_report
        return d


def extract_features(
    protein: ProteinSequence,
    ada_db: ADAReferenceDB | None = None,
    t_result: TCellResult | None = None,
    mhc1_result: MHC1Result | None = None,
    b_result: BCellResult | None = None,
    mhc_backend: str | None = None,
    drug_class: str = "",
) -> tuple[FeatureSet, TCellResult, MHC1Result, BCellResult, dict[str, str]]:
    from .mhc_backend import resolve_backend

    requested = resolve_backend(mhc_backend)
    mhc2_meta: dict[str, str] | None = None
    mhc1_meta: dict[str, str] | None = None

    if t_result is None:
        t_result, mhc2_meta = predict_mhc2_backend(protein, backend=mhc_backend)
    if mhc1_result is None:
        mhc1_result, mhc1_meta = predict_mhc1_backend(protein, backend=mhc_backend)
    b_result = b_result or predict_b_cell_epitopes(protein)

    # Prefer MHC-II meta; surface MHC-I fallback if MHC-II succeeded but MHC-I did not.
    if mhc2_meta and mhc1_meta:
        if mhc2_meta.get("fallback_reason"):
            mhc_meta = mhc2_meta
        elif mhc1_meta.get("fallback_reason"):
            mhc_meta = mhc1_meta
        else:
            mhc_meta = mhc2_meta
    elif mhc2_meta:
        mhc_meta = mhc2_meta
    elif mhc1_meta:
        mhc_meta = mhc1_meta
    else:
        mhc_meta = {
            "backend_requested": requested,
            "backend_used": requested,
            "fallback_reason": "",
            "note_zh": "",
            "note_en": "",
        }
    phys = profile_physchem(protein.sequence)
    regions = annotate_antibody_regions(protein.sequence)

    cdr_density = None
    fw_density = None
    non_gl = None
    cdr_n = None
    fw_n = None
    early_epi = None
    ab_cal = None
    prot_cal = None
    # CDR heuristics can false-positive on long enzymes; trust INN/class first.
    pclass_hint = infer_protein_class(protein.id, drug_class, protein.length)
    strong_non_mab = pclass_hint in {
        "insulin",
        "hormone",
        "peptide",
        "cytokine",
        "cytokine_ifn",
        "coagulation",
        "plasma",
        "ert_lysosomal",
        "enzyme_xeno",
        "enzyme_other",
        "fusion",
        "coagulation_xeno",
    }
    treat_as_antibody = bool(
        regions.is_antibody_like
        and not strong_non_mab
        and (
            _is_mab_row(protein.id, "")
            or protein.length <= 480
        )
    )
    gl_report: GermlineMismatchReport | None = None
    gl_dict = None
    if treat_as_antibody:
        starts = list(t_result.strong_starts or [])
        burden = region_epitope_burden(starts, PEPTIDE_LENGTH, regions.cdr_mask)
        cdr_density = burden["cdr_epitope_density"]
        fw_density = burden["framework_epitope_density"]
        non_gl = burden["non_germline_burden"]
        cdr_n = burden["cdr_strong_count"]
        fw_n = burden["framework_strong_count"]
        gl_report = analyze_germline_mismatches(protein.sequence, regions)
        gl_dict = gl_report.to_dict() if gl_report else None
        b_spans = [(e.start, e.end) for e in (b_result.epitopes or [])]
        cal = calibrated_ab_epitope_score(
            sequence=protein.sequence,
            regions=regions,
            non_germline_burden=float(non_gl),
            charge_density=float(phys.charge_density if phys else 0.0),
            agg_risk_hint=float(phys.agg_risk_hint if phys else 0.0),
            inn=protein.id,
            drug_class=drug_class or "",
            b_coverage=float(b_result.coverage),
            cdr_epitope_density=float(cdr_density),
            germline_report=gl_report,
            t_strong_starts=starts,
            b_epitope_spans=b_spans,
        )
        ab_cal = cal.to_dict()
        early_epi = float(cal.score)
    else:
        pcal = calibrated_protein_epitope_score(
            protein.sequence,
            inn=protein.id,
            drug_class=drug_class or "",
            t_epitope_density=float(t_result.epitope_density),
            t_mean_top_score=float(t_result.mean_top_score),
            b_coverage=float(b_result.coverage),
            mhc1_epitope_density=float(mhc1_result.epitope_density),
            charge_density=float(phys.charge_density if phys else 0.0),
            agg_risk_hint=float(phys.agg_risk_hint if phys else 0.0),
        )
        prot_cal = pcal.to_dict()
        early_epi = float(pcal.score)

    neighbor_ada = None
    neighbor_sim = 0.0
    neighbor_trust = 0.0
    neighbor_nab = None
    neighbor_nab_note = ""
    neighbor_nab_note_zh = ""
    neighbors: list[NeighborHit] = []
    if ada_db is not None:
        neighbor_ada, neighbors = ada_db.similarity_weighted_ada(protein.sequence)
        if neighbors:
            neighbor_sim = neighbors[0].similarity
            primary = neighbors[0].method if neighbors[0].method in {"vh", "vl", "fv"} else None
            neighbor_trust = max(
                _trust(n) * _modality_weight(n, protein.length, primary_method=primary)
                for n in neighbors
            )
            nab_vals = []
            nab_w = []
            notes_en = []
            notes_zh = []
            for n in neighbors:
                w = _trust(n) * _modality_weight(n, protein.length, primary_method=primary)
                if w <= 0:
                    continue
                if n.nab_rate_pct is not None and not (
                    isinstance(n.nab_rate_pct, float) and np.isnan(n.nab_rate_pct)
                ):
                    nab_vals.append(float(n.nab_rate_pct))
                    nab_w.append(w)
                    if n.nab_note:
                        notes_en.append(f"{n.inn}: {n.nab_note}")
                    if n.nab_note_zh:
                        notes_zh.append(f"{n.inn}: {n.nab_note_zh}")
            if nab_vals:
                neighbor_nab = float(np.average(nab_vals, weights=nab_w))
                neighbor_nab_note = " | ".join(notes_en[:2])
                neighbor_nab_note_zh = " | ".join(notes_zh[:2])

    # Refine non-mAb class using closest marketed neighbor's drug_class when INN is opaque
    if (
        prot_cal is not None
        and neighbors
        and neighbor_sim >= 0.85
        and str(prot_cal.get("protein_class") or "") in {"other", "peptide", "fusion", "enzyme_other"}
    ):
        hint = f"{neighbors[0].inn} {neighbors[0].drug_class}"
        pcal = calibrated_protein_epitope_score(
            protein.sequence,
            inn=protein.id or str(neighbors[0].inn or ""),
            drug_class=hint,
            t_epitope_density=float(t_result.epitope_density),
            t_mean_top_score=float(t_result.mean_top_score),
            b_coverage=float(b_result.coverage),
            mhc1_epitope_density=float(mhc1_result.epitope_density),
            charge_density=float(phys.charge_density if phys else 0.0),
            agg_risk_hint=float(phys.agg_risk_hint if phys else 0.0),
        )
        prot_cal = pcal.to_dict()
        early_epi = float(pcal.score)

    features = FeatureSet(
        length=protein.length,
        t_epitope_density=t_result.epitope_density,
        t_strong_count=t_result.strong_count,
        t_mean_top_score=t_result.mean_top_score,
        mhc1_epitope_density=mhc1_result.epitope_density,
        mhc1_strong_count=mhc1_result.strong_count,
        mhc1_mean_top_score=mhc1_result.mean_top_score,
        b_coverage=b_result.coverage,
        b_epitope_count=len(b_result.epitopes),
        b_mean_antigenicity=b_result.mean_antigenicity,
        neighbor_ada_pct=neighbor_ada,
        neighbor_similarity=neighbor_sim,
        neighbor_trust=neighbor_trust,
        neighbor_nab_pct=neighbor_nab,
        neighbor_nab_note=neighbor_nab_note,
        neighbor_nab_note_zh=neighbor_nab_note_zh,
        neighbors=neighbors,
        physchem=phys,
        antibody_regions=regions,
        cdr_epitope_density=cdr_density,
        framework_epitope_density=fw_density,
        non_germline_burden=non_gl,
        cdr_strong_count=cdr_n,
        framework_strong_count=fw_n,
        early_screen_epitope=early_epi,
        ab_epitope_calibration=ab_cal,
        protein_epitope_calibration=prot_cal,
        germline_report=gl_dict,
    )
    return features, t_result, mhc1_result, b_result, mhc_meta
