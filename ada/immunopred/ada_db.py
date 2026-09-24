"""
Clinical ADA reference database utilities.

Bundles a curated starter table of marketed biologics with published ADA
incidence ranges, and can ingest the public Immunogenicity Database
Collaborative (IDC DS V1) release when downloaded locally.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .sequence import kmer_jaccard, sequence_ratio

DEFAULT_DATA = Path(__file__).resolve().parents[2] / "data" / "ada_reference.csv"


@dataclass
class NeighborHit:
    name: str
    inn: str
    drug_class: str
    ada_rate_pct: float
    similarity: float
    source: str
    method: str = "kmer"
    ada_citation: str = ""
    sequence_citation: str = ""
    nab_rate_pct: float | None = None
    nab_note: str = ""
    nab_note_zh: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "inn": self.inn,
            "drug_class": self.drug_class,
            "ada_rate_pct": self.ada_rate_pct,
            "similarity": self.similarity,
            "source": self.source,
            "method": self.method,
            "ada_citation": self.ada_citation,
            "sequence_citation": self.sequence_citation,
            "nab_rate_pct": self.nab_rate_pct,
            "nab_note": self.nab_note,
            "nab_note_zh": self.nab_note_zh,
        }


class ADAReferenceDB:
    """Load and query clinical ADA priors for protein therapeutics."""

    def __init__(self, df: pd.DataFrame):
        required = {
            "name",
            "inn",
            "drug_class",
            "ada_rate_pct",
            "sequence",
            "source",
        }
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"ADA reference table missing columns: {sorted(missing)}")
        self.df = df.copy()
        self.df["sequence"] = self.df["sequence"].fillna("").astype(str).str.upper()
        for col in (
            "heavy",
            "light",
            "ada_citation",
            "sequence_citation",
            "nab_note",
            "nab_note_zh",
            "source",
        ):
            if col not in self.df.columns:
                self.df[col] = ""
            self.df[col] = self.df[col].fillna("").astype(str)
        if "nab_rate_pct" not in self.df.columns:
            self.df["nab_rate_pct"] = np.nan
        self.df["nab_rate_pct"] = pd.to_numeric(self.df["nab_rate_pct"], errors="coerce")
        for col in ("heavy", "light"):
            self.df[col] = self.df[col].fillna("").astype(str).str.upper()
        self.df["ada_rate_pct"] = pd.to_numeric(self.df["ada_rate_pct"], errors="coerce")
        self.df = self.df.dropna(subset=["ada_rate_pct"])

    @classmethod
    def from_csv(cls, path: str | Path | None = None) -> "ADAReferenceDB":
        path = Path(path) if path else DEFAULT_DATA
        if not path.exists():
            raise FileNotFoundError(
                f"ADA reference file not found: {path}. "
                "Place ada_reference.csv under data/ or pass an explicit path."
            )
        return cls(pd.read_csv(path))

    @classmethod
    def from_idc_directory(cls, idc_dir: str | Path) -> "ADAReferenceDB":
        """
        Build a simplified reference table from a local IDC DS V1 checkout.

        Expects CSV files that include therapeutic metadata, sequences, and
        clinical ADA frequencies. Column names are matched flexibly because
        IDC releases may evolve.
        """
        idc_dir = Path(idc_dir)
        if not idc_dir.exists():
            raise FileNotFoundError(idc_dir)

        csvs = list(idc_dir.rglob("*.csv"))
        if not csvs:
            raise FileNotFoundError(f"No CSV files under {idc_dir}")

        frames = []
        for p in csvs:
            try:
                frames.append(pd.read_csv(p, low_memory=False))
            except Exception:
                continue
        if not frames:
            raise ValueError("Could not read any IDC CSV files")

        big = pd.concat(frames, ignore_index=True, sort=False)

        def pick(*names: str) -> str | None:
            lower = {c.lower(): c for c in big.columns}
            for n in names:
                if n.lower() in lower:
                    return lower[n.lower()]
            # partial contains
            for n in names:
                for lc, orig in lower.items():
                    if n.lower() in lc:
                        return orig
            return None

        name_col = pick("Molecule Assessed for ADA INN Name", "INN", "Therapeutic", "name")
        rate_col = pick(
            "Frequency of ADA+ patients",
            "ADA frequency",
            "ada_rate",
            "Frequency of ADA",
        )
        seq_col = pick("Sequence", "Amino Acid Sequence", "sequence")
        class_col = pick("Modality", "Molecule Type", "drug_class", "Therapeutic Type")
        source_col = pick("Source", "Reference", "source")

        if not name_col or not rate_col:
            raise ValueError(
                "Could not map IDC columns for drug name and ADA frequency. "
                "Inspect the CSV headers and update from_idc_directory()."
            )

        out = pd.DataFrame(
            {
                "name": big[name_col].astype(str),
                "inn": big[name_col].astype(str),
                "drug_class": big[class_col].astype(str) if class_col else "unknown",
                "ada_rate_pct": pd.to_numeric(big[rate_col], errors="coerce"),
                "sequence": big[seq_col].astype(str) if seq_col else "",
                "source": big[source_col].astype(str) if source_col else "IDC-DS-V1",
            }
        )
        # Aggregate multiple cohorts: median ADA rate per molecule
        out = (
            out.dropna(subset=["ada_rate_pct"])
            .groupby(["name", "inn", "drug_class"], as_index=False)
            .agg(
                ada_rate_pct=("ada_rate_pct", "median"),
                sequence=("sequence", "first"),
                source=("source", "first"),
            )
        )
        out["notes"] = "Aggregated median ADA frequency from local IDC export"
        return cls(out)

    def with_sequences(self) -> "ADAReferenceDB":
        return ADAReferenceDB(self.df[self.df["sequence"].str.len() > 20].copy())

    def nearest_neighbors(
        self,
        sequence: str,
        k: int = 5,
        min_similarity: float = 0.05,
    ) -> list[NeighborHit]:
        seq = sequence.upper()
        antibody_like = _looks_like_antibody(seq)
        hits: list[NeighborHit] = []
        for _, row in self.df.iterrows():
            ref = str(row["sequence"])
            heavy = str(row["heavy"])
            light = str(row["light"])
            if len(ref) < 8 and len(heavy) < 20:
                continue
            sim, method = _pair_similarity(seq, heavy, light, ref)
            if antibody_like and method not in {"vh", "vl", "fv"}:
                continue
            floor = 0.85 if method in {"vh", "vl", "fv"} else max(min_similarity, 0.12)
            if sim < floor:
                continue
            nab_raw = row.get("nab_rate_pct")
            nab_val = None if pd.isna(nab_raw) else float(nab_raw)
            hit = NeighborHit(
                name=str(row["name"]),
                inn=str(row["inn"]),
                drug_class=str(row["drug_class"]),
                ada_rate_pct=float(row["ada_rate_pct"]),
                similarity=float(sim),
                source=str(row.get("source", "")),
                method=method,
                ada_citation=str(row.get("ada_citation", "") or ""),
                sequence_citation=str(row.get("sequence_citation", "") or ""),
                nab_rate_pct=nab_val,
                nab_note=str(row.get("nab_note", "") or ""),
                nab_note_zh=str(row.get("nab_note_zh", "") or ""),
            )
            # Keep near-exact marketed hits for display even if blend-trust is tiny
            hits.append(hit)
        hits.sort(key=lambda h: (_trust(h), h.similarity), reverse=True)
        return hits[:k]

    def similarity_weighted_ada(
        self,
        sequence: str,
        k: int = 5,
        min_similarity: float = 0.05,
    ) -> tuple[float | None, list[NeighborHit]]:
        neighbors = self.nearest_neighbors(sequence, k=max(k, 10), min_similarity=min_similarity)
        if not neighbors:
            return None, []
        qlen = len(sequence)
        primary_method = None
        ab_hits = [n for n in neighbors if n.method in {"vh", "vl", "fv"}]
        if ab_hits:
            primary_method = max(ab_hits, key=lambda n: (_trust(n), n.similarity)).method

        scored: list[tuple[float, NeighborHit]] = []
        for n in neighbors:
            w = _trust(n) * _modality_weight(n, qlen, primary_method=primary_method)
            if w <= 0:
                continue
            scored.append((w, n))
        if not scored:
            return None, neighbors[:k]
        scored.sort(key=lambda x: x[0], reverse=True)
        kept = [n for _, n in scored[:k]]
        warr = np.array([w for w, _ in scored[:k]], dtype=float)
        rarr = np.array([n.ada_rate_pct for n in kept], dtype=float)
        return float(np.dot(warr, rarr) / warr.sum()), kept

    def summary(self) -> dict[str, Any]:
        return {
            "n_entries": int(len(self.df)),
            "n_with_sequence": int((self.df["sequence"].str.len() > 20).sum()),
            "ada_rate_median": float(self.df["ada_rate_pct"].median()),
            "ada_rate_mean": float(self.df["ada_rate_pct"].mean()),
            "classes": sorted(self.df["drug_class"].dropna().unique().tolist()),
        }


def _looks_like_antibody(seq: str) -> bool:
    if len(seq) < 90:
        return False
    head = seq[:12]
    heavy_starts = ("EVQL", "QVQL", "QIQL", "DVQL", "AVQL", "QLQL", "QMQL", "QVQV", "QVTL")
    light_starts = (
        "DIQM",
        "EIVL",
        "DIVM",
        "AIQL",
        "DIQL",
        "DVVMT",
        "EIVMT",
        "DIQMT",
        "DIVLT",
        "EIVLT",
        "DILL",
        "DIVL",
        "DIQMT",
        "EIVMT",
        "SSEL",
        "SYEL",
    )
    return any(head.startswith(p) for p in heavy_starts + light_starts)


def _pair_similarity(query: str, heavy: str, light: str, full: str) -> tuple[float, str]:
    """Compare query to VH, VL, and Fv (VH+VL). Prefer the best chain-aware hit."""
    candidates: list[tuple[float, str]] = []
    if len(heavy) >= 80:
        candidates.append((sequence_ratio(query, heavy), "vh"))
        if len(query) > len(heavy) + 40:
            candidates.append((sequence_ratio(query[: len(heavy)], heavy), "vh"))
    if len(light) >= 80:
        candidates.append((sequence_ratio(query, light), "vl"))
        if len(query) > len(light) + 40:
            candidates.append((sequence_ratio(query[-len(light) :], light), "vl"))
            # Also try after VH-sized prefix when H+L concat without constants
            if len(heavy) >= 80 and len(query) >= len(heavy) + len(light) - 5:
                mid = query[len(heavy) : len(heavy) + len(light)]
                candidates.append((sequence_ratio(mid, light), "vl"))
    if len(heavy) >= 80 and len(light) >= 80:
        fv = heavy + light
        candidates.append((sequence_ratio(query, fv), "fv"))
        if len(query) >= len(fv):
            candidates.append((sequence_ratio(query[: len(fv)], fv), "fv"))
    if candidates:
        # Prefer fv over single-chain when similarities tie (clearer marketed-drug call)
        prefer = {"fv": 3, "vh": 2, "vl": 1}
        return max(candidates, key=lambda x: (x[0], prefer.get(x[1], 0)))
    if len(full) < 20:
        return 0.0, "kmer"
    return kmer_jaccard(query, full, k=5), "kmer"


def _trust(hit: NeighborHit) -> float:
    """How much a neighbor should move the clinical prior. 0 = ignore for blending."""
    if hit.method in {"vh", "vl", "fv"}:
        return max(0.0, min(1.0, (hit.similarity - 0.85) / 0.15))
    return max(0.0, min(1.0, (hit.similarity - 0.15) / 0.35))


def _modality_weight(
    hit: NeighborHit,
    query_len: int,
    primary_method: str | None = None,
) -> float:
    """
    Down-weight conjugates / mismatched chains that reuse related sequences
    but have different clinical ADA priors.
    """
    cls = (hit.drug_class or "").lower()
    inn = (hit.inn or "").lower()
    w = 1.0

    # Near-perfect self / same-variable-domain hits: keep clinical prior authoritative
    if hit.similarity >= 0.985 and hit.method in {"vh", "vl", "fv"}:
        return 1.0

    if any(
        k in cls or k in inn
        for k in (
            "conjugate",
            "adc",
            "emtansine",
            "deruxtecan",
            "vedotin",
            "pegylat",
            "pegvisomant",
            "fusion",
            "albumin",
            "antagonist",
            "ctp",
            "lonapeg",
            "somapacitan",
            "somatrogon",
            "efpeg",
        )
    ):
        w *= 0.25 if hit.similarity >= 0.95 else 0.12

    hormone_mod = any(
        k in cls or k in inn
        for k in ("pegylat", "fusion", "albumin", "antagonist", "ctp", "lonapeg", "somatrogon")
    )
    if hormone_mod and (abs(query_len - 191) <= 8 or abs(query_len - 217) <= 8):
        w *= 0.5

    if primary_method in {"vh", "vl", "fv"} and hit.method in {"vh", "vl", "fv"}:
        if primary_method == "vh" and hit.method == "vl" and hit.similarity < 0.98:
            w *= 0.15
        elif primary_method == "vl" and hit.method == "vh" and hit.similarity < 0.98:
            w *= 0.15

    return w
