"""
Sequence-level foreignness vs human self / endogenous reference pool.

Replaces pure class-label foreignness for non-Ab proteins when a usable
human reference pool is available. Class labels remain a soft prior / floor
for pathogen enzymes and a tolerance driver (e.g. ERT CRIM-negative).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from .sequence import kmer_jaccard, sequence_ratio

_ROOT = Path(__file__).resolve().parents[2]
_EXTRA_CSV = _ROOT / "data" / "extra_protein_sequences.csv"
_UNIPROT_DIR = _ROOT / "data" / "cache_uniprot"

# drug_class / notes tokens that are NOT reliable human self references
_XENO_RE = (
    "bacterial",
    "streptococ",
    "escherichia",
    "erwinia",
    "botulinum",
    "exendin",
    "salmon",
    "urate-oxidase",
    "asparaginase",
    "pathogen",
    "imlifidase",
    "ides",
    "endos",
    "pegvaliase",
    "glucarpidase",
    "rasburicase",
    "pegloticase",
    "crisantaspase",
    "hirudin",
    "glatiramer",
)


@dataclass
class HomologyReport:
    best_identity: float
    best_label: str
    foreignness: float  # 0–100 from 1 − identity
    n_refs: int
    method_note_zh: str
    method_note_en: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "best_identity": round(float(self.best_identity), 4),
            "best_label": self.best_label,
            "foreignness": round(float(self.foreignness), 3),
            "n_refs": int(self.n_refs),
            "method_note_zh": self.method_note_zh,
            "method_note_en": self.method_note_en,
        }


def _is_xeno_blob(blob: str) -> bool:
    s = (blob or "").lower()
    return any(k in s for k in _XENO_RE)


def _clean(seq: str) -> str:
    return "".join(c for c in (seq or "").upper() if c.isalpha())


@lru_cache(maxsize=1)
def _human_self_pool() -> tuple[tuple[str, str], ...]:
    """
    (label, sequence) pairs for human / near-self endogenous proteins.

    Sources: extra_protein_sequences.csv (non-xeno rows) + unique UniProt
    cache entries referenced by those rows.
    """
    pool: dict[str, str] = {}

    if _EXTRA_CSV.is_file():
        with _EXTRA_CSV.open(encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                inn = str(row.get("inn") or "")
                cls = str(row.get("drug_class") or "")
                notes = str(row.get("notes") or "")
                if _is_xeno_blob(f"{inn} {cls} {notes}"):
                    continue
                seq = _clean(str(row.get("sequence") or ""))
                if len(seq) < 20:
                    continue
                label = inn or str(row.get("uniprot") or "ref")
                # Keep longest per label; prefer mature-looking lengths
                if label not in pool or len(seq) > len(pool[label]):
                    pool[label] = seq
                uid = str(row.get("uniprot") or "").strip()
                if uid and uid not in pool:
                    fa = _UNIPROT_DIR / f"{uid}.fa"
                    if fa.is_file():
                        body = "".join(
                            ln.strip()
                            for ln in fa.read_text(encoding="utf-8").splitlines()
                            if ln.strip() and not ln.startswith(">")
                        )
                        body = _clean(body)
                        if len(body) >= 20:
                            pool[uid] = body

    # Deduplicate near-identical sequences by keeping first of each 12-mer signature
    out: list[tuple[str, str]] = []
    seen: set[str] = set()
    for lab, seq in pool.items():
        sig = seq[:40] + f"#{len(seq)}"
        if sig in seen:
            continue
        seen.add(sig)
        out.append((lab, seq))
    return tuple(out)


def best_human_homology(sequence: str) -> HomologyReport:
    """
    Best similarity to the human self pool → foreignness = 100*(1-sim).

    Uses k-mer Jaccard for screening, then SequenceMatcher.ratio on top hits.
    """
    q = _clean(sequence)
    refs = _human_self_pool()
    if len(q) < 15 or not refs:
        return HomologyReport(
            best_identity=0.0,
            best_label="",
            foreignness=50.0,
            n_refs=len(refs),
            method_note_zh="人源同源库不可用或序列过短；退回中性异源度。",
            method_note_en="Human self pool unavailable or sequence too short; neutral foreignness.",
        )

    # Length-gated candidates
    cands: list[tuple[float, str, str]] = []
    for lab, ref in refs:
        ratio = min(len(q), len(ref)) / max(len(q), len(ref))
        if ratio < 0.35:
            continue
        # Fast screen
        kj = kmer_jaccard(q, ref, k=4 if min(len(q), len(ref)) < 80 else 5)
        if kj < 0.08 and ratio < 0.85:
            continue
        cands.append((kj, lab, ref))
    cands.sort(reverse=True)

    best_sim = 0.0
    best_lab = ""
    for kj, lab, ref in cands[:24]:
        # Blend window identity-ish ratio with k-mer for short peptides
        sim = float(sequence_ratio(q, ref))
        if len(q) < 60 or len(ref) < 60:
            sim = 0.55 * sim + 0.45 * float(kj)
        else:
            sim = 0.70 * sim + 0.30 * float(kj)
        if sim > best_sim:
            best_sim = sim
            best_lab = lab

    # If nothing screened in, still try a few length-matched refs
    if best_sim <= 0 and refs:
        by_len = sorted(refs, key=lambda t: abs(len(t[1]) - len(q)))[:8]
        for lab, ref in by_len:
            sim = float(sequence_ratio(q, ref))
            if sim > best_sim:
                best_sim = sim
                best_lab = lab

    foreign = float(np.clip(100.0 * (1.0 - best_sim), 0.0, 100.0))
    note_zh = (
        f"对人源内源参考库比对：最佳相似={best_sim:.2f}（{best_lab or '无'}），"
        f"序列异源度={foreign:.0f}；库规模={len(refs)}。"
    )
    note_en = (
        f"vs human self pool: best similarity={best_sim:.2f} ({best_lab or 'none'}), "
        f"sequence foreignness={foreign:.0f}; pool n={len(refs)}."
    )
    return HomologyReport(
        best_identity=float(best_sim),
        best_label=best_lab,
        foreignness=foreign,
        n_refs=len(refs),
        method_note_zh=note_zh,
        method_note_en=note_en,
    )
