"""Shared helpers for epitope scanning."""
import math

from ..config import STRONG_AFFINITY_NM, WEAK_AFFINITY_NM


def affinity_level(affinity_nm,
                   strong=STRONG_AFFINITY_NM,
                   weak=WEAK_AFFINITY_NM):
    if affinity_nm <= strong:
        return "strong"
    if affinity_nm <= weak:
        return "intermediate"
    return "weak"


def heat(affinity_nm):
    """Map affinity (nM) to a 0..1 track intensity (higher = better binder)."""
    v = (4.3 - math.log10(max(affinity_nm, 1.0))) / 2.6
    return max(0.0, min(1.0, round(v, 3)))


def merge_hits(peptide_hits, max_gap=1):
    """Merge overlapping hit peptides into epitope segments.

    peptide_hits: list of dicts with start,end,peptide,affinity,allele_breaks
    Returns sorted segments with representative peptide & allele union.
    """
    hits = sorted(peptide_hits, key=lambda h: (h["start"], h["end"]))
    segments = []
    for hit in hits:
        if segments and hit["start"] <= segments[-1]["end"] + max_gap:
            seg = segments[-1]
            seg["end"] = max(seg["end"], hit["end"])
            if hit["affinity"] < seg["affinity"]:
                seg["affinity"] = hit["affinity"]
                seg["peptide"] = hit["peptide"]
                seg["peptide_start"] = hit["start"]
            for a in hit["strong_alleles"]:
                if a not in seg["strong_alleles"]:
                    seg["strong_alleles"].append(a)
            seg["promiscuity"] = max(seg["promiscuity"], hit["promiscuity"])
            seg["n_peptides"] += 1
        else:
            segments.append({
                "start": hit["start"],
                "end": hit["end"],
                "peptide": hit["peptide"],
                "peptide_start": hit["start"],
                "affinity": hit["affinity"],
                "strong_alleles": list(hit["strong_alleles"]),
                "promiscuity": hit["promiscuity"],
                "n_peptides": 1,
            })
    for seg in segments:
        seg["level"] = affinity_level(seg["affinity"])
    return segments
