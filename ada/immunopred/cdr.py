"""
Heuristic antibody CDR / framework annotation.

Uses conserved Cys/Trp (and WGXG / FGXG) anchors — Chothia/Kabat-style —
without requiring IMGT germline databases. Good enough to separate CDR vs
framework epitope burden for early screening.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from typing import Any
try:
    from typing import Literal
except ImportError:  # Python 3.7
    from typing_extensions import Literal

ChainKind = Literal["vh", "vl", "unknown"]


@dataclass
class CDRSpan:
    name: str  # e.g. H1, L3
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    sequence: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class DomainAnnotation:
    kind: ChainKind
    offset: int  # 0-based start of this domain in the full query
    length: int
    cdrs: list[CDRSpan]
    framework_residues: int
    cdr_residues: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "offset": self.offset,
            "length": self.length,
            "cdrs": [c.to_dict() for c in self.cdrs],
            "framework_residues": self.framework_residues,
            "cdr_residues": self.cdr_residues,
        }


@dataclass
class AntibodyRegionMap:
    is_antibody_like: bool
    domains: list[DomainAnnotation]
    cdr_mask: list[bool]  # per residue, True = CDR
    cdr_fraction: float
    method_note_zh: str
    method_note_en: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "is_antibody_like": self.is_antibody_like,
            "domains": [d.to_dict() for d in self.domains],
            "cdr_fraction": round(self.cdr_fraction, 4),
            "cdr_residue_count": int(sum(self.cdr_mask)),
            "method_note_zh": self.method_note_zh,
            "method_note_en": self.method_note_en,
            "cdrs": [
                {**c.to_dict(), "domain": d.kind, "abs_start": c.start + d.offset, "abs_end": c.end + d.offset}
                for d in self.domains
                for c in d.cdrs
            ],
        }


_VH_HEAD = re.compile(r"^(EVQL|QVQL|QIQL|DVQL|AVQL|QLQL|QMQL|QITL|QVTL)")
_VL_HEAD = re.compile(
    r"^(DIQM|DIQL|DIVM|DIVL|EIVL|EIVM|AIQL|DVVMT|EIVMT|DIQMT|DIVLT|EIVLT|SSELT|QSALT|QSVLT)"
)


def _find_all(seq: str, aa: str) -> list[int]:
    return [i for i, c in enumerate(seq) if c == aa]


def _annotate_vh(seq: str) -> list[CDRSpan] | None:
    """Return H1/H2/H3 as 1-based inclusive spans, or None if anchors fail."""
    n = len(seq)
    if n < 90:
        return None
    # Conserved Cys22-ish
    cys = [i for i in _find_all(seq, "C") if 15 <= i <= 35]
    if not cys:
        return None
    c1 = cys[0]
    # Trp36-ish after Cys
    trps = [i for i in _find_all(seq, "W") if c1 + 8 <= i <= c1 + 22]
    if not trps:
        return None
    w36 = trps[0]
    # CDR-H1: roughly after Cys+4 .. before Trp
    h1_s = min(c1 + 4, w36 - 1)
    h1_e = max(h1_s, w36 - 1)
    if h1_e - h1_s < 3:
        h1_s, h1_e = max(0, w36 - 8), w36 - 1

    # CDR-H2: after FR2 (~WVRQ…), typically starts ~15 aa after W36; keep compact
    h2_s = min(n - 8, w36 + 14)
    cys2 = [i for i in _find_all(seq, "C") if w36 + 40 <= i <= n - 15]
    if not cys2:
        h2_e = min(n - 20, h2_s + 16)
        h3_s = min(n - 12, h2_e + 25)
        h3_e = n - 1
    else:
        c92 = cys2[0]
        h2_e = max(h2_s + 5, min(c92 - 10, h2_s + 16))
        m = re.search(r"WG.G", seq[c92:])
        if m:
            h3_e = c92 + m.start() - 1
        else:
            h3_e = min(n - 1, c92 + 20)
        h3_s = c92 + 2
        if h3_e < h3_s:
            h3_e = min(n - 1, h3_s + 8)
        # Cap extreme H3
        if h3_e - h3_s > 30:
            h3_e = h3_s + 30

    spans = [
        ("H1", h1_s, h1_e),
        ("H2", h2_s, h2_e),
        ("H3", h3_s, h3_e),
    ]
    out: list[CDRSpan] = []
    for name, s, e in spans:
        s = max(0, min(s, n - 1))
        e = max(s, min(e, n - 1))
        out.append(CDRSpan(name=name, start=s + 1, end=e + 1, sequence=seq[s : e + 1]))
    return out


def _annotate_vl(seq: str) -> list[CDRSpan] | None:
    n = len(seq)
    if n < 80:
        return None
    cys = [i for i in _find_all(seq, "C") if 16 <= i <= 30]
    if not cys:
        return None
    c1 = cys[0]
    trps = [i for i in _find_all(seq, "W") if c1 + 6 <= i <= c1 + 20]
    if not trps:
        # some VL use other FR2 aromatics — soft fail
        trps = [i for i in range(c1 + 10, min(n, c1 + 18))]
        if not trps:
            return None
    w35 = trps[0]
    l1_s = min(c1 + 3, w35 - 1)
    l1_e = max(l1_s, w35 - 1)

    cys2 = [i for i in _find_all(seq, "C") if w35 + 35 <= i <= n - 10]
    if not cys2:
        l2_s = min(n - 20, w35 + 14)
        l2_e = min(n - 12, l2_s + 8)
        l3_s = min(n - 10, l2_e + 20)
        l3_e = n - 1
    else:
        c88 = cys2[0]
        l2_s = min(c88 - 30, w35 + 14)
        l2_e = max(l2_s + 4, min(c88 - 20, l2_s + 9))
        m = re.search(r"FG.G", seq[c88:])
        if m:
            l3_e = c88 + m.start() - 1
        else:
            l3_e = min(n - 1, c88 + 14)
        l3_s = c88 + 1
        if l3_e < l3_s:
            l3_e = min(n - 1, l3_s + 7)

    out: list[CDRSpan] = []
    for name, s, e in [("L1", l1_s, l1_e), ("L2", l2_s, l2_e), ("L3", l3_s, l3_e)]:
        s = max(0, min(int(s), n - 1))
        e = max(s, min(int(e), n - 1))
        out.append(CDRSpan(name=name, start=s + 1, end=e + 1, sequence=seq[s : e + 1]))
    return out


def _guess_kind(seq: str) -> ChainKind:
    head = seq[:8]
    if _VH_HEAD.match(head):
        return "vh"
    if _VL_HEAD.match(head):
        return "vl"
    # weaker: Cys22 + Trp36 pattern density
    if len(seq) >= 100 and "C" in seq[15:35] and "W" in seq[30:50]:
        return "vh"
    if len(seq) >= 90 and "C" in seq[16:30]:
        return "vl"
    return "unknown"


def _split_domains(seq: str) -> list[tuple[int, str, ChainKind]]:
    """
    Split possible VH+VL (or VL+VH) concatenations into domains.
    Returns list of (offset, subdomain_seq, kind).
    """
    n = len(seq)
    if n < 160:
        kind = _guess_kind(seq)
        return [(0, seq, kind)]

    best: tuple[int, ChainKind] | None = None
    for i in range(100, min(n - 80, 160)):
        kind = _guess_kind(seq[i:])
        if kind != "unknown":
            best = (i, kind)
            if kind == "vl" and _guess_kind(seq[:i]) in {"vh", "unknown"}:
                break
            if kind == "vh" and _guess_kind(seq[:i]) == "vl":
                break
    if best is None:
        # Forced mid split if sequence looks long enough for Fv
        cut = 120 if n >= 200 else n // 2
        k1 = _guess_kind(seq[:cut]) or "vh"
        k2 = _guess_kind(seq[cut:])
        if k1 == "vh" and k2 != "vl":
            k2 = "vl"
        elif k1 == "vl" and k2 != "vh":
            k2 = "vh"
        elif k1 == "unknown":
            k1, k2 = "vh", "vl"
        return [(0, seq[:cut], k1), (cut, seq[cut:], k2)]

    cut, k2 = best
    k1 = _guess_kind(seq[:cut])
    if k1 == "unknown":
        k1 = "vh" if k2 == "vl" else "vl"
    # Prefer VH then VL ordering for typical H+L paste
    if k1 == "vh" and k2 == "vh":
        k2 = "vl"
    if k1 == "vl" and k2 == "vl":
        k2 = "vh"
    return [(0, seq[:cut], k1), (cut, seq[cut:], k2)]


def annotate_antibody_regions(sequence: str) -> AntibodyRegionMap:
    seq = "".join(c for c in sequence.upper() if c.isalpha())
    n = len(seq)
    mask = [False] * n
    domains: list[DomainAnnotation] = []

    parts = _split_domains(seq)
    for offset, sub, kind in parts:
        cdrs: list[CDRSpan] | None = None
        if kind == "vh":
            cdrs = _annotate_vh(sub)
        elif kind == "vl":
            cdrs = _annotate_vl(sub)
        else:
            # try both
            cdrs = _annotate_vh(sub) or _annotate_vl(sub)
            if cdrs and cdrs[0].name.startswith("H"):
                kind = "vh"
            elif cdrs:
                kind = "vl"

        if not cdrs:
            domains.append(
                DomainAnnotation(
                    kind=kind,
                    offset=offset,
                    length=len(sub),
                    cdrs=[],
                    framework_residues=len(sub),
                    cdr_residues=0,
                )
            )
            continue

        cdr_res = 0
        for c in cdrs:
            for i in range(c.start - 1, c.end):
                j = offset + i
                if 0 <= j < n:
                    mask[j] = True
                    cdr_res += 1
        domains.append(
            DomainAnnotation(
                kind=kind,
                offset=offset,
                length=len(sub),
                cdrs=cdrs,
                framework_residues=max(0, len(sub) - cdr_res),
                cdr_residues=cdr_res,
            )
        )

    n_cdr = sum(mask)
    is_ab = n_cdr >= 12 and any(d.cdrs for d in domains)
    return AntibodyRegionMap(
        is_antibody_like=is_ab,
        domains=domains,
        cdr_mask=mask,
        cdr_fraction=(n_cdr / n) if n else 0.0,
        method_note_zh=(
            "CDR 由保守 Cys/Trp 与 WGXG/FGXG 锚点 Heuristic 划分（类 Chothia），"
            "非 IMGT 正式编号；用于早期筛选区分 CDR 与框架表位负担。"
        ),
        method_note_en=(
            "CDRs are inferred from conserved Cys/Trp and WGXG/FGXG anchors "
            "(Chothia-like), not official IMGT numbering; used to separate CDR vs "
            "framework epitope burden for early screening."
        ),
    )


def peptide_overlaps_cdr(
    start_1based: int,
    length: int,
    cdr_mask: list[bool],
    min_cdr_residues: int = 5,
) -> bool:
    """True if enough peptide residues lie in CDR (avoids FR/CDR border inflation)."""
    s = start_1based - 1
    e = s + length
    n = len(cdr_mask)
    hit = 0
    for i in range(max(0, s), min(n, e)):
        if cdr_mask[i]:
            hit += 1
            if hit >= min_cdr_residues:
                return True
    return False


def region_epitope_burden(
    strong_starts: list[int],
    peptide_length: int,
    cdr_mask: list[bool],
) -> dict[str, float]:
    """
    Split strong MHC binder starts into CDR-overlapping vs framework-only.

    Densities for reporting: per 100 aa of the respective region.
    non_germline_burden: CDR hits full weight + framework hits ×0.15, per 100 aa
    of the full sequence (early-screen primary metric).
    """
    n = max(len(cdr_mask), 1)
    cdr_len = max(sum(cdr_mask), 1)
    fw_len = max(n - sum(cdr_mask), 1)
    cdr_hits = 0
    fw_hits = 0
    for st in strong_starts:
        if peptide_overlaps_cdr(st, peptide_length, cdr_mask):
            cdr_hits += 1
        else:
            fw_hits += 1
    per100 = n / 100.0
    return {
        "cdr_strong_count": float(cdr_hits),
        "framework_strong_count": float(fw_hits),
        "cdr_epitope_density": (cdr_hits / cdr_len) * 100.0,
        "framework_epitope_density": (fw_hits / fw_len) * 100.0,
        "non_germline_burden": (cdr_hits * 1.0 + fw_hits * 0.15) / max(per100, 0.01),
    }
