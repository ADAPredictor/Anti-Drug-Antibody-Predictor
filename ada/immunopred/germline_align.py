"""
Best-hit human V-gene germline mismatch analysis (no HMMER / hmmscan).

Aligns truncated VH/VL domains to a bundled IMGT-like human V FASTA (or ANARCI
package germlines / tiny embedded fallback) with Biopython PairwiseAligner.
CDR3 is largely excluded by truncating the query to the V-region length.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

from .cdr import AntibodyRegionMap, annotate_antibody_regions

# V-gene typically ends before full CDR3; truncate to avoid CDR3 junk mismatches.
_VH_V_MAX = 110
_VL_V_MAX = 107

_METHOD_ZH = (
    "Biopython 最佳命中对比捆绑的人源 IMGT 样 V 基因（非 ANARCI/HMMER）；"
    "仅比对 V 区约前 100–110（VH）/ 95–107（VL）aa，CDR3 大多排除或预期高错配。"
)
_METHOD_EN = (
    "Biopython best-hit vs bundled human IMGT-like V genes (not ANARCI/HMMER); "
    "aligns only V-region approx. first ~100–110 aa (VH) / ~95–107 (VL) — "
    "CDR3 largely excluded or expected high mismatch."
)

# Tiny embedded set if FASTA and ANARCI are both unavailable.
_EMBEDDED: list[tuple[str, str, str]] = [
    (
        "H",
        "IGHV1-69*01",
        "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCAR",
    ),
    (
        "H",
        "IGHV3-23*01",
        "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK",
    ),
    (
        "H",
        "IGHV3-30*01",
        "QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR",
    ),
    (
        "H",
        "IGHV4-34*01",
        "QVQLQQWGAGLLKPSETLSLTCAVYGGSFSGYYWSWIRQPPGKGLEWIGEINHSGSTNYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAR",
    ),
    (
        "H",
        "IGHV4-59*01",
        "QVQLQESGPGLVKPSETLSLTCTVSGGSISSYYWSWIRQPPGKGLEWIGYIYYSGSTNYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAR",
    ),
    (
        "K",
        "IGKV1-39*01",
        "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTP",
    ),
    (
        "K",
        "IGKV3-20*01",
        "EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSP",
    ),
    (
        "L",
        "IGLV2-14*01",
        "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTLV",
    ),
]


@dataclass
class GermlineHit:
    chain: str  # H or L (K/L → L for reporting)
    v_gene: str
    identity: float
    mismatch_abs_positions: list[int]  # 0-based on full concatenated sequence
    n_mismatch: int
    n_aligned: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "chain": self.chain,
            "v_gene": self.v_gene,
            "identity": round(float(self.identity), 4),
            "mismatch_abs_positions": list(self.mismatch_abs_positions),
            "n_mismatch": int(self.n_mismatch),
            "n_aligned": int(self.n_aligned),
        }


@dataclass
class GermlineMismatchReport:
    hits: list[GermlineHit]
    mismatch_abs: set[int] = field(default_factory=set)
    cdr_mismatch_abs: set[int] = field(default_factory=set)
    fr_mismatch_abs: set[int] = field(default_factory=set)
    hetero_frac: float = 0.0
    cdr_hetero_frac: float = 0.0
    fr_hetero_frac: float = 0.0
    method_note_zh: str = _METHOD_ZH
    method_note_en: str = _METHOD_EN

    def to_dict(self) -> dict[str, Any]:
        return {
            "hits": [h.to_dict() for h in self.hits],
            "n_mismatch": len(self.mismatch_abs),
            "n_cdr_mismatch": len(self.cdr_mismatch_abs),
            "n_fr_mismatch": len(self.fr_mismatch_abs),
            "hetero_frac": round(float(self.hetero_frac), 4),
            "cdr_hetero_frac": round(float(self.cdr_hetero_frac), 4),
            "fr_hetero_frac": round(float(self.fr_hetero_frac), 4),
            "method_note_zh": self.method_note_zh,
            "method_note_en": self.method_note_en,
        }


def _package_data_fasta() -> Path:
    # src/immunopred/ → project root data/
    here = Path(__file__).resolve()
    return here.parents[2] / "data" / "germlines" / "human_v_genes.fasta"


def _parse_fasta(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    name = ""
    chain = "?"
    parts: list[str] = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name and parts:
                    rows.append((chain, name, "".join(parts).upper()))
                header = line[1:].split()
                name = header[0] if header else "unknown"
                chain = "?"
                for tok in header[1:]:
                    if tok.startswith("chain="):
                        chain = tok.split("=", 1)[1]
                if chain == "?":
                    if name.startswith("IGHV"):
                        chain = "H"
                    elif name.startswith("IGKV"):
                        chain = "K"
                    elif name.startswith("IGLV"):
                        chain = "L"
                parts = []
            else:
                parts.append("".join(c for c in line if c.isalpha()))
        if name and parts:
            rows.append((chain, name, "".join(parts).upper()))
    return rows


def _from_anarci() -> list[tuple[str, str, str]]:
    from anarci.germlines import all_germlines

    V = all_germlines["V"]
    out: list[tuple[str, str, str]] = []
    for ch in ("H", "K", "L"):
        hum = V[ch]["human"]
        by_gene: dict[str, tuple[str, str, str]] = {}
        for name, seq in hum.items():
            gene = name.split("*")[0]
            allele = name.split("*")[1] if "*" in name else "99"
            aa = "".join(c for c in seq if c.isalpha()).upper()
            if len(aa) < 70:
                continue
            cur = by_gene.get(gene)
            if cur is None or (allele == "01" and cur[2] != "01") or (
                cur[2] != "01" and allele < cur[2]
            ):
                by_gene[gene] = (name, aa, allele)
        for gene in sorted(by_gene):
            name, aa, _ = by_gene[gene]
            out.append((ch, name, aa))
    return out


@lru_cache(maxsize=1)
def load_human_v_germlines() -> tuple[tuple[str, str, str], ...]:
    """Return ((chain, name, aa), ...) preferring bundled FASTA."""
    path = _package_data_fasta()
    if path.is_file():
        rows = _parse_fasta(path)
        if rows:
            return tuple(rows)
    try:
        rows = _from_anarci()
        if rows:
            return tuple(rows)
    except Exception:  # noqa: BLE001
        pass
    return tuple(_EMBEDDED)


def _report_chain(db_chain: str) -> str:
    return "H" if db_chain == "H" else "L"


def _v_truncate(seq: str, kind: str) -> str:
    if kind == "vh":
        return seq[:_VH_V_MAX]
    return seq[:_VL_V_MAX]


def _aligner():
    from Bio.Align import PairwiseAligner

    aln = PairwiseAligner()
    aln.mode = "global"
    aln.match_score = 2
    aln.mismatch_score = -1
    aln.open_gap_score = -5
    aln.extend_gap_score = -1
    return aln


def _identity_and_query_mismatches(
    query: str, subject: str, aligner
) -> tuple[float, list[int], int]:
    """
    Global alignment; return (identity, 0-based mismatch indices on query, n_aligned).
    Gaps in query are skipped for mismatch positions; identity = matches / aligned cols
    where both residues are letters.
    """
    if not query or not subject:
        return 0.0, [], 0
    alignments = aligner.align(query, subject)
    if not alignments:
        return 0.0, [], 0
    aln = next(iter(alignments))
    # Biopython 1.80+: aln[0], aln[1] are aligned strings with gaps
    try:
        q_aln = str(aln[0])
        s_aln = str(aln[1])
    except Exception:  # noqa: BLE001
        # older API
        q_aln, s_aln = aln
        q_aln, s_aln = str(q_aln), str(s_aln)

    matches = 0
    aligned = 0
    mism_q: list[int] = []
    qi = 0
    for a, b in zip(q_aln, s_aln):
        if a == "-":
            continue
        if b == "-":
            # insertion relative to germline — count as mismatch at query pos
            mism_q.append(qi)
            aligned += 1
            qi += 1
            continue
        aligned += 1
        if a == b:
            matches += 1
        else:
            mism_q.append(qi)
        qi += 1
    ident = (matches / aligned) if aligned else 0.0
    return float(ident), mism_q, int(aligned)


def _best_hit_for_domain(
    domain_seq: str,
    kind: str,
    offset: int,
    library: Iterable[tuple[str, str, str]],
    aligner,
) -> GermlineHit | None:
    q = _v_truncate(domain_seq, kind)
    if len(q) < 60:
        return None
    want = "H" if kind == "vh" else {"K", "L"}
    best: GermlineHit | None = None
    for ch, name, aa in library:
        if kind == "vh" and ch != "H":
            continue
        if kind == "vl" and ch not in want:
            continue
        # Subject also truncated to V length
        subj = aa[: len(q) + 5]
        ident, mism_local, n_aln = _identity_and_query_mismatches(q, subj, aligner)
        if n_aln < 50:
            continue
        abs_pos = [offset + i for i in mism_local if 0 <= i < len(domain_seq)]
        hit = GermlineHit(
            chain=_report_chain(ch),
            v_gene=name,
            identity=ident,
            mismatch_abs_positions=abs_pos,
            n_mismatch=len(abs_pos),
            n_aligned=n_aln,
        )
        if best is None or hit.identity > best.identity:
            best = hit
    return best


def analyze_germline_mismatches(
    sequence: str,
    regions: AntibodyRegionMap | None = None,
) -> GermlineMismatchReport | None:
    """
    Best-hit V-gene mismatches mapped to absolute coordinates.

    Returns None if no antibody-like domains / no usable alignment.
    """
    seq = "".join(c for c in (sequence or "").upper() if c.isalpha())
    if len(seq) < 80:
        return None
    reg = regions if regions is not None else annotate_antibody_regions(seq)
    if not reg.domains:
        return None

    library = load_human_v_germlines()
    aligner = _aligner()
    hits: list[GermlineHit] = []

    for dom in reg.domains:
        kind = dom.kind
        if kind == "unknown":
            # try VH then VL on subdomain
            sub = seq[dom.offset : dom.offset + dom.length]
            for try_kind in ("vh", "vl"):
                hit = _best_hit_for_domain(sub, try_kind, dom.offset, library, aligner)
                if hit and hit.identity >= 0.55:
                    hits.append(hit)
                    break
            continue
        if kind not in {"vh", "vl"}:
            continue
        sub = seq[dom.offset : dom.offset + dom.length]
        hit = _best_hit_for_domain(sub, kind, dom.offset, library, aligner)
        if hit:
            hits.append(hit)

    if not hits:
        return None

    mismatch_abs: set[int] = set()
    for h in hits:
        mismatch_abs.update(h.mismatch_abs_positions)

    cdr_mask = reg.cdr_mask
    cdr_mm: set[int] = set()
    fr_mm: set[int] = set()
    for i in mismatch_abs:
        if 0 <= i < len(cdr_mask) and cdr_mask[i]:
            cdr_mm.add(i)
        else:
            fr_mm.add(i)

    n_aligned = sum(h.n_aligned for h in hits)
    n_mm = sum(h.n_mismatch for h in hits)
    hetero = (n_mm / n_aligned) if n_aligned else 0.0

    # CDR / FR fractions relative to residues of that class within the aligned V window
    cdr_aligned = 0
    fr_aligned = 0
    for h, dom in zip(hits, [d for d in reg.domains if d.kind in {"vh", "vl", "unknown"}]):
        # approximate: count mask over truncated domain
        start = dom.offset
        end = start + min(dom.length, _VH_V_MAX if dom.kind == "vh" else _VL_V_MAX)
        for i in range(start, min(end, len(cdr_mask))):
            if cdr_mask[i]:
                cdr_aligned += 1
            else:
                fr_aligned += 1
    # If zip length mismatch (unknown domains), recompute simply from all hits' span
    if cdr_aligned + fr_aligned == 0:
        for i in range(len(cdr_mask)):
            # only positions that appear in any hit's aligned window approx
            pass
        for h in hits:
            for p in range(min(h.n_aligned, len(cdr_mask))):
                # weak fallback — use mismatch denominators only
                break
        cdr_aligned = max(len(cdr_mm), 1)
        fr_aligned = max(len(fr_mm), 1)

    cdr_hetero = (len(cdr_mm) / max(cdr_aligned, 1)) if cdr_aligned else 0.0
    fr_hetero = (len(fr_mm) / max(fr_aligned, 1)) if fr_aligned else 0.0

    return GermlineMismatchReport(
        hits=hits,
        mismatch_abs=mismatch_abs,
        cdr_mismatch_abs=cdr_mm,
        fr_mismatch_abs=fr_mm,
        hetero_frac=float(hetero),
        cdr_hetero_frac=float(cdr_hetero),
        fr_hetero_frac=float(fr_hetero),
        method_note_zh=_METHOD_ZH,
        method_note_en=_METHOD_EN,
    )


def epitope_overlap_mismatches(
    strong_starts: list[int],
    peptide_len: int,
    mismatch_abs: set[int],
    min_hit: int = 3,
) -> int:
    """Count peptides (1-based starts) that cover >= min_hit mismatch residues."""
    if not strong_starts or not mismatch_abs:
        return 0
    n = 0
    for st in strong_starts:
        s0 = int(st) - 1
        e0 = s0 + int(peptide_len)
        hit = sum(1 for i in range(s0, e0) if i in mismatch_abs)
        if hit >= min_hit:
            n += 1
    return n


def span_overlap_mismatches(
    spans: list[tuple[int, int]],
    mismatch_abs: set[int],
    min_hit: int = 3,
    *,
    one_based: bool = True,
) -> int:
    """Count inclusive (start, end) spans overlapping >= min_hit mismatch residues."""
    if not spans or not mismatch_abs:
        return 0
    n = 0
    for a, b in spans:
        if one_based:
            s0 = int(a) - 1
            e0 = int(b)  # exclusive in 0-based
        else:
            s0 = int(a)
            e0 = int(b) + 1
        hit = sum(1 for i in range(s0, e0) if i in mismatch_abs)
        if hit >= min_hit:
            n += 1
    return n
