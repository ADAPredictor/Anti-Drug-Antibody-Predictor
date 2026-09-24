"""Sequence parsing and validation utilities."""

from __future__ import annotations

import re
from dataclasses import dataclass
from difflib import SequenceMatcher
from io import StringIO
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.Seq import Seq

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS = set("BXZJUO*")


@dataclass(frozen=True)
class ProteinSequence:
    """A named amino-acid sequence."""

    id: str
    sequence: str
    description: str = ""
    chains: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        seq = self.sequence.upper().replace(" ", "").replace("\n", "")
        object.__setattr__(self, "sequence", seq)
        bad = sorted({aa for aa in seq if aa not in VALID_AA})
        if bad:
            raise ValueError(
                f"Sequence '{self.id}' contains non-standard residues: {''.join(bad)}"
            )

    @property
    def length(self) -> int:
        return len(self.sequence)

    def peptides(self, length: int, step: int = 1) -> list[tuple[int, str]]:
        """Return (1-based start, peptide) sliding windows."""
        if length <= 0:
            raise ValueError("peptide length must be positive")
        if self.length < length:
            return []
        return [
            (i + 1, self.sequence[i : i + length])
            for i in range(0, self.length - length + 1, step)
        ]


def clean_sequence(raw: str) -> str:
    seq = re.sub(r"[^A-Za-z]", "", raw.upper())
    if any(aa in AMBIGUOUS for aa in seq):
        raise ValueError("Ambiguous or non-standard amino acids are not supported")
    bad = sorted({aa for aa in seq if aa not in VALID_AA})
    if bad:
        raise ValueError(f"Invalid amino acids: {''.join(bad)}")
    if not seq:
        raise ValueError("Empty amino-acid sequence")
    return seq


def load_fasta(path: str | Path) -> list[ProteinSequence]:
    records: list[ProteinSequence] = []
    for rec in SeqIO.parse(str(path), "fasta"):
        records.append(
            ProteinSequence(
                id=rec.id,
                sequence=str(rec.seq),
                description=rec.description,
            )
        )
    if not records:
        raise ValueError(f"No sequences found in {path}")
    return records


def from_string(sequence: str, seq_id: str = "query") -> ProteinSequence:
    return ProteinSequence(id=seq_id, sequence=clean_sequence(sequence))


def parse_query(text: str, seq_id: str = "query") -> ProteinSequence:
    """
    Accept raw amino acids or multi-FASTA (e.g. heavy + light).

    Multi-FASTA chains are concatenated in file order (recommend heavy then light).
    Header lines starting with '>' are ignored for residue content.
    """
    raw = (text or "").strip()
    if not raw:
        raise ValueError("Empty sequence")

    looks_fasta = raw.startswith(">") or "\n>" in raw
    if looks_fasta:
        records = list(SeqIO.parse(StringIO(raw), "fasta"))
        if not records:
            raise ValueError("FASTA contains no sequences")
        chains = [clean_sequence(str(r.seq)) for r in records]
        joined = "".join(chains)
        auto_id = "+".join(r.id for r in records[:4])
        return ProteinSequence(
            id=seq_id if seq_id and seq_id != "query" else auto_id,
            sequence=joined,
            description=f"{len(chains)} chain(s) concatenated",
            chains=tuple(chains),
        )

    lines = [
        ln.strip()
        for ln in raw.splitlines()
        if ln.strip() and not ln.strip().startswith(">")
    ]
    return ProteinSequence(id=seq_id, sequence=clean_sequence("".join(lines)))


def identity(a: str, b: str) -> float:
    """Global identity after length-aware pairwise comparison."""
    if not a or not b:
        return 0.0
    short, long = (a, b) if len(a) <= len(b) else (b, a)
    n = len(short)
    if n == 0:
        return 0.0
    best = 0
    for i in range(len(long) - n + 1):
        window = long[i : i + n]
        matches = sum(x == y for x, y in zip(short, window))
        if matches > best:
            best = matches
    return best / n


def sequence_ratio(a: str, b: str) -> float:
    """Identity-like ratio. Best for variable domains of similar length."""
    if not a or not b:
        return 0.0
    if a == b:
        return 1.0
    return SequenceMatcher(None, a, b).ratio()


def kmer_jaccard(a: str, b: str, k: int = 5) -> float:
    """Jaccard similarity of k-mer sets."""
    if len(a) < k or len(b) < k:
        return identity(a, b)

    def kmers(s: str) -> set[str]:
        return {s[i : i + k] for i in range(len(s) - k + 1)}

    sa, sb = kmers(a), kmers(b)
    inter = len(sa & sb)
    union = len(sa | sb)
    return inter / union if union else 0.0


def concatenate_chains(chains: Iterable[str], linker: str = "") -> str:
    return linker.join(clean_sequence(c) for c in chains)


def to_biopython(seq: ProteinSequence) -> Seq:
    return Seq(seq.sequence)
