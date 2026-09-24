"""Protein sequence input / FASTA parsing / validation."""
import re
from collections import namedtuple

Chain = namedtuple("Chain", ["chain_id", "name", "seq"])

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
ALLOWED_AA = CANONICAL_AA | {"X"}  # X = unknown / mixture

_FASTA_HEADER = re.compile(r"^>(.*)$")


class SequenceError(ValueError):
    pass


def clean_sequence(raw):
    """Uppercase, strip whitespace / numbering / FASTA gap characters."""
    seq = re.sub(r"\s+", "", raw).upper()
    seq = seq.replace(".", "").replace("-", "").replace("*", "")
    return seq


def parse_fasta(text):
    """Parse FASTA text into chains.

    Tolerant of two common input styles:
    * a headerless blob is treated as a single chain;
    * a ``>SEQUENCE`` line whose "name" is itself bare amino acids (a very
      common user slip) starts a new chain instead of swallowing the
      sequence as a header and leaving the record empty.

    Returns list[Chain].
    """
    text = text.strip()
    if not text:
        raise SequenceError("提供的序列为空。")

    if not text.startswith(">"):
        # headerless input: blank-line-separated blocks are separate chains
        # (e.g. heavy / light chains pasted one after another); a single
        # block stays a single chain even if soft-wrapped
        blocks = [b for b in re.split(r"\n\s*\n", text) if b.strip()]
        letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        chains = []
        for i, block in enumerate(blocks):
            seq = clean_sequence(block)
            _validate(seq)
            cid = letters[i] if i < 26 else "Z"
            chains.append(Chain(cid, "query", seq))
        return chains

    chains = []
    header = None
    buf = []
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    def _chain_id():
        return letters[len(chains)] if len(chains) < 26 else "Z"

    def _flush():
        seq = clean_sequence("".join(buf))
        del buf[:]
        if not seq:
            return
        if header is None:
            chains.append(Chain(_chain_id(), "query", seq))
        else:
            parts = header.split(None, 1)
            cid = parts[0] if parts and parts[0] else _chain_id()
            name = parts[1] if len(parts) > 1 else cid
            chains.append(Chain(cid, name, seq))

    for line in text.splitlines():
        line = line.strip()
        m = _FASTA_HEADER.match(line)
        if not m:
            buf.append(line)
            continue
        content = m.group(1).strip()
        _flush()
        if content and _looks_like_sequence(content):
            header = None
            buf.append(content)
        else:
            header = content or None
    _flush()

    if not chains:
        raise SequenceError("未能从 FASTA 中解析出序列。")
    for chain in chains:
        _validate(chain.seq)
    return chains


def _looks_like_sequence(text):
    """True when a '>' line is a bare AA string (>= one class-I peptide),
    i.e. the user wrote ``>SEQUENCE`` instead of ``>name``."""
    return len(text) >= MIN_CHAIN_LENGTH and all(
        aa in ALLOWED_AA for aa in text)


MIN_CHAIN_LENGTH = 8  # shortest class-I MHC peptide; below this epitope scanning is empty


def _validate(seq):
    if len(seq) < MIN_CHAIN_LENGTH:
        raise SequenceError(
            "序列过短（%d 个残基，少于 %d），至少需要一条 %d 残基以上的链才能进行表位分析。"
            % (len(seq), MIN_CHAIN_LENGTH, MIN_CHAIN_LENGTH)
        )
    invalid = sorted({aa for aa in seq if aa not in ALLOWED_AA})
    if invalid:
        raise SequenceError(
            "序列中存在非法氨基酸字符: %s（仅允许 20 种标准氨基酸及 X）" % ", ".join(invalid)
        )


def total_length(chains):
    return sum(len(c.seq) for c in chains)
