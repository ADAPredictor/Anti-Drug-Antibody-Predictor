"""T-cell epitope prediction with a pluggable MHC-binding backend.

Backends
--------
* ``HeuristicClass1`` - always available, zero-dependency screening baseline.
  Allele supertype anchor motifs (P2 / C-terminal preferences, SYFPEITHI-
  style) are scored and mapped to a pseudo-affinity in nM.  Intended for
  triage, NOT as a replacement for pan-allele predictors.
* ``MHCflurryClass1`` - used automatically when the free, Apache-licensed
  ``mhcflurry`` package and its models are installed.  Implements real
  pan-allele MHC-I binding prediction; alleles unsupported by the local
  models transparently fall back to the heuristic.

MHC-II scanning uses a pocket-hydrophobicity baseline (P1/P4/P6/P9) and
exposes the same interface so a NetMHCIIPan adapter can be added later.
"""
import atexit
import json
import os
import re
import subprocess
import threading
from pathlib import Path

from ..config import STRONG_AFFINITY_NM, WEAK_AFFINITY_NM
from .base import affinity_level, heat, merge_hits

_WORKER_PATH = Path(__file__).with_name("mhcflurry_worker.py")

# --------------------------------------------------------------------------
# MHC-I heuristic
# --------------------------------------------------------------------------
# supertype -> (matching numbers, preferred P2, preferred C-terminal residue)
SUPERTYPE_TABLE = {
    "A01": (("01",), "TS", "YF"),
    "A02": (("02",), "LMIVT", "VLIA"),
    "A03": (("03", "11", "30", "31", "32", "33"), "LVI", "KRY"),
    "A24": (("23", "24"), "YFW", "FLIWM"),
    "B07": (("07", "35", "51", "53", "54"), "P", "VLIFYA"),
    "B08": (("08",), "KR", "LF"),
    "B44": (("18", "40", "44", "45"), "E", "YFLW"),
}
_GENERIC_P2 = "LVIMA"
_GENERIC_PO = "LVIMFWY"
_ALLELE_RE = re.compile(r"HLA-?([ABC])\*?(\d{1,2})")

_FAVORABLE = set("AVLIMFYW")
_ACIDIC = set("DE")


def parse_supertype(allele):
    m = _ALLELE_RE.search(allele.upper().replace(":", "*"))
    if not m:
        return "generic"
    gene, num = m.group(1), m.group(2)
    for name, (nums, _p2, _po) in SUPERTYPE_TABLE.items():
        if name[0] == gene and num in nums:
            return name
    return "generic"


def _anchor_score(residue, preferred):
    if residue in preferred:
        return 1.0
    return 0.0


def heuristic_class1_affinity(peptide, allele):
    """Pseudo-affinity (nM) from anchor + interior preferences."""
    supertype = parse_supertype(allele)
    if supertype == "generic":
        p2, po = _GENERIC_P2, _GENERIC_PO
    else:
        _, p2, po = SUPERTYPE_TABLE[supertype]
    s_p2 = _anchor_score(peptide[1], p2)
    s_po = _anchor_score(peptide[-1], po)
    interior = peptide[2:-1]
    fav = sum(1.0 for aa in interior if aa in _FAVORABLE) / len(interior)
    pen = sum(1.0 for aa in interior if aa in _ACIDIC) / len(interior)
    s_int = max(0.0, fav - 0.5 * pen)
    score = 0.42 * s_p2 + 0.42 * s_po + 0.16 * s_int
    return 10.0 ** (4.4 - 2.8 * score)


# --------------------------------------------------------------------------
# MHCflurry adapter (optional, isolated Python >= 3.10 via subprocess)
# --------------------------------------------------------------------------
def _find_mhcflurry_python():
    """Resolve the interpreter that has mhcflurry installed.

    Override with ADA_MHCFLURRY_PYTHON; otherwise use the embeddable runtime
    provisioned by tools/setup_mhcflurry.ps1.
    """
    env_path = os.environ.get("ADA_MHCFLURRY_PYTHON")
    if env_path and Path(env_path).exists():
        return env_path
    if os.environ.get("ADA_DISABLE_MHCFLURRY"):
        return None
    candidate = Path(__file__).resolve().parents[2] / "tools" / "py310" / "python.exe"
    return str(candidate) if candidate.exists() else None


class MHCflurrySubprocess:
    """Persistent worker: the mhcflurry models stay loaded in one child
    process and serve many scans over a JSON-lines channel (thread-safe)."""

    name = "mhcflurry"

    def __init__(self, python_exe, startup_timeout=300):
        self.python_exe = python_exe
        self._lock = threading.Lock()
        self.proc = subprocess.Popen(
            [python_exe, str(_WORKER_PATH), "--serve"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True, bufsize=1, encoding="utf-8",
        )
        try:
            # the worker answers its first probe line only after loading the
            # models; that first round-trip is the one-time warm-up cost
            data = self._roundtrip({"probe": True}, startup_timeout)
            self.version = data.get("version", "unknown")
        except Exception:
            self.close()
            raise
        atexit.register(self.close)

    def _kill_worker(self):
        proc = self.proc
        if proc is not None and proc.poll() is None:
            try:
                proc.kill()
            except Exception:
                pass

    def _roundtrip(self, payload, timeout):
        line = json.dumps(payload, ensure_ascii=True) + "\n"
        with self._lock:
            if self.proc is None or self.proc.poll() is not None:
                raise RuntimeError("mhcflurry worker is not running")
            self.proc.stdin.write(line)
            self.proc.stdin.flush()
            watchdog = threading.Timer(timeout, self._kill_worker)
            watchdog.start()
            try:
                raw = self.proc.stdout.readline()
            finally:
                watchdog.cancel()
        if not raw:
            raise RuntimeError("mhcflurry worker closed or timed out")
        data = json.loads(raw)
        if not data.get("ok"):
            raise RuntimeError(data.get("error", "worker reported failure"))
        return data

    def close(self):
        with self._lock:
            if self.proc is not None and self.proc.poll() is None:
                try:
                    self.proc.stdin.close()
                    self.proc.wait(timeout=5)
                except Exception:
                    try:
                        self.proc.kill()
                    except Exception:
                        pass
            self.proc = None

    def predict_groups(self, groups, alleles, timeout=900):
        """Flatten every length group into one peptide list per allele, send
        a single request, and slice the returned columns back per group.
        Returns ``{tag: {allele: nM-list-or-None}}``."""
        flat_peptides, slices = [], []
        for tag, peptides in groups:
            start = len(flat_peptides)
            flat_peptides.extend(peptides)
            slices.append((tag, start, len(flat_peptides)))
        requests = [{"allele": a, "peptides": flat_peptides}
                    for a in alleles]
        data = self._roundtrip({"requests": requests}, timeout)
        out = dict((tag, {}) for tag, _, _ in slices)
        for allele, col in zip(alleles, data["results"]):
            for tag, start, end in slices:
                out[tag][allele] = (col[start:end] if col is not None
                                    else None)
        return out


def _try_load_mhcflurry():
    python_exe = _find_mhcflurry_python()
    if not python_exe:
        return None
    try:
        return MHCflurrySubprocess(python_exe)
    except Exception:
        return None


_BACKEND = None


def get_class1_backend():
    """Singleton: mhcflurry worker when usable, else False (heuristic)."""
    global _BACKEND
    if _BACKEND is None:
        _BACKEND = _try_load_mhcflurry() or False
    return _BACKEND


def disable_class1_backend():
    """Retire a dead/timed-out persistent worker (stays heuristic until
    process restart)."""
    global _BACKEND
    if _BACKEND:
        try:
            _BACKEND.close()
        except Exception:
            pass
    _BACKEND = False


# --------------------------------------------------------------------------
# Scanning
# --------------------------------------------------------------------------
def _affinity_matrix(peptides, alleles, mhcflurry_cols):
    """Return ({allele: (values_nM, source_name)}, n_from_mhcflurry)."""
    columns = {}
    n_mhcflurry = 0
    for allele in alleles:
        col = (mhcflurry_cols or {}).get(allele)
        if col is not None and len(col) == len(peptides):
            columns[allele] = (col, "mhcflurry")
            n_mhcflurry += 1
        else:
            columns[allele] = (
                [heuristic_class1_affinity(p, allele) for p in peptides],
                "heuristic",
            )
    return columns, n_mhcflurry


def scan_class1(seq, alleles, lengths=None,
                strong=STRONG_AFFINITY_NM, weak=WEAK_AFFINITY_NM):
    backend = get_class1_backend()
    if lengths is None:
        lengths = (8, 9, 10, 11) if backend else (9,)

    groups = []
    for length in lengths:
        if len(seq) >= length:
            groups.append(
                (length, [seq[i:i + length]
                          for i in range(len(seq) - length + 1)]))

    # one persistent worker round-trip for the whole scan; if the worker is
    # dead/timed out it is retired and this scan falls back to the heuristic
    predicted = {}
    if backend and groups:
        try:
            predicted = backend.predict_groups(groups, alleles)
        except Exception:
            predicted = {}
            disable_class1_backend()
            backend = None
    backend_name = ("mhcflurry-%s-subprocess" % backend.version
                    if backend else "heuristic-v1")

    track = [0.0] * len(seq)
    peptide_hits = []
    n_scanned = 0
    n_windows = 0
    best = None
    any_fallback = False

    for length, peptides in groups:
        matrix, n_mhcflurry = _affinity_matrix(
            peptides, alleles, predicted.get(length, {}))
        if backend and n_mhcflurry < len(alleles):
            any_fallback = True
        n_scanned += len(peptides) * len(alleles)
        n_windows += len(peptides)

        values_by_allele = {a: col[0] for a, col in matrix.items()}
        for i, pep in enumerate(peptides):
            strong_alleles, n_bind = [], 0
            row_best = None
            for allele in alleles:
                aff = values_by_allele[allele][i]
                if best is None or aff < best:
                    best = aff
                if row_best is None or aff < row_best:
                    row_best = aff
                if aff <= weak:
                    n_bind += 1
                if aff <= strong:
                    strong_alleles.append(allele)
            for pos in range(i, i + length):
                h = heat(row_best)
                if h > track[pos]:
                    track[pos] = h
            # burden counting only uses strong binders; weak promiscuity is
            # retained for display but never promotes a peptide to a hit
            if strong_alleles:
                peptide_hits.append({
                    "start": i,
                    "end": i + length,
                    "peptide": pep,
                    "affinity": row_best,
                    "strong_alleles": strong_alleles,
                    "promiscuity": n_bind,
                })

    promiscuous = sum(1 for h in peptide_hits if len(h["strong_alleles"]) >= 2)
    used = "mhcflurry+heuristic-fallback" if any_fallback else backend_name
    return {
        "backend": used,
        "lengths": list(lengths),
        "peptides_scanned": n_scanned,
        "windows": n_windows,
        "strong_hit_count": len(peptide_hits),
        "alleles": list(alleles),
        "track": track,
        "segments": merge_hits(peptide_hits),
        "summary": {
            "strong_hit_peptides": len(peptide_hits),
            "promiscuous_strong_peptides": promiscuous,
            "best_affinity_nm": round(best, 1) if best is not None else None,
        },
    }


# --------------------------------------------------------------------------
# MHC-II heuristic
# --------------------------------------------------------------------------
_II_POCKETS = (0, 3, 5, 8)  # P1, P4, P6, P9 of the 9-aa core inside a 15-mer
_II_HYDROPHOBIC = set("AVLIMFYW")


def heuristic_class2_affinity(core9):
    # steeper curve: all 4 pockets hydrophobic -> ~40 nM (strong),
    # 3/4 -> ~500 nM (intermediate), 2/4 -> ~6300 nM (weak)
    matches = sum(1 for p in _II_POCKETS if core9[p] in _II_HYDROPHOBIC)
    score = matches / len(_II_POCKETS)
    return 10.0 ** (5.2 - 3.6 * score)


def scan_class2(seq, alleles=None, core_len=15,
                strong=STRONG_AFFINITY_NM, weak=WEAK_AFFINITY_NM):
    if len(seq) < core_len:
        return _empty_class2(seq, alleles)
    track = [0.0] * len(seq)
    peptide_hits = []
    best = None
    for i in range(len(seq) - core_len + 1):
        pep = seq[i:i + core_len]
        aff = heuristic_class2_affinity(pep)
        best = aff if best is None else min(best, aff)
        for pos in range(i, i + core_len):
            h = heat(aff)
            if h > track[pos]:
                track[pos] = h
        if aff <= strong:
            peptide_hits.append({
                "start": i,
                "end": i + core_len,
                "peptide": pep,
                "affinity": aff,
                "strong_alleles": ["generic MHC-II"],
                "promiscuity": 1,
            })
    n_windows = len(seq) - core_len + 1
    return {
        "backend": "heuristic-classII-v1",
        "lengths": [core_len],
        "peptides_scanned": n_windows,
        "windows": n_windows,
        "strong_hit_count": len(peptide_hits),
        "alleles": list(alleles or []),
        "track": track,
        "segments": merge_hits(peptide_hits),
        "summary": {
            "strong_hit_peptides": len(peptide_hits),
            "promiscuous_strong_peptides": 0,
            "best_affinity_nm": round(best, 1) if best is not None else None,
        },
    }


def _empty_class2(seq, alleles):
    return {
        "backend": "heuristic-classII-v1",
        "lengths": [15],
        "peptides_scanned": 0,
        "windows": 0,
        "strong_hit_count": 0,
        "alleles": list(alleles or []),
        "track": [0.0] * len(seq),
        "segments": [],
        "summary": {"strong_hit_peptides": 0, "promiscuous_strong_peptides": 0,
                    "best_affinity_nm": None},
    }
