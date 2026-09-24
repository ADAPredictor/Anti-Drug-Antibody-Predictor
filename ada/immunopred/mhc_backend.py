"""
Pluggable MHC prediction backends: heuristic (default), IEDB, local NetMHC.

On any remote/local failure, auto-fallback to the built-in heuristic predictors.
Never raises network/subprocess failures to callers.
"""

from __future__ import annotations

import csv
import io
import os
import shutil
import subprocess
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

import numpy as np

from .mhc_i import (
    DEFAULT_MHC1_ALLELES,
    MHC1Epitope,
    MHC1Result,
    predict_mhc1_epitopes,
)
from .sequence import ProteinSequence
from .t_cell import (
    DEFAULT_ALLELES,
    TCellEpitope,
    TCellResult,
    predict_t_cell_epitopes,
)

VALID_BACKENDS = frozenset({"heuristic", "iedb", "local_netmhc"})
_ALIASES = {
    "local": "local_netmhc",
    "netmhc": "local_netmhc",
    "netmhcpan": "local_netmhc",
}

IEDB_MHCII_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
IEDB_MHCI_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
IEDB_TIMEOUT_S = 45
IEDB_MHCII_METHOD = "netmhciipan_el-4.3"
IEDB_MHCI_METHOD = "netmhcpan_el"

_NOTES = {
    "heuristic": {
        "zh": "使用内置 MHC Heuristic（口袋偏好矩阵）。",
        "en": "Using built-in MHC heuristic (pocket-preference matrices).",
    },
    "iedb": {
        "zh": "使用 IEDB NetMHCpan / NetMHCIIpan EL。",
        "en": "Using IEDB NetMHCpan / NetMHCIIpan EL.",
    },
    "local_netmhc": {
        "zh": "使用本地 NetMHC 可执行文件。",
        "en": "Using local NetMHC executable.",
    },
    "fallback": {
        "zh": "已回退到内置 Heuristic。",
        "en": "Fell back to built-in heuristic.",
    },
}


def resolve_backend(requested: str | None = None) -> str:
    """Resolve backend name from argument or IMMUNOPRED_MHC_BACKEND (default heuristic)."""
    raw = (requested if requested is not None else os.environ.get("IMMUNOPRED_MHC_BACKEND")) or "heuristic"
    key = str(raw).strip().lower().replace("-", "_")
    key = _ALIASES.get(key, key)
    if key not in VALID_BACKENDS:
        return "heuristic"
    return key


def _meta(
    requested: str,
    used: str,
    fallback_reason: str = "",
    note_zh: str | None = None,
    note_en: str | None = None,
) -> dict[str, str]:
    base = _NOTES.get(used, _NOTES["heuristic"])
    if fallback_reason:
        note_zh = note_zh or f"{_NOTES['fallback']['zh']}（{fallback_reason}）"
        note_en = note_en or f"{_NOTES['fallback']['en']} ({fallback_reason})"
    else:
        note_zh = note_zh or base["zh"]
        note_en = note_en or base["en"]
    return {
        "backend_requested": requested,
        "backend_used": used,
        "fallback_reason": fallback_reason or "",
        "note_zh": note_zh,
        "note_en": note_en,
    }


def _public_meta(meta: dict[str, str]) -> dict[str, str]:
    """Shape expected by PredictResult.to_dict() / API."""
    return {
        "requested": meta.get("backend_requested", ""),
        "used": meta.get("backend_used", ""),
        "fallback_reason": meta.get("fallback_reason", ""),
        "note_zh": meta.get("note_zh", ""),
        "note_en": meta.get("note_en", ""),
    }


def meta_for_result(meta: dict[str, str]) -> dict[str, str]:
    return _public_meta(meta)


def _normalize_mhc2_allele(allele: str) -> str:
    a = allele.strip()
    if a.upper().startswith("HLA-"):
        return a
    if a.upper().startswith("DRB") or a.upper().startswith("DQA") or a.upper().startswith("DP"):
        return f"HLA-{a}"
    return a


def _normalize_mhc1_allele(allele: str) -> str:
    a = allele.strip()
    if a.upper().startswith("HLA-"):
        return a
    if a.upper().startswith(("A*", "B*", "C*")):
        return f"HLA-{a}"
    return a


def _iedb_post(url: str, fields: dict[str, str], timeout: float = IEDB_TIMEOUT_S) -> str:
    data = urllib.parse.urlencode(fields).encode("utf-8")
    req = urllib.request.Request(url, data=data, method="POST")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def _parse_tsv(text: str) -> list[dict[str, str]]:
    text = text.strip()
    if not text:
        return []
    # IEDB sometimes prefixes errors; require a header line with tabs
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or "\t" not in lines[0]:
        raise ValueError(f"Unexpected IEDB response: {text[:200]}")
    reader = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
    return [dict(row) for row in reader if row]


def _rank_to_percentile(rank_val: float) -> float:
    """IEDB rank/%rank: lower = stronger. Map to percentile where higher = stronger."""
    return float(max(0.0, min(100.0, 100.0 - rank_val)))


def _is_strong(rank_val: float, percentile: float) -> bool:
    return rank_val <= 10.0 or percentile >= 90.0


def _summarize_alleles(
    alleles: list[str],
    hits: list[Any],
    score_attr: str = "score",
) -> dict[str, dict[str, float]]:
    summary: dict[str, dict[str, float]] = {}
    for allele in alleles:
        allele_hits = [h for h in hits if h.allele == allele or h.allele.endswith(allele.replace("HLA-", ""))]
        if not allele_hits:
            # also match by normalized display
            allele_hits = [h for h in hits if allele in h.allele or h.allele in allele]
        strong = sum(1 for h in allele_hits if h.strong)
        scores = [float(getattr(h, score_attr)) for h in allele_hits]
        summary[allele] = {
            "n_peptides": float(len(allele_hits)),
            "strong": float(strong),
            "mean_score": float(np.mean(scores)) if scores else 0.0,
            "max_score": float(np.max(scores)) if scores else 0.0,
        }
    return summary


def _heuristic_mhc2(
    protein: ProteinSequence,
    alleles: list[str] | None = None,
    top_n: int = 30,
    **kwargs: Any,
) -> TCellResult:
    return predict_t_cell_epitopes(protein, alleles=alleles, top_n=top_n, **kwargs)


def _heuristic_mhc1(
    protein: ProteinSequence,
    alleles: list[str] | None = None,
    top_n: int = 30,
    **kwargs: Any,
) -> MHC1Result:
    return predict_mhc1_epitopes(protein, alleles=alleles, top_n=top_n, **kwargs)


def _iedb_mhc2(
    protein: ProteinSequence,
    alleles: list[str] | None = None,
    top_n: int = 30,
) -> TCellResult:
    alleles = alleles or DEFAULT_ALLELES
    iedb_alleles = [_normalize_mhc2_allele(a) for a in alleles]
    # Prefer one request with comma-separated alleles; fall back per-allele.
    rows: list[dict[str, str]] = []
    try:
        raw = _iedb_post(
            IEDB_MHCII_URL,
            {
                "method": IEDB_MHCII_METHOD,
                "sequence_text": protein.sequence,
                "allele": ",".join(iedb_alleles),
            },
        )
        rows = _parse_tsv(raw)
    except (urllib.error.URLError, TimeoutError, OSError, ValueError):
        rows = []
        for iedb_allele in iedb_alleles:
            raw = _iedb_post(
                IEDB_MHCII_URL,
                {
                    "method": IEDB_MHCII_METHOD,
                    "sequence_text": protein.sequence,
                    "allele": iedb_allele,
                },
            )
            rows.extend(_parse_tsv(raw))

    allele_alias = { _normalize_mhc2_allele(a): a for a in alleles }
    all_hits: list[TCellEpitope] = []
    for row in rows:
        try:
            start = int(float(row.get("start") or 0))
            peptide = str(row.get("peptide") or "")
            core = str(row.get("core_peptide") or row.get("core") or "")
            score = float(row.get("score") or 0.0)
            rank_raw = row.get("rank") or row.get("percentile_rank") or row.get("%rank") or "100"
            rank_val = float(rank_raw)
            row_allele = str(row.get("allele") or "")
        except (TypeError, ValueError):
            continue
        if not peptide:
            continue
        display = allele_alias.get(row_allele) or allele_alias.get(_normalize_mhc2_allele(row_allele))
        if display is None:
            for orig, norm in ((a, _normalize_mhc2_allele(a)) for a in alleles):
                if norm == row_allele or orig in row_allele or row_allele.endswith(orig.replace("HLA-", "")):
                    display = orig
                    break
        display = display or row_allele
        pct = _rank_to_percentile(rank_val)
        all_hits.append(
            TCellEpitope(
                start=start,
                peptide=peptide,
                core=core or (peptide[3:12] if len(peptide) >= 12 else peptide[:9]),
                allele=display,
                score=score,
                percentile=pct,
                strong=_is_strong(rank_val, pct),
            )
        )

    strong_epitopes = [e for e in all_hits if e.strong]
    strong_epitopes.sort(key=lambda e: (e.percentile, e.score), reverse=True)
    unique_strong_starts = {e.start for e in strong_epitopes}
    density = (len(unique_strong_starts) / max(protein.length, 1)) * 100.0
    top = strong_epitopes[:top_n]
    mean_top = float(np.mean([e.score for e in top])) if top else 0.0
    return TCellResult(
        epitopes=top,
        strong_count=len(unique_strong_starts),
        epitope_density=density,
        mean_top_score=mean_top,
        allele_summary=_summarize_alleles(alleles, all_hits),
        strong_starts=sorted(unique_strong_starts),
    )


def _iedb_mhc1(
    protein: ProteinSequence,
    alleles: list[str] | None = None,
    top_n: int = 30,
) -> MHC1Result:
    alleles = alleles or DEFAULT_MHC1_ALLELES
    iedb_alleles = [_normalize_mhc1_allele(a) for a in alleles]
    rows: list[dict[str, str]] = []
    try:
        raw = _iedb_post(
            IEDB_MHCI_URL,
            {
                "method": IEDB_MHCI_METHOD,
                "sequence_text": protein.sequence,
                "allele": ",".join(iedb_alleles),
                "length": ",".join(["9"] * len(iedb_alleles)),
            },
        )
        rows = _parse_tsv(raw)
    except (urllib.error.URLError, TimeoutError, OSError, ValueError):
        rows = []
        for iedb_allele in iedb_alleles:
            raw = _iedb_post(
                IEDB_MHCI_URL,
                {
                    "method": IEDB_MHCI_METHOD,
                    "sequence_text": protein.sequence,
                    "allele": iedb_allele,
                    "length": "9",
                },
            )
            rows.extend(_parse_tsv(raw))

    allele_alias = { _normalize_mhc1_allele(a): a for a in alleles }
    all_hits: list[MHC1Epitope] = []
    for row in rows:
        try:
            start = int(float(row.get("start") or 0))
            peptide = str(row.get("peptide") or "")
            score = float(row.get("score") or 0.0)
            rank_raw = (
                row.get("percentile_rank")
                or row.get("rank")
                or row.get("%rank")
                or "100"
            )
            rank_val = float(rank_raw)
            row_allele = str(row.get("allele") or "")
        except (TypeError, ValueError):
            continue
        if not peptide:
            continue
        display = allele_alias.get(row_allele) or allele_alias.get(_normalize_mhc1_allele(row_allele))
        if display is None:
            for orig in alleles:
                if _normalize_mhc1_allele(orig) == row_allele or orig == row_allele:
                    display = orig
                    break
        display = display or row_allele
        pct = _rank_to_percentile(rank_val)
        all_hits.append(
            MHC1Epitope(
                start=start,
                peptide=peptide,
                allele=display,
                score=score,
                percentile=pct,
                strong=_is_strong(rank_val, pct),
            )
        )

    strong_epitopes = [e for e in all_hits if e.strong]
    strong_epitopes.sort(key=lambda e: (e.percentile, e.score), reverse=True)
    unique_starts = {e.start for e in strong_epitopes}
    density = (len(unique_starts) / max(protein.length, 1)) * 100.0
    top = strong_epitopes[:top_n]
    mean_top = float(np.mean([e.score for e in top])) if top else 0.0
    return MHC1Result(
        epitopes=top,
        strong_count=len(unique_starts),
        epitope_density=density,
        mean_top_score=mean_top,
        allele_summary=_summarize_alleles(alleles, all_hits),
    )


def _find_netmhc_binary(kind: str) -> str | None:
    """kind: 'mhc2' or 'mhc1'."""
    if kind == "mhc2":
        env_keys = ("IMMUNOPRED_NETMHCIIPAN", "NETMHCIIPAN")
        names = ("netMHCIIpan", "netMHCIIpan.exe", "netmhciipan")
    else:
        env_keys = ("IMMUNOPRED_NETMHCPAN", "NETMHCPAN")
        names = ("netMHCpan", "netMHCpan.exe", "netmhcpan")
    for key in env_keys:
        path = os.environ.get(key)
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return path
        if path and os.path.isfile(path):
            return path  # Windows may lack X_OK semantics
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None


def _local_netmhc_mhc2(
    protein: ProteinSequence,
    alleles: list[str] | None = None,
    top_n: int = 30,
) -> TCellResult:
    binary = _find_netmhc_binary("mhc2")
    if not binary:
        raise RuntimeError("local NetMHC binary not configured")
    # Minimal stub: invoke once; parsing is site-specific — if output unusable, raise.
    alleles = alleles or DEFAULT_ALLELES
    allele = _normalize_mhc2_allele(alleles[0])
    proc = subprocess.run(
        [binary, "-a", allele, "-f", "-"],
        input=f">query\n{protein.sequence}\n",
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"netMHCIIpan exited {proc.returncode}: {(proc.stderr or '')[:200]}")
    # Without a stable parse path, treat as not configured for production use.
    raise RuntimeError("local NetMHC output parsing not configured")


def _local_netmhc_mhc1(
    protein: ProteinSequence,
    alleles: list[str] | None = None,
    top_n: int = 30,
) -> MHC1Result:
    binary = _find_netmhc_binary("mhc1")
    if not binary:
        raise RuntimeError("local NetMHC binary not configured")
    alleles = alleles or DEFAULT_MHC1_ALLELES
    allele = _normalize_mhc1_allele(alleles[0])
    proc = subprocess.run(
        [binary, "-a", allele, "-l", "9", "-f", "-"],
        input=f">query\n{protein.sequence}\n",
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"netMHCpan exited {proc.returncode}: {(proc.stderr or '')[:200]}")
    raise RuntimeError("local NetMHC output parsing not configured")


def predict_mhc2(
    protein: ProteinSequence,
    backend: str | None = None,
    alleles: list[str] | None = None,
    top_n: int = 30,
    **kwargs: Any,
) -> tuple[TCellResult, dict[str, str]]:
    requested = resolve_backend(backend)
    if requested == "heuristic":
        result = _heuristic_mhc2(protein, alleles=alleles, top_n=top_n, **kwargs)
        return result, _meta(requested, "heuristic")

    try:
        if requested == "iedb":
            result = _iedb_mhc2(protein, alleles=alleles, top_n=top_n)
            return result, _meta(requested, "iedb")
        if requested == "local_netmhc":
            result = _local_netmhc_mhc2(protein, alleles=alleles, top_n=top_n)
            return result, _meta(requested, "local_netmhc")
    except (urllib.error.URLError, TimeoutError, OSError, ValueError, RuntimeError, KeyError) as exc:
        reason = str(exc) or exc.__class__.__name__
        result = _heuristic_mhc2(protein, alleles=alleles, top_n=top_n, **kwargs)
        return result, _meta(requested, "heuristic", fallback_reason=reason)
    except Exception as exc:  # noqa: BLE001 — never raise to caller
        reason = str(exc) or exc.__class__.__name__
        result = _heuristic_mhc2(protein, alleles=alleles, top_n=top_n, **kwargs)
        return result, _meta(requested, "heuristic", fallback_reason=reason)

    result = _heuristic_mhc2(protein, alleles=alleles, top_n=top_n, **kwargs)
    return result, _meta(requested, "heuristic")


def predict_mhc1(
    protein: ProteinSequence,
    backend: str | None = None,
    alleles: list[str] | None = None,
    top_n: int = 30,
    **kwargs: Any,
) -> tuple[MHC1Result, dict[str, str]]:
    requested = resolve_backend(backend)
    if requested == "heuristic":
        result = _heuristic_mhc1(protein, alleles=alleles, top_n=top_n, **kwargs)
        return result, _meta(requested, "heuristic")

    try:
        if requested == "iedb":
            result = _iedb_mhc1(protein, alleles=alleles, top_n=top_n)
            return result, _meta(requested, "iedb")
        if requested == "local_netmhc":
            result = _local_netmhc_mhc1(protein, alleles=alleles, top_n=top_n)
            return result, _meta(requested, "local_netmhc")
    except (urllib.error.URLError, TimeoutError, OSError, ValueError, RuntimeError, KeyError) as exc:
        reason = str(exc) or exc.__class__.__name__
        result = _heuristic_mhc1(protein, alleles=alleles, top_n=top_n, **kwargs)
        return result, _meta(requested, "heuristic", fallback_reason=reason)
    except Exception as exc:  # noqa: BLE001
        reason = str(exc) or exc.__class__.__name__
        result = _heuristic_mhc1(protein, alleles=alleles, top_n=top_n, **kwargs)
        return result, _meta(requested, "heuristic", fallback_reason=reason)

    result = _heuristic_mhc1(protein, alleles=alleles, top_n=top_n, **kwargs)
    return result, _meta(requested, "heuristic")
