"""
Benchmark: ADA internal MHC scoring vs IEDB NetMHCpan / NetMHCIIpan.

Runs both predictors on the same real antibody variable-region chains
(heavy + light) pulled from data/ada_reference.csv, then compares:
  1) per-peptide ranking agreement (Spearman rho) by chain / allele
  2) top-decile epitope overlap (Jaccard)
  3) strong-binder counts under each method's threshold
  4) ADA composite sequence-risk vs clinical ADA across the panel
     and lineage-stratified group means

All results are real API outputs; raw responses cached to paper/_iedb_cache.json.

Usage:  python paper/benchmark_iedb.py
"""
from __future__ import annotations

import json
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from ada.immunopred.sequence import ProteinSequence
from ada.immunopred.mhc_i import _mhc1_matrices, _score_peptide as score_mhc1_peptide, PEPTIDE_LENGTH as LEN1
from ada.immunopred.t_cell import (
    load_hla_matrices, _best_core, PEPTIDE_LENGTH as LEN2,
)

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(__file__).parent / "_iedb_cache.json"

PANEL = {
    # inn: known lineage
    "infliximab": "chimeric",
    "rituximab": "chimeric",
    "cetuximab": "chimeric",
    "trastuzumab": "humanized",
    "bevacizumab": "humanized",
    "omalizumab": "humanized",
    "adalimumab": "human",
    "panitumumab": "human",
    "golimumab": "human",
}

MHC1_ALLELES = ["HLA-A*02:01", "HLA-B*07:02", "HLA-A*01:01"]
MHC2_ALLELES = ["HLA-DRB1*01:01", "HLA-DRB1*04:01"]
# ADA internal allele names (strip HLA- prefix for MHC-II already plain)
ADA_MHC1 = MHC1_ALLELES
ADA_MHC2 = [a.replace("HLA-DRB1", "DRB1") for a in MHC2_ALLELES]

MHC1_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
MHC2_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"


def load_cache() -> dict:
    if CACHE.exists():
        return json.loads(CACHE.read_text(encoding="utf-8"))
    return {}


def save_cache(c: dict) -> None:
    CACHE.write_text(json.dumps(c, ensure_ascii=False, indent=1), encoding="utf-8")


def iedb_call(url: str, params: dict, cache: dict, key: str) -> str:
    if key in cache:
        return cache[key]
    data = urllib.parse.urlencode(params).encode()
    last = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data)
            txt = urllib.request.urlopen(req, timeout=120).read().decode()
            cache[key] = txt
            save_cache(cache)
            time.sleep(1.0)
            return txt
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(5 + attempt * 5)
    raise RuntimeError(f"IEDB call failed: {last}")


def parse_tsv(txt: str) -> list[dict]:
    lines = [ln for ln in txt.strip().splitlines() if ln.strip()]
    header = [h.strip() for h in lines[0].split("\t")]
    out = []
    for ln in lines[1:]:
        d = dict(zip(header, ln.split("\t")))
        # normalise column names across MHC-I / MHC-II endpoints
        if "percentile_rank" not in d and "rank" in d:
            d["percentile_rank"] = d["rank"]
        out.append(d)
    return out


def rankdata(xs: list[float]) -> list[float]:
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> float:
    if len(x) < 4:
        return float("nan")
    rx, ry = rankdata(x), rankdata(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = sum((a - mx) ** 2 for a in rx) ** 0.5
    vy = sum((b - my) ** 2 for b in ry) ** 0.5
    return cov / (vx * vy) if vx and vy else float("nan")


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a or b) else float("nan")


def _engine_percentiles(scores: list[float]) -> list[float]:
    """Replicate immunopred in-protein percentile: linspace(0,100,n) by ascending score order."""
    n = len(scores)
    order = sorted(range(n), key=lambda i: scores[i])
    ranks = [0.0] * n
    for r, idx in enumerate(order):
        ranks[idx] = (r / (n - 1) * 100.0) if n > 1 else 100.0
    return ranks


def ada_mhc1_windows(prot: ProteinSequence) -> dict[str, dict]:
    """All 9-mer windows scored by ADA internal pocket matrices (per allele)."""
    mats = _mhc1_matrices()
    windows = prot.peptides(LEN1, step=1)
    out: dict[str, dict] = {}
    for allele in ADA_MHC1:
        pssm = mats[allele]
        scored = [score_mhc1_peptide(pep, pssm) for _, pep in windows]
        pcts = _engine_percentiles(scored)
        out[allele] = {
            (start, pep): (s, bool(p >= 90.0 and s >= 0.15))
            for (start, pep), s, p in zip(windows, scored, pcts)
        }
    return out


def ada_mhc2_windows(prot: ProteinSequence) -> dict[str, dict]:
    """All 15-mer windows scored by ADA internal MHC-II PSSM (per allele)."""
    data = load_hla_matrices(None)
    allele_mats = data["alleles"]
    core_offset = int(data.get("core_offset", 3))
    windows = prot.peptides(LEN2, step=1)
    out: dict[str, dict] = {}
    for allele in ADA_MHC2:
        pssm = allele_mats[allele]["pssm"]
        thr = float(allele_mats[allele].get("threshold_strong", 1.8))
        scored = [_best_core(pep, pssm, core_offset=core_offset)[1] for _, pep in windows]
        pcts = _engine_percentiles(scored)
        out[allele] = {
            (start, pep): (s, bool(p >= 90.0 and s >= thr * 0.5))
            for (start, pep), s, p in zip(windows, scored, pcts)
        }
    return out


def compare_allele(ada_map: dict, iedb_rows: list[dict]) -> dict | None:
    common = []
    for ir in iedb_rows:
        key = (int(ir["start"]), ir["peptide"])
        if key in ada_map and ir.get("percentile_rank") not in (None, ""):
            common.append((ada_map[key][0], float(ir["percentile_rank"]), ada_map[key][1]))
    if len(common) < 10:
        return None
    x = [c[0] for c in common]
    y = [-c[1] for c in common]  # lower IEDB percentile rank = stronger
    rho = spearman(x, y)
    n10 = max(1, len(common) // 10)
    top_ada = {k for k, _ in sorted(ada_map.items(), key=lambda kv: kv[1][0], reverse=True)[:n10]}
    top_iedb = {(int(r["start"]), r["peptide"])
                for r in sorted(iedb_rows, key=lambda r: float(r["percentile_rank"]))[:n10]}
    return {
        "n": len(common),
        "spearman": round(rho, 3),
        "top10_jaccard": round(jaccard(top_ada, top_iedb), 3),
        "strong_ada": sum(1 for c in common if c[2]),
        "strong_iedb_pct2": sum(1 for c in common if c[1] <= 2.0),
    }


def main() -> None:
    df = pd.read_csv(ROOT / "data" / "ada_reference.csv")
    df["inn_l"] = df["inn"].astype(str).str.lower()
    rows = []
    for inn, lineage in PANEL.items():
        sub = df[df["inn_l"] == inn]
        if sub.empty:
            print(f"[skip] {inn} not in reference CSV")
            continue
        r = sub.iloc[0]
        if not isinstance(r["heavy"], str) or not isinstance(r["light"], str):
            print(f"[skip] {inn} missing chains")
            continue
        rows.append({"inn": inn, "lineage": lineage, "name": r["name"],
                     "heavy": r["heavy"].strip(), "light": r["light"].strip(),
                     "ada_pct": float(r["ada_rate_pct"]) if pd.notna(r["ada_rate_pct"]) else None})
    print(f"Panel size: {len(rows)} mAbs")

    cache = load_cache()
    pairwise: list[dict] = []

    for drug in rows:
        for chain_key, chain_seq in (("H", drug["heavy"]), ("L", drug["light"])):
            prot = ProteinSequence(id=f"{drug['inn']}_{chain_key}", sequence=chain_seq)

            # ADA internal all-window scores
            a1 = ada_mhc1_windows(prot)
            a2 = ada_mhc2_windows(prot)

            # ---- MHC-I vs NetMHCpan (one allele per call; IEDB requirement) ----
            for allele in MHC1_ALLELES:
                txt = iedb_call(
                    MHC1_URL,
                    {"method": "netmhcpan", "sequence_text": chain_seq,
                     "allele": allele, "length": "9"},
                    cache, f"mhc1|{drug['inn']}|{chain_key}|{allele}",
                )
                iedb1 = parse_tsv(txt)
                stats = compare_allele(a1[allele], iedb1)
                if stats:
                    pairwise.append({"drug": drug["inn"], "chain": chain_key, "mhc": "I",
                                     "allele": allele, **stats})

            # ---- MHC-II vs NetMHCIIpan ----
            txt = iedb_call(
                MHC2_URL,
                {"method": "netmhciipan", "sequence_text": chain_seq,
                 "allele": ",".join(MHC2_ALLELES)},
                cache, f"mhc2|{drug['inn']}|{chain_key}",
            )
            iedb2 = parse_tsv(txt)
            for iedb_allele, ada_allele in zip(MHC2_ALLELES, ADA_MHC2):
                stats = compare_allele(a2[ada_allele],
                                       [r for r in iedb2 if r.get("allele") == iedb_allele])
                if stats:
                    pairwise.append({"drug": drug["inn"], "chain": chain_key, "mhc": "II",
                                     "allele": iedb_allele, **stats})
        print(f"[done] {drug['inn']} ({drug['lineage']})", flush=True)

    comp = pd.DataFrame(pairwise)
    out = Path(__file__).parent / "benchmark_peptide_comparison.csv"
    comp.to_csv(out, index=False, encoding="utf-8-sig")

    print("\n=== Per-peptide ranking agreement (ADA vs IEDB) ===")
    print(comp.groupby("mhc")[["spearman", "top10_jaccard"]].agg(["mean", "median", "std"]).round(3))
    print("\nStrong-binder counts (sum across chains/alleles):")
    print(comp.groupby("mhc")[["strong_ada", "strong_iedb_pct2"]].sum())
    print(f"\nDetailed table saved: {out}")

    # ---- composite sequence-risk vs clinical ADA ----
    print("\n=== Composite sequence risk vs clinical ADA ===")
    comp_rows = []
    for drug in rows:
        if drug["ada_pct"] is None:
            continue
        fasta = f">{drug['inn']} heavy\n{drug['heavy']}\n>{drug['inn']} light\n{drug['light']}\n"
        body = json.dumps({"sequence": fasta, "name": drug["inn"],
                           "drug_class": "antibody"}).encode()
        req = urllib.request.Request(
            "http://127.0.0.1:8000/api/predict", data=body,
            headers={"Content-Type": "application/json"})
        pred = json.loads(urllib.request.urlopen(req, timeout=120).read())
        comp_rows.append({
            "inn": drug["inn"], "lineage": drug["lineage"],
            "ada_seq_risk": round(float(pred["risk"]["score"]), 1),
            "tier": pred["risk"]["tier"],
            "clinical_ada_pct": drug["ada_pct"],
            "nab_risk": round(float(pred["nab"]["score"]), 1),
        })
        time.sleep(0.5)
    cdf = pd.DataFrame(comp_rows)
    cout = Path(__file__).parent / "benchmark_composite.csv"
    cdf.to_csv(cout, index=False, encoding="utf-8-sig")
    print(cdf.to_string(index=False))
    if len(cdf) >= 5:
        rho = spearman(cdf["ada_seq_risk"].tolist(), cdf["clinical_ada_pct"].tolist())
        print(f"\nSpearman rho (ADA seq-risk vs clinical ADA), n={len(cdf)}: {rho:.3f}")
        print("\nLineage-stratified means:")
        print(cdf.groupby("lineage")[["ada_seq_risk", "clinical_ada_pct"]]
              .mean().round(1).reindex(["chimeric", "humanized", "human"]))
    print(f"Composite table saved: {cout}")


if __name__ == "__main__":
    main()
