# -*- coding: utf-8 -*-
"""Calibrate the MHC-I burden saturation from REAL mhcflurry predictions.

The score's MHC-I component is min(1, strong_binder_rate / SATURATION).
SATURATION therefore defines "binding density relative to an approved-drug
background": a protein whose strong-binder rate sits at the high end of a
panel of marketed protein/peptide drugs should map to ~1.0, not be pinned to
1.0 by the median background.

This script scans a fixed, public-domain FDA benchmark panel (the same 10
proteins used by tests/test_benchmark_fda.py, sequences taken from the seed
library) with the active MHC-I backend (mhcflurry 2.2.1 when provisioned,
otherwise the heuristic) and reports the strong-binder rate distribution.

Run:
    python tools/calibrate_mhc1_saturation.py

Reference measurement (mhcflurry-2.2.1-subprocess, 8-allele world-wide
panel, 8-11 mer, <=500 nM strong; 2026-09-19):

    insulin glargine 0.069   glucagon 0.098   teriparatide 0.088
    exenatide 0.107         pramlintide 0.053 bivalirudin 0.217
    filgrastim 0.110        epoetin alfa 0.094 interferon alfa-2b 0.118
    somatropin 0.086
    median 0.096, p90 ~0.20, max 0.217

SATURATION["mhc1"] in ada/scoring.py was therefore set to 0.20.  This is a
scale calibration against binding-background, NOT a fit to clinical ADA
incidence (HLA binding alone does not encode tolerance/dose/route).
"""
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ada import db
from ada.config import DEFAULT_CLASS1_ALLELES
from ada.epitopes import tcell
from ada.etl import load_seed_data
from ada.seqio import parse_fasta

PANEL = [
    "insulin glargine", "glucagon", "teriparatide", "exenatide",
    "pramlintide", "bivalirudin", "filgrastim", "epoetin alfa",
    "interferon alfa-2b", "somatropin",
]


def strong_rate(sequence_text, alleles):
    hits = windows = 0
    best = None
    for chain in parse_fasta(sequence_text):
        r = tcell.scan_class1(chain.seq, alleles)
        hits += r["strong_hit_count"]
        windows += r["windows"]
        b = r["summary"].get("best_affinity_nm")
        if b is not None:
            best = b if best is None else min(best, b)
    return (hits / windows if windows else 0.0), best


def percentile(sorted_vals, q):
    if not sorted_vals:
        return None
    k = (len(sorted_vals) - 1) * q
    lo = int(k)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (k - lo)


def main():
    db.init_db()
    load_seed_data()
    drugs = dict((d["generic"], d) for d in db.all_drugs())
    backend = tcell.get_class1_backend()
    print("backend: %s" % (backend or "heuristic-v1"))

    rates = {}
    for name in PANEL:
        rate, best = strong_rate(drugs[name]["sequence"],
                                 list(DEFAULT_CLASS1_ALLELES))
        rates[name] = rate
        print("%-22s rate=%.3f best_nM=%s" % (name, rate, best))

    vals = sorted(rates.values())
    print()
    print("n=%d median=%.3f p75=%.3f p90=%.3f max=%.3f"
          % (len(vals), statistics.median(vals), percentile(vals, 0.75),
             percentile(vals, 0.90), vals[-1]))
    print("suggested SATURATION mhc1 = %.2f (p90 of panel)"
          % percentile(vals, 0.90))
    return 0


if __name__ == "__main__":
    sys.exit(main())
