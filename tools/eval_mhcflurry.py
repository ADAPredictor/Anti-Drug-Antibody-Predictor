# -*- coding: utf-8 -*-
"""Evaluate the real mhcflurry MHC-I backend against the heuristic baseline.

Reproducible evaluation (NOT an ADA-incidence validation):

    python tools/eval_mhcflurry.py

Runs every benchmark protein twice with use_cache=False - first pinned to
the heuristic anchor-motif baseline, then through the mhcflurry 2.2.1
subprocess models (tools/py310 + tools/mhcflurry_data) - and reports the
MHC-I strong-binder rates and total scores side by side.  Also asserts that
the canonical influenza epitope GILGFVFTL is predicted as a strong
HLA-A*02:01 binder (<=500 nM) by the neural model.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ada import db
from ada.epitopes import tcell
from ada.etl import load_seed_data
from ada.pipeline import analyze

BENCHMARK = [
    "insulin glargine", "glucagon", "teriparatide", "exenatide",
    "pramlintide", "bivalirudin", "filgrastim", "epoetin alfa",
    "interferon alfa-2b", "somatropin",
]


def run_pass(seq, modality, real_backend):
    if real_backend:
        if tcell._BACKEND is None or tcell._BACKEND is False:
            tcell._BACKEND = None  # (re-)probe the persistent worker once
    else:
        tcell._BACKEND = False  # pin heuristic baseline
    return analyze(seq, modality=modality, use_cache=False, run_class2=True)


def summarize(r):
    cov = r["score"]["coverage"]
    mhc1 = r["epitopes"]["mhc1"]
    return {
        "score": r["score"]["total"],
        "rate": cov["mhc1_strong_rate"],
        "pep": cov["mhc1_strong_peptides"],
        "best": mhc1.get("best_affinity_nm"),
        "backend": mhc1["backend"],
    }


def main():
    db.init_db()
    load_seed_data()
    drugs = dict((d["generic"], d) for d in db.all_drugs())

    print("== heuristic pass ==")
    heur = {}
    for name in BENCHMARK:
        d = drugs[name]
        heur[name] = summarize(run_pass(d["sequence"], d["modality"], False))

    print("== mhcflurry pass (one persistent worker, models loaded once) ==")
    real = {}
    for name in BENCHMARK:
        d = drugs[name]
        r = run_pass(d["sequence"], d["modality"], True)
        real[name] = summarize(r)
        print("  done: %s [%s]" % (name, real[name]["backend"]))

    print()
    print("%-22s %-27s %-27s %s" % ("drug", "heur: score/rate/pep/best",
                                    "mhcflurry: score/rate/pep/best", "dScore"))
    for name in BENCHMARK:
        h, m = heur[name], real[name]
        print("%-22s %5.1f %5.3f %4d %8.0f   %5.1f %5.3f %4d %8.0f   %+5.1f"
              % (name, h["score"], h["rate"], h["pep"], h["best"] or -1,
                 m["score"], m["rate"], m["pep"], m["best"] or -1,
                 m["score"] - h["score"]))

    # canonical positive control: influenza M1 58-66, HLA-A*02:01
    ctrl = analyze("GILGFVFTL", modality="other",
                   class1_alleles=["HLA-A*02:01"],
                   run_class2=False, use_cache=False)
    best = ctrl["epitopes"]["mhc1"]["best_affinity_nm"]
    backend = ctrl["epitopes"]["mhc1"]["backend"]
    print()
    print("GILGFVFTL / HLA-A*02:01 -> %.1f nM [%s]" % (best, backend))
    if "mhcflurry" not in backend:
        print("WARNING: mhcflurry backend not active; cannot validate")
        return 2
    assert best is not None and best <= 500.0, \
        "canonical epitope not predicted as strong binder: %s" % best
    print("PASS: canonical influenza epitope <= 500 nM")
    return 0


if __name__ == "__main__":
    sys.exit(main())
