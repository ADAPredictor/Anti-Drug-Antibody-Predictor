# -*- coding: utf-8 -*-
"""Calibrate continuous B-cell epitope thresholds against real 3D data.

Labels: residue exposed (RSA > 0.25) from Shrake-Rupley SASA computed on
whole crystallographic assemblies / AlphaFold models.  Candidate rules
(Parker hydrophilicity, Chou-Fasman turns) are scored for precision,
enrichment over base rate, informedness and coverage on pooled residues.

Run from the repository root:
    python tools/calibrate_bcell_rsa.py

Structures are fetched once and cached under ./data/structures.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ada import structure
from ada.epitopes.bcell import parker_profile, cf_turn_profile

TARGETS = [
    ("insulin", "pdb", "1MSO"),
    ("hGH", "pdb", "1HGU"),
    ("IFNa2", "pdb", "1RH2"),
    ("GCSF", "af", "P09919"),
    ("EPO", "af", "P01588"),
]

pool = []  # (parker, cf_turn, exposed_flag)
for name, kind, acc in TARGETS:
    try:
        if kind == "pdb":
            text, src = structure.fetch_pdb_text(pdb_id=acc)
        else:
            text, src = structure.fetch_pdb_text(uniprot_acc=acc)
        residues = structure.parse_pdb(text)
        chain_order = []
        for r in residues:
            if r["chain"] not in chain_order:
                chain_order.append(r["chain"])
        # SASA over the ENTIRE assembly so crystal/biological-unit neighbours
        # are included for every residue.
        sasa = structure._shrake_rupley(residues)
        by_chain = {cid: [] for cid in chain_order}
        for r, sv in zip(residues, sasa):
            by_chain[r["chain"]].append((r, sv))
        before = len(pool)
        for cid in chain_order:
            items = by_chain[cid]
            seq = "".join(r["aa"] for r, _ in items)
            if len(seq) < 15:
                continue
            pk = parker_profile(seq)
            cf = cf_turn_profile(seq)
            for i, (r, sv) in enumerate(items):
                ref = structure.MAX_ASA.get(r["aa"], 200.0)
                exposed = 1 if (sv / ref) > structure.EXPOSED_CUTOFF else 0
                pool.append((pk[i], cf[i], exposed))
        n_added = len(pool) - before
        base = sum(l for _, _, l in pool[before:]) / n_added
        print("OK %-8s %-24s chains=%d n=%-4d exposed=%.2f" % (
            name, src, len(chain_order), n_added, base))
    except Exception as exc:
        print("FAIL %-8s %s: %r" % (name, acc, exc))

n = len(pool)
base = sum(l for _, _, l in pool) / n
print("\npooled residues: %d, base exposed rate: %.3f\n" % (n, base))
print("%-28s %6s %6s %6s %8s %8s" % ("rule", "recall", "prec", "enrich",
                                     "informed", "coverage"))


def report(label, flags):
    tp = sum(1 for f, (_, _, l) in zip(flags, pool) if f and l)
    fp = sum(1 for f, (_, _, l) in zip(flags, pool) if f and not l)
    pos = sum(l for _, _, l in pool)
    neg = n - pos
    pred = tp + fp
    rec = tp / pos if pos else 0.0
    prec = tp / pred if pred else 0.0
    informed = (tp / pos if pos else 0.0) - (fp / neg if neg else 0.0)
    print("%-28s %6.3f %6.3f %6.2f %8.3f %8.3f" % (
        label, rec, prec, prec / base if base else 0.0, informed, pred / n))


for t in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5):
    report("parker>=%.1f" % t, [p >= t for p, _, _ in pool])
for t in (1.0, 1.5, 2.0, 2.5):
    report("parker>=%.1f & CF>=1.0" % t,
           [p >= t and c >= 1.0 for p, c, _ in pool])
report("CF>=1.0", [c >= 1.0 for _, c, _ in pool])
report("parker>=1.5 | CF>=1.05", [p >= 1.5 or c >= 1.05 for p, c, _ in pool])
