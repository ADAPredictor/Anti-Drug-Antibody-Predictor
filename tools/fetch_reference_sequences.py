# -*- coding: utf-8 -*-
"""Fetch real reference sequences for seed drugs that lack them.

Sources of truth (no hand-typed sequences):
* UniProt REST   - wild-type human enzymes (ERT drugs); the mature chain is
                   sliced using the entry's own "Chain" feature annotation.
* RCSB PDB       - therapeutic-antibody Fab structures; sequences come from
                   experimentally determined SEQRES records and the PDB ID is
                   stored for provenance.

Usage:
    python tools/fetch_reference_sequences.py           # dry run, prints plan
    python tools/fetch_reference_sequences.py --apply   # writes seed_drugs.json

Sequences that cannot be sourced authoritatively are left null on purpose.
"""
import json
import re
import sys

import requests

import os
HERE = os.path.dirname(os.path.abspath(__file__))
SEED_PATH = os.path.join(HERE, "..", "ada", "data", "seed_drugs.json")

# generic -> authoritative source spec.
# PDB entries below were verified in-session: each entry title/FASTA header
# names the exact therapeutic (e.g. "OKT3 Fab", "MEDI-493 (Synagis)"), and
# mutant constructs / crystallization aids were explicitly rejected.
PLAN = {
    # ERT enzymes: UniProt entry, mature chain sliced per annotation
    "laronidase":         {"uniprot": "P35475"},
    "agalsidase beta":    {"uniprot": "P06280"},
    "alglucosidase alfa": {"uniprot": "P10253"},
    # antibodies: verified PDB entry + keyword that must appear in the entry
    # title or FASTA headers; chains selected per spec
    "muromonab-CD3":      {"pdb_id": "1SY6", "expect": "OKT3",
                           "chain_keyword": "OKT3"},
    "infliximab":         {"pdb_id": "5VH3", "expect": "Infliximab"},
    "rituximab":          {"pdb_id": "4KAQ", "expect": "Rituximab"},
    "bevacizumab":        {"pdb_id": "7V5N", "expect": "bevacizumab"},
    "omalizumab":         {"pdb_id": "2XA8", "expect": "omalizumab"},
    "palivizumab":        {"pdb_id": "2HWZ", "expect": "MEDI-493"},
    "certolizumab pegol": {"pdb_id": "5WUV", "expect": "Certolizumab"},
    "nivolumab":          {"pdb_id": "5GGQ", "expect": "nivolumab"},
    "natalizumab":        {"pdb_id": "4IRZ", "expect": "Natalizumab"},
    # secukinumab is only deposited as an Fv; pick its own chains and avoid
    # the XAB5 crystallization-aid Fab also present in the entry
    "secukinumab":        {"pdb_id": "9SFX", "expect": "Secukinumab",
                           "chain_keyword": "Secukinumab"},
    # NOTE: etanercept / abatacept stay null on purpose - no PDB entry of the
    # full fusion construct is verified, and a UniProt ectodomain + Fc
    # assembly would require asserting an unverified linker sequence.
}


def fetch_uniprot_chain(acc):
    """(mature sequence, chain name, annotated?) from UniProt JSON."""
    r = requests.get("https://rest.uniprot.org/uniprotkb/%s.json" % acc,
                     timeout=30)
    r.raise_for_status()
    entry = r.json()
    seq = entry["sequence"]["value"]
    for feat in entry.get("features", []):
        if feat.get("type") == "Chain":
            loc = feat["location"]
            start = loc["start"]["value"]
            end = loc["end"]["value"]
            return (seq[start - 1:end], feat.get("description") or acc, True)
    return seq, acc, False  # no Chain annotation: full precursor


def fetch_pdb_fasta(entry_id):
    """[(header, seq)] from RCSB FASTA, one per entity."""
    r = requests.get("https://www.rcsb.org/fasta/entry/%s" % entry_id,
                     timeout=30)
    r.raise_for_status()
    out = []
    header, buf = None, []
    for line in r.text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                out.append((header, "".join(buf)))
            header, buf = line[1:], []
        elif line:
            buf.append(line)
    if header is not None:
        out.append((header, "".join(buf)))
    return out


def fab_pair(entry_id, chain_keyword=None):
    """Two Fab chains (200-280 aa) of an entry; with a keyword, pick the two
    chains whose headers contain it and relax the length window (Fv)."""
    chains = fetch_pdb_fasta(entry_id)
    if chain_keyword:
        picked = [(h, s) for h, s in chains
                  if chain_keyword.lower() in h.lower()]
        lo, hi = 100, 280
    else:
        picked = chains
        lo, hi = 200, 280
    picked = [(h, s) for h, s in picked if lo <= len(s) <= hi]
    return picked if len(picked) == 2 else None


def _entry_title(entry_id):
    r = requests.get(
        "https://data.rcsb.org/rest/v1/core/entry/%s" % entry_id, timeout=30)
    r.raise_for_status()
    return r.json().get("struct", {}).get("title", "")


def fill_pdb(rec, generic):
    """Verified-entry lookup: expect keyword in entry title or FASTA headers."""
    spec = PLAN[generic]
    entry_id = spec["pdb_id"]
    expect = spec["expect"].lower()
    title = _entry_title(entry_id)
    headers = [h for h, _ in fetch_pdb_fasta(entry_id)]
    if not (expect in title.lower()
            or any(expect in h.lower() for h in headers)):
        raise RuntimeError(
            "entry %s does not name %r (title: %s)"
            % (entry_id, expect, title[:80]))
    pair = fab_pair(entry_id, spec.get("chain_keyword"))
    if not pair:
        raise RuntimeError("entry %s has no usable 2-chain pair" % entry_id)
    fasta = "\n".join(
        ">%s\n%s" % (_chain_label(h, entry_id, i), s)
        for i, (h, s) in enumerate(pair)) + "\n"
    return fasta, entry_id


def _chain_label(header, entry_id, idx):
    """Human-readable chain name from an RCSB FASTA header."""
    parts = header.split("|")
    for p in parts[1:]:
        if "chain" in p.lower():
            return re.sub(r"[^A-Za-z0-9_]", "_", p.strip()).strip("_")
    return "%s_entity%d" % (entry_id, idx + 1)


def main():
    apply_mode = "--apply" in sys.argv
    with open(SEED_PATH, "r", encoding="utf-8") as fh:
        records = json.load(fh)
    by_generic = dict((r["generic"], r) for r in records)

    filled, skipped = [], []
    for generic, spec in PLAN.items():
        rec = by_generic[generic]
        if rec.get("sequence"):
            continue  # already sourced

        if "uniprot" in spec:
            acc = spec["uniprot"]
            try:
                seq, name, annotated = fetch_uniprot_chain(acc)
            except Exception as exc:
                skipped.append((generic, "uniprot %s: %s" % (acc, exc)))
                continue
            if not annotated:
                skipped.append((generic,
                                "uniprot %s: no Chain feature" % acc))
                continue
            print("[uniprot %s] %-18s %d aa  mature=%s"
                  % (acc, generic, len(seq), name))
            if apply_mode:
                rec["sequence"] = ">%s\n%s\n" % (
                    name.replace(" ", "_").lower(), seq)
                rec["uniprot_acc"] = acc
            filled.append(generic)
            continue

        try:
            got = fill_pdb(rec, generic)
        except Exception as exc:
            got = None
            print("[pdb] %-18s error: %s" % (generic, exc))
        if not got:
            skipped.append((generic, "no 2-chain Fab entry found on PDB"))
            continue
        fasta, entry_id = got
        lens = [len(s) for s in
                [seg for seg in fasta.split("\n") if seg and not
                 seg.startswith(">")]]
        print("[pdb %s] %-18s chains=%s" % (entry_id, generic, lens))
        if apply_mode:
            rec["sequence"] = fasta
            rec["pdb_id"] = entry_id
        filled.append(generic)

    print("\nfilled: %d, skipped: %d" % (len(filled), len(skipped)))
    for g, why in skipped:
        print("  SKIP %-18s %s" % (g, why))
    if apply_mode:
        with open(SEED_PATH, "w", encoding="utf-8") as fh:
            json.dump(records, fh, ensure_ascii=False, indent=2)
            fh.write("\n")
        print("seed_drugs.json updated")


if __name__ == "__main__":
    main()
