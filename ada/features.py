"""Physicochemical and liability feature calculation (stdlib only).

All functions work on a single uppercase sequence string.  Use
``compute_features(chains)`` for the aggregate record used by the scorer.
"""
import math

# --- Scales ---------------------------------------------------------------
# Kyte-Doolittle hydropathy (J. Mol. Biol. 1982, 157, 105-132)
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
    "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
    "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
    "Y": -1.3, "V": 4.2,
}

# Chou-Fasman secondary-structure propensities (Adv. Enzymol. 1978)
CF_HELIX = {
    "A": 1.42, "R": 0.98, "N": 0.67, "D": 1.01, "C": 0.70, "Q": 1.11,
    "E": 1.51, "G": 0.57, "H": 1.00, "I": 1.08, "L": 1.21, "K": 1.16,
    "M": 1.45, "F": 1.13, "P": 0.57, "S": 0.77, "T": 0.83, "W": 1.08,
    "Y": 0.69, "V": 1.06,
}
CF_SHEET = {
    "A": 0.83, "R": 0.93, "N": 0.89, "D": 0.54, "C": 1.19, "Q": 1.10,
    "E": 0.37, "G": 0.75, "H": 0.87, "I": 1.60, "L": 1.30, "K": 0.74,
    "M": 1.05, "F": 1.38, "P": 0.55, "S": 0.75, "T": 1.19, "W": 1.37,
    "Y": 1.47, "V": 1.70,
}

# Average residue masses (residue = amino acid - H2O), monoisotopic averages
RESIDUE_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
    "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
    "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
    "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
WATER_MASS = 18.01528

# pKa values (Lehninger / standard ExPASy ProtParam set)
PK_N_TERM = 9.69
PK_C_TERM = 2.34
PK_SIDE = {"K": 10.54, "R": 12.48, "H": 6.00, "D": 3.86, "E": 4.25,
           "C": 8.33, "Y": 10.07}

POSITIVE = set("KRH")
NEGATIVE = set("DE")
POLAR = set("STNQ")
HYDROPHOBIC = set("AVLIMFYW")
AROMATIC = set("FWY")


def molecular_weight(seq):
    return sum(RESIDUE_MASS.get(aa, RESIDUE_MASS["A"]) for aa in seq) + WATER_MASS


def _charge(seq, pH, n_termini=1):
    pos = n_termini * 1.0 / (1.0 + 10 ** (pH - PK_N_TERM))
    neg = n_termini * 1.0 / (1.0 + 10 ** (PK_C_TERM - pH))
    for aa in seq:
        pka = PK_SIDE.get(aa)
        if pka is None:
            continue
        if aa in POSITIVE:
            pos += 1.0 / (1.0 + 10 ** (pH - pka))
        else:
            neg += 1.0 / (1.0 + 10 ** (pka - pH))
    return pos - neg


def isoelectric_point(seq):
    """Bisection search for the pH of zero net charge (single chain)."""
    return isoelectric_point_multi([seq])


def isoelectric_point_multi(seqs):
    """pI of a multi-chain molecule, counting one terminus pair per chain."""
    lo, hi = 0.0, 14.0
    for _ in range(80):
        mid = (lo + hi) / 2.0
        net = sum(_charge(s, mid, n_termini=1) for s in seqs)
        if net > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 2)


def gravy(seq):
    vals = [KD[aa] for aa in seq if aa in KD]
    return round(sum(vals) / len(vals), 3) if vals else 0.0


def extinction_coefficient(seq, all_cystines=True):
    """Pace (1995): e280 = 5500*Trp + 1490*Tyr + 125*cystine."""
    n_w = seq.count("W")
    n_y = seq.count("Y")
    n_c = seq.count("C")
    n_cystine = (n_c // 2) if all_cystines else 0
    return 5500 * n_w + 1490 * n_y + 125 * n_cystine


def count_glycosylation(seq):
    """N-X-S/T sequons where X != P."""
    n = 0
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST":
            n += 1
    return n


def count_motif(seq, motifs):
    return sum(seq.count(m) for m in motifs)


def hydrophobic_patch_positions(seq, window=9, cutoff=1.6):
    """Positions inside strongly hydrophobic windows (aggregation-prone)."""
    flagged = [0] * len(seq)
    if len(seq) < window:
        return flagged
    for i in range(len(seq) - window + 1):
        frag = seq[i:i + window]
        mean = sum(KD.get(aa, 0.0) for aa in frag) / window
        if mean >= cutoff:
            for j in range(i, i + window):
                flagged[j] = 1
    return flagged


def _window_mean(scale, seq, i, half):
    lo = max(0, i - half)
    hi = min(len(seq), i + half + 1)
    vals = [scale.get(aa, 0.0) for aa in seq[lo:hi]]
    return sum(vals) / len(vals)


def secondary_propensity(seq, window=7):
    """Per-residue mean Chou-Fasman helix/sheet propensities."""
    half = window // 2
    helix, sheet = [], []
    for i in range(len(seq)):
        helix.append(round(_window_mean(CF_HELIX, seq, i, half), 3))
        sheet.append(round(_window_mean(CF_SHEET, seq, i, half), 3))
    return helix, sheet


def per_chain(seq):
    n = len(seq)
    ext_red = extinction_coefficient(seq, all_cystines=False)
    ext_cys = extinction_coefficient(seq, all_cystines=True)
    mw = molecular_weight(seq)
    return {
        "length": n,
        "molecular_weight_da": round(mw, 1),
        "pI": isoelectric_point(seq),
        "gravy": gravy(seq),
        "extinction_280_reduced": ext_red,
        "extinction_280_cystines": ext_cys,
        "absorbance_1gpl_reduced": round(ext_red / mw, 3) if mw else 0.0,
        "n_glycosylation": count_glycosylation(seq),
        "n_deamidation_NG_NS_NN": count_motif(seq, ("NG", "NS", "NN")),
        "n_isomerization_DG_DS_DD": count_motif(seq, ("DG", "DS", "DD")),
        "n_acid_clipping_DP": seq.count("DP"),
        "n_methionine": seq.count("M"),
        "n_cysteine": seq.count("C"),
        "n_tryptophan": seq.count("W"),
        "frac_negative": round(sum(aa in NEGATIVE for aa in seq) / n, 4),
        "frac_positive": round(sum(aa in POSITIVE for aa in seq) / n, 4),
        "frac_aromatic": round(sum(aa in AROMATIC for aa in seq) / n, 4),
        "net_charge_pH7": round(_charge(seq, 7.0), 2),
        "hydrophobic_patch": hydrophobic_patch_positions(seq),
    }


def compute_features(chains):
    """Aggregate per-chain features into a record-level feature dict."""
    per = []
    total_len = sum(len(c.seq) for c in chains)
    mw = 0.0
    ext_r = ext_c = 0
    gravy_num = 0.0
    counts = {
        "n_glycosylation": 0,
        "n_deamidation_NG_NS_NN": 0,
        "n_isomerization_DG_DS_DD": 0,
        "n_acid_clipping_DP": 0,
        "n_methionine": 0,
        "n_cysteine": 0,
        "n_tryptophan": 0,
    }
    patches = []  # (offset, flags)
    offset = 0
    for c in chains:
        f = per_chain(c.seq)
        per.append({"chain_id": c.chain_id, "name": c.name,
                    "length": f["length"], "pI": f["pI"],
                    "molecular_weight_da": f["molecular_weight_da"]})
        mw += f["molecular_weight_da"]
        ext_r += f["extinction_280_reduced"]
        ext_c += f["extinction_280_cystines"]
        gravy_num += f["gravy"] * f["length"]
        for k in counts:
            counts[k] += f[k]
        patches.append((offset, c.chain_id, f["hydrophobic_patch"]))
        offset += len(c.seq)

    patch_flags_global = [0] * total_len
    for off, _cid, flags in patches:
        for j, v in enumerate(flags):
            if v:
                patch_flags_global[off + j] = 1

    seqs = [c.seq for c in chains]
    net_charge7 = sum(_charge(s, 7.0, n_termini=1) for s in seqs)
    out = {
        "n_chains": len(chains),
        "length_total": total_len,
        "molecular_weight_da": round(mw, 1),
        "pI": isoelectric_point_multi(seqs),
        "gravy_mean": round(gravy_num / total_len, 3),
        "extinction_280_reduced": ext_r,
        "extinction_280_all_cystines": ext_c,
        "absorbance_1gpl_reduced": round(ext_r / mw, 3) if mw else 0.0,
        "net_charge_pH7": round(net_charge7, 2),
        "frac_hydrophobic_patch": round(sum(patch_flags_global) / total_len, 4),
        "chains": per,
        **counts,
    }
    return out, patch_flags_global
