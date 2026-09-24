"""Continuous B-cell epitope prediction.

Two published propensity scales are used (values reproduced verbatim from
the IEDB B-cell tutorial, https://tools.iedb.org/bcell/help/):

* Parker hydrophilicity (Parker, Guo & Hodges, Biochemistry 1986, 25:5425),
  7-residue arithmetic window mean; hydrophilic surface peaks correlate
  with antigenic sites.
* Chou & Fasman beta-turn propensity (Adv Enzymol 1978;47:45-148), window
  mean, values > 1.0 denote above-average turn propensity (surface loops).

Candidate epitope residues must be both hydrophilic and turn/flexible;
runs >= 5 residues are reported (IEDB convention).  The Parker cut-off was
calibrated against residue RSA from real 3D structures (see calibration
note on PARKER_THRESHOLD).  When an experimental structure is supplied,
:mod:`ada.structure` RSA evidence is folded into the final ensemble score.
"""
WINDOW = 7
MIN_LEN = 5

# Parker et al. 1986 hydrophilicity (IEDB table, order ACDEFGHIKLMNPQRSTVWY)
PARKER = {
    "A": 2.1, "C": 1.4, "D": 10.0, "E": 7.8, "F": -9.2, "G": 5.7,
    "H": 2.1, "I": -8.0, "K": 5.7, "L": -9.2, "M": -4.2, "N": 7.0,
    "P": 2.1, "Q": 6.0, "R": 4.2, "S": 6.5, "T": 5.2, "V": -3.7,
    "W": -10.0, "Y": -1.9,
}

# Chou & Fasman beta-turn propensity (IEDB table)
CF_TURN = {
    "A": 0.66, "C": 1.19, "D": 1.46, "E": 0.74, "F": 0.60, "G": 1.56,
    "H": 0.95, "I": 0.47, "K": 1.01, "L": 0.59, "M": 0.60, "N": 1.56,
    "P": 1.52, "Q": 0.98, "R": 0.95, "S": 1.43, "T": 0.96, "V": 0.50,
    "W": 0.96, "Y": 1.14,
}

# Thresholds calibrated against residue exposure labels (RSA > 0.25,
# Shrake-Rupley, 60-point sphere) pooled from real structures: insulin
# PDB 1MSO, hGH 1HGU, IFN-alpha2 1RH2 and AlphaFold models of G-CSF/EPO
# (1476 residues; base exposed rate 0.69; reproducible via
# tools/calibrate_bcell_rsa.py).  The AND rule Parker >= 1.5 with CF-turn
# >= 1.0 achieved precision 0.83 / 1.20-fold enrichment / 0.23 coverage
# (informedness 0.148), trading sensitivity for specificity.
PARKER_THRESHOLD = 1.5
CF_TURN_THRESHOLD = 1.0


def _window_mean(scale, seq, half):
    out = []
    for i in range(len(seq)):
        lo = max(0, i - half)
        hi = min(len(seq), i + half + 1)
        out.append(sum(scale.get(aa, 0.0) for aa in seq[lo:hi]) / (hi - lo))
    return out


def parker_profile(seq, window=WINDOW):
    """Per-residue Parker hydrophilicity (7-aa centered window mean)."""
    return [round(v, 3) for v in _window_mean(PARKER, seq, window // 2)]


def cf_turn_profile(seq, window=WINDOW):
    """Per-residue Chou-Fasman beta-turn propensity (window mean)."""
    return [round(v, 3) for v in _window_mean(CF_TURN, seq, window // 2)]


def _parker_heat(value):
    """Map raw Parker scale to 0..1 for track rendering (0 at neutral)."""
    v = (value + 2.0) / 7.0  # -2 -> 0, +5 -> 1 (covers relevant range)
    return max(0.0, min(1.0, round(v, 3)))


def predict_linear(seq, parker_threshold=PARKER_THRESHOLD,
                   turn_threshold=CF_TURN_THRESHOLD, min_len=MIN_LEN):
    """Return per-residue heat scores and merged epitope segments.

    Candidate residues satisfy BOTH Parker hydrophilicity >= threshold and
    Chou-Fasman turn propensity >= threshold (independent published lines
    of evidence), which controls the high false-positive rate of either
    scale used alone.
    """
    parker = parker_profile(seq)
    turn = cf_turn_profile(seq)
    scores = [_parker_heat(v) for v in parker]
    flags = [p >= parker_threshold and t >= turn_threshold
             for p, t in zip(parker, turn)]

    segments = []
    run = []
    for i, flag in enumerate(flags + [False]):  # sentinel flushes last run
        if flag:
            run.append(i)
        elif run:
            if len(run) >= min_len:
                start, end = run[0], run[-1] + 1
                mean_p = round(sum(parker[start:end]) / (end - start), 2)
                mean_t = round(sum(turn[start:end]) / (end - start), 2)
                segments.append({
                    "start": start,
                    "end": end,
                    "seq": seq[start:end],
                    "mean_parker": mean_p,
                    "mean_turn": mean_t,
                })
            run = []
    return {
        "scores": scores,
        "parker": parker,
        "turn": turn,
        "segments": segments,
        "window": WINDOW,
        "method": "parker+cf-turn-v2",
        "parker_threshold": parker_threshold,
        "turn_threshold": turn_threshold,
    }
