"""Optional ML layer: logistic regression on the curated FDA dataset.

Rule-based scoring (``ada.scoring``) is always available.  Once the local
database contains enough labelled protein sequences (seed + FDA catalog +
user data), ``train_model`` fits a regularized logistic regression mapping
the same mechanistic features to observed ADA category (high vs rest) and
persists the coefficients as JSON.  Install scikit-learn to enable it.
"""
import json

from .. import db
from ..config import DATA_DIR
from ..epitopes import bcell, tcell
from ..features import compute_features
from ..scoring import liability_rate, MODALITY_PRIOR
from ..seqio import parse_fasta

MODEL_PATH = DATA_DIR / "model.json"

FEATURE_NAMES = [
    "mhc1_cov_frac", "mhc2_cov_frac", "bcell_frac",
    "hydrophobic_patch_frac", "liability_rate", "modality_prior",
]

MIN_SAMPLES = 10


def build_design_row(sequence_fasta, modality="unknown"):
    chains = parse_fasta(sequence_fasta)
    length = sum(len(c.seq) for c in chains)
    features, _ = compute_features(chains)
    mhc1_cov = mhc2_cov = bcell_cov = 0
    for chain in chains:
        r1 = tcell.scan_class1(chain.seq)
        mhc1_cov += sum(1 for s in r1["segments"] for _ in range(s["start"], s["end"]))
        r2 = tcell.scan_class2(chain.seq)
        mhc2_cov += sum(1 for s in r2["segments"] for _ in range(s["start"], s["end"]))
        rb = bcell.predict_linear(chain.seq)
        bcell_cov += sum(len(x["seq"]) for x in rb["segments"])
    return [
        mhc1_cov / length,
        mhc2_cov / length,
        bcell_cov / length,
        features["frac_hydrophobic_patch"],
        liability_rate(features, length),
        MODALITY_PRIOR.get(modality, MODALITY_PRIOR["unknown"]),
    ]


def train_model():
    try:
        from sklearn.linear_model import LogisticRegression
    except ImportError:
        return {"status": "skipped",
                "reason": "scikit-learn 未安装；运行 pip install scikit-learn 后重试。"}

    X, y, names = [], [], []
    for drug in db.all_drugs():
        if not drug.get("sequence") or not drug.get("ada_category"):
            continue
        try:
            row = build_design_row(drug["sequence"], drug.get("modality") or "unknown")
        except Exception:
            continue
        X.append(row)
        y.append(1 if drug["ada_category"] == "high" else 0)
        names.append(drug["brand"])

    if len(X) < MIN_SAMPLES:
        return {"status": "insufficient_data",
                "reason": "仅有 %d 条带序列的标注药物（至少需要 %d）。"
                          % (len(X), MIN_SAMPLES),
                "usable_drugs": names}
    if len(set(y)) < 2:
        return {"status": "insufficient_data",
                "reason": "标注中只有一个类别，无法训练分类器。",
                "usable_drugs": names}

    model = LogisticRegression(max_iter=1000, class_weight="balanced")
    model.fit(X, y)
    payload = {
        "feature_names": FEATURE_NAMES,
        "coefficients": model.coef_.tolist()[0],
        "intercept": float(model.intercept_[0]),
        "n_samples": len(X),
        "training_drugs": names,
        "training_accuracy": float(model.score(X, y)),
    }
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(MODEL_PATH, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    payload["status"] = "trained"
    return payload
