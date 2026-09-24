# -*- coding: utf-8 -*-
"""Subprocess worker: run mhcflurry class-I prediction under any Python.

The application runs on Python 3.7 but mhcflurry 2.x requires Python >=3.10,
so predictions are delegated to an isolated interpreter (e.g. the embeddable
runtime provisioned by tools/setup_mhcflurry.ps1).  Two modes:

One-shot (stdin = one JSON document, stdout = one JSON document):

    stdin  : {"probe": true}
             or {"requests": [{"allele": "HLA-A*02:01",
                                "peptides": ["SIINFEKL", ...]}, ...]}
    stdout : {"ok": true, "version": "2.2.1",
              "results": [[ic50_nM, ...], ...]}

Persistent (``python mhcflurry_worker.py --serve``): models are loaded once,
then the worker reads one JSON request per stdin line and writes one JSON
response per stdout line, until stdin closes.

* If mhcflurry cannot handle an allele, the corresponding result entry is
  null and the adapter falls back to the heuristic for that column.
* On error: {"ok": false, "error": "..."}.
"""
import contextlib
import io
import json
import os
import sys
from pathlib import Path

# Keep mhcflurry's model downloads inside the project tree by default (the
# upstream default is a per-user AppData directory). Must be set before
# ``import mhcflurry`` because downloads.py reads it at import time.
_DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "tools" / "mhcflurry_data"
os.environ.setdefault("MHCFLURRY_DATA_DIR", str(_DEFAULT_DATA_DIR))

_AFFINITY_COLUMNS = (
    "affinity",
    "mhcflurry_affinity",
    "mhcflurry_presentation_affinity",
    "mhcflurry_presentation_score",
)


def _silenced(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()), \
            contextlib.redirect_stderr(io.StringIO()):
        return fn(*args, **kwargs)


def models_present():
    from mhcflurry import downloads
    downloads.get_path("models_class1_presentation")


def load_predictor():
    """Load the best available trained predictor; raise if none usable."""
    from mhcflurry import Class1PresentationPredictor
    try:
        predictor = _silenced(Class1PresentationPredictor.load)
        return predictor, "presentation"
    except Exception:
        from mhcflurry import Class1AffinityPredictor
        predictor = _silenced(Class1AffinityPredictor.load)
        return predictor, "affinity"


def predict_one(predictor, kind, allele, peptides):
    """IC50-like nM values for one allele (raises on unsupported allele)."""
    if kind == "presentation":
        # In Class1PresentationPredictor the ``alleles`` list is a
        # per-individual genotype (<=6 entries); prediction returns, for
        # every peptide, the best value across the genotype. A single-element
        # genotype therefore gives one row per peptide for THIS allele.
        df = _silenced(predictor.predict,
                       peptides=peptides, alleles=[allele],
                       include_affinity_percentile=False, verbose=0)
        column = next((c for c in _AFFINITY_COLUMNS if c in df.columns), None)
        if column is None:
            raise RuntimeError("unexpected predictor columns: %s"
                               % list(df.columns))
        values = df[column].values
    else:
        values = _silenced(predictor.predict, peptides=peptides, allele=allele)
    if hasattr(values, "tolist"):
        values = values.tolist()
    return [float(v) for v in values]


def handle(predictor, kind, version, payload):
    if payload.get("probe"):
        return {"ok": True, "version": version}
    results = []
    for req in payload.get("requests", []):
        try:
            results.append(predict_one(predictor, kind,
                                       req["allele"], list(req["peptides"])))
        except Exception as exc:  # unsupported allele -> column null
            results.append(None)
            sys.stderr.write("allele %s failed: %s\n"
                             % (req["allele"], exc))
    return {"ok": True, "version": version, "results": results}


def serve():
    import mhcflurry
    version = getattr(mhcflurry, "__version__", "unknown")
    try:
        models_present()
        predictor, kind = load_predictor()
    except Exception as exc:
        sys.stdout.write(json.dumps(
            {"ok": False, "version": version,
             "error": "models not loaded: %s" % exc}) + "\n")
        sys.stdout.flush()
        return 1
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            payload = json.loads(line)
            response = handle(predictor, kind, version, payload)
        except Exception as exc:
            response = {"ok": False, "version": version,
                        "error": "request failed: %s" % exc}
        sys.stdout.write(json.dumps(response) + "\n")
        sys.stdout.flush()
    return 0


def oneshot():
    raw = sys.stdin.read()
    try:
        payload = json.loads(raw)
    except Exception as exc:
        print(json.dumps({"ok": False, "error": "bad json: %s" % exc}))
        return 1

    import mhcflurry
    version = getattr(mhcflurry, "__version__", "unknown")

    if payload.get("probe"):
        try:
            models_present()
        except Exception as exc:
            print(json.dumps({"ok": False, "version": version,
                              "error": "models not fetched: %s" % exc}))
            return 1
        print(json.dumps({"ok": True, "version": version}))
        return 0

    try:
        predictor, kind = load_predictor()
    except Exception as exc:
        print(json.dumps({"ok": False, "version": version,
                          "error": "models not loaded: %s" % exc}))
        return 1

    print(json.dumps(handle(predictor, kind, version, payload)))
    return 0


if __name__ == "__main__":
    sys.exit(serve() if len(sys.argv) > 1 and sys.argv[1] == "--serve"
             else oneshot())
