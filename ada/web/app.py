"""FastAPI application: REST API + single-page web UI.

Prediction engine: immunopred (sequence-intrinsic immunogenicity risk).
Reference database: immunopred ada_reference.csv (289 marketed biologics).
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .. import db
from ..etl import load_seed_data, run_online_enrichment
from ..immunopred import predict_immunogenicity
from ..immunopred.ada_db import ADAReferenceDB
from ..immunopred.sequence import parse_query
from ..immunopred.t_cell import DEFAULT_ALLELES as DEFAULT_MHC2_ALLELES
from ..immunopred.mhc_i import DEFAULT_MHC1_ALLELES

STATIC_DIR = Path(__file__).resolve().parent / "static"

app = FastAPI(
    title="ADA - 蛋白药物免疫原性预测平台",
    description="immunopred 引擎：T/B 细胞表位 + 抗体 CDR/种系校准 + 289 条上市药物 ADA 参照库",
    version="0.2.0",
)


@app.on_event("startup")
def _startup():
    db.init_db()
    load_seed_data()


@app.middleware("http")
async def _no_cache_static(request, call_next):
    response = await call_next(request)
    if request.url.path.startswith("/static") or request.url.path == "/":
        response.headers["Cache-Control"] = "no-cache"
    return response


# ---------------------------------------------------------------- models ---
class PredictRequest(BaseModel):
    sequence: str = Field(..., description="原始氨基酸序列或 FASTA（可多链，如 VH+VL）")
    name: Optional[str] = Field(default=None, description="可选序列名称")
    drug_class: str = Field(
        default="",
        description="药物类别提示：antibody / enzyme / hormone_peptide / cytokine / fc_fusion 等；"
                    "留空则自动推断。")


class ETLRequest(BaseModel):
    limit: int = 25
    overwrite: bool = False


# --------------------------------------------------------- cached reference --
@lru_cache(maxsize=1)
def _ada_db() -> Optional[ADAReferenceDB]:
    try:
        return ADAReferenceDB.from_csv()
    except FileNotFoundError:
        return None


@lru_cache(maxsize=1)
def _ada_summary() -> Dict[str, Any]:
    db_obj = _ada_db()
    if db_obj is None:
        return {"n_entries": 0, "note": "ada_reference.csv 未找到"}
    return db_obj.summary()


# ------------------------------------------------------ CSV drug library ---
CSV_DRUG_PATH = Path(__file__).resolve().parents[2] / "data" / "ada_reference.csv"


@lru_cache(maxsize=1)
def _csv_drugs() -> List[Dict[str, Any]]:
    """Load 289-entry ada_reference.csv as list of dicts for the drug library tab."""
    if not CSV_DRUG_PATH.exists():
        return []
    df = pd.read_csv(CSV_DRUG_PATH)
    drugs = []
    for _, row in df.iterrows():
        name = str(row.get("name", "")).strip()
        if not name:
            continue
        ada_raw = row.get("ada_rate_pct")
        ada_val = float(ada_raw) if pd.notna(ada_raw) else None
        nab_raw = row.get("nab_rate_pct")
        nab_val = float(nab_raw) if pd.notna(nab_raw) else None
        drugs.append({
            "id": _slug(name),
            "name": name,
            "inn": str(row.get("inn", "") or "").strip(),
            "drug_class": str(row.get("drug_class", "") or "").strip(),
            "ada_rate_pct": ada_val,
            "nab_rate_pct": nab_val,
            "nab_note": str(row.get("nab_note", "") or "").strip(),
            "nab_note_zh": str(row.get("nab_note_zh", "") or "").strip(),
            "route": str(row.get("route", "") or "").strip(),
            "approved": str(row.get("approved", "") or "").strip(),
            "source": str(row.get("source", "") or "").strip(),
            "notes": str(row.get("notes", "") or "").strip(),
            "ada_citation": str(row.get("ada_citation", "") or "").strip(),
            "sequence_citation": str(row.get("sequence_citation", "") or "").strip(),
            "heavy": str(row.get("heavy", "") or "").strip(),
            "light": str(row.get("light", "") or "").strip(),
        })
    return drugs


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


@lru_cache(maxsize=1)
def _csv_drug_classes() -> List[str]:
    return sorted({d["drug_class"] for d in _csv_drugs() if d["drug_class"]})


# ----------------------------------------------------------------- pages ---
@app.get("/", include_in_schema=False)
def index():
    return FileResponse(str(STATIC_DIR / "index.html"))


app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


# ------------------------------------------------------------------- API ---
@app.get("/api/health")
def health():
    return {"status": "ok", "reference_entries": _ada_summary().get("n_entries", 0)}


@app.get("/api/meta")
def meta():
    return {
        "engine": "immunopred-0.1",
        "mhc1_backend": "heuristic-pocket-pssm",
        "mhc2_backend": "heuristic-dr-pssm",
        "default_mhc1_alleles": list(DEFAULT_MHC1_ALLELES),
        "default_mhc2_alleles": list(DEFAULT_MHC2_ALLELES),
        "reference_db": _ada_summary(),
        "drug_count": len(_csv_drugs()),
        "drug_classes": _csv_drug_classes(),
    }


@app.post("/api/predict")
def predict(req: PredictRequest):
    try:
        protein = parse_query(req.sequence, seq_id=(req.name or "query"))
    except ValueError as exc:
        return JSONResponse(status_code=400, content={"error": str(exc)})
    try:
        ada_db = _ada_db()
        result = predict_immunogenicity(
            protein,
            ada_db=ada_db,
            use_builtin_db=False,
            mhc_backend="heuristic",
            drug_class=req.drug_class or "",
        )
        return result.to_dict()
    except ValueError as exc:
        return JSONResponse(status_code=400, content={"error": str(exc)})
    except Exception as exc:
        return JSONResponse(status_code=500, content={"error": "预测失败：%s" % exc})


@app.get("/api/drugs")
def drugs(q: Optional[str] = None, drug_class: Optional[str] = None, limit: int = 500):
    """List drugs from the 289-entry CSV reference library."""
    items = _csv_drugs()
    if q:
        ql = q.lower()
        items = [d for d in items if ql in d["name"].lower() or ql in d["inn"].lower()]
    if drug_class:
        items = [d for d in items if d["drug_class"] == drug_class]
    # Sort by ADA rate descending (None last)
    items.sort(key=lambda d: (d["ada_rate_pct"] is not None, d["ada_rate_pct"] or 0), reverse=True)
    return items[:limit]


@app.get("/api/drugs/{drug_id}")
def drug_detail(drug_id: str):
    """Get full detail for a drug from the CSV reference library."""
    for d in _csv_drugs():
        if d["id"] == drug_id:
            return d
    raise HTTPException(status_code=404, detail="药物不存在")


@app.post("/api/etl/online")
def etl_online(req: ETLRequest):
    try:
        return run_online_enrichment(limit=req.limit, overwrite=req.overwrite)
    except Exception as exc:
        return JSONResponse(status_code=502, content={"error": str(exc)})
