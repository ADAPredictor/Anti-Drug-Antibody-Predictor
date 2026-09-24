"""ETL orchestration: seed data, openFDA enrichment, FDA catalog import."""
import json

from . import db
from .config import SEED_DATA_PATH


def load_seed_data():
    """Upsert the bundled curated FDA drug dataset (idempotent)."""
    with open(SEED_DATA_PATH, "r", encoding="utf-8") as fh:
        records = json.load(fh)
    ids = []
    for rec in records:
        rec = dict(rec)
        rec["id"] = db.slugify(rec["brand"])
        db.upsert_drug(rec)
        ids.append(rec["id"])
    return ids


def _ada_category(values):
    if not values:
        return None
    top = max(values)
    if top < 5:
        return "low"
    if top < 30:
        return "moderate"
    return "high"


def run_online_enrichment(limit=None, overwrite=False):
    """Enrich curated drugs with live FDA SPL label immunogenicity data."""
    from .sources import opfda

    drugs = db.all_drugs()
    if limit:
        drugs = drugs[:limit]
    summary = {"queried": 0, "enriched": 0, "failed": 0, "errors": []}
    for drug in drugs:
        summary["queried"] += 1
        try:
            patch = opfda.enrich_from_label(drug["brand"])
        except Exception as exc:  # network / parse issues must not abort ETL
            summary["failed"] += 1
            summary["errors"].append("%s: %s" % (drug["brand"], exc))
            continue
        if not patch:
            summary["failed"] += 1
            continue
        update = {
            "brand": drug["brand"],
            "label_url": patch["label_url"],
            "label_evidence": patch["label_evidence"] or drug.get("label_evidence"),
            "source": "curated+fda_label",
        }
        values = patch["ada_values"]
        if values and (overwrite or drug.get("ada_min_pct") is None):
            update["ada_min_pct"] = min(values)
            update["ada_max_pct"] = max(values)
            update["ada_category"] = _ada_category(values)
        nab = patch["neutralizing_values"]
        if nab and (overwrite or drug.get("nab_pct") is None):
            update["nab_pct"] = max(nab)
        db.upsert_drug(update)
        summary["enriched"] += 1
    db.log_source("openFDA/SPL labels", "ok" if not summary["failed"] else "partial",
                  summary["enriched"], json.dumps(summary, ensure_ascii=False))
    return summary


def import_fda_catalog(url=None):
    """Download and import the FDA protein-drug catalog workbook."""
    from .sources import catalog

    rows = catalog.fetch_and_parse(url)
    imported = 0
    for row in rows:
        values = row.pop("ada_values")
        ada_text = row.pop("ada_text")
        rec = dict(row)
        if values:
            rec["ada_min_pct"] = min(values)
            rec["ada_max_pct"] = max(values)
            rec["ada_category"] = _ada_category(values)
        rec["label_evidence"] = ada_text
        db.upsert_drug(rec)
        imported += 1
    db.log_source("FDA protein drug catalog", "ok", imported)
    return imported
