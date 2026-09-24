"""Importer for the FDA "Protein Drug Products" catalog workbook (.xlsx).

FDA occasionally republishes this dataset (trade name, application type,
approval date and an immunogenicity column).  Because both the URL and the
sheet layout change over time, the importer is deliberately generic: it
scans every sheet for a header row whose cells mention "immunogenicity"
and maps columns by fuzzy name matching.

Activate by exporting ADA_FDA_CATALOG_URL pointing at the current file.
"""
import io
import re

from ..config import FDA_CATALOG_URL, HTTP_TIMEOUT, USER_AGENT

_PCT_RE = re.compile(r"(\d+(?:\.\d+)?)\s*%")

COLUMN_HINTS = {
    "brand": ("trade name", "brand", "product name", "drug name"),
    "generic": ("proper name", "nonproprietary", "generic", "active ingredient",
                "ingredient"),
    "application": ("bla", "nda", "application number"),
    "year": ("approval date", "approval year", "date of approval"),
    "ada": ("immunogenicity", "anti-drug", "ada", "antibody"),
    "modality": ("product type", "modality", "route"),
}


def _http_get(url):
    import requests
    return requests.get(url, headers={"User-Agent": USER_AGENT},
                        timeout=max(HTTP_TIMEOUT, 60))


def _norm(cell):
    return re.sub(r"\s+", " ", str(cell or "").lower()).strip()


def _find_header(rows):
    for ri, row in enumerate(rows[:15]):
        joined = " ".join(_norm(c) for c in row)
        if "immunogenicity" in joined or "anti-drug" in joined:
            mapping = {}
            for ci, cell in enumerate(row):
                name = _norm(cell)
                for key, hints in COLUMN_HINTS.items():
                    if any(h in name for h in hints):
                        mapping.setdefault(key, ci)
            if "brand" in mapping and "ada" in mapping:
                return ri, mapping
    return None, None


def parse_catalog_workbook(data):
    """Best-effort parse -> list of dict rows."""
    import openpyxl
    wb = openpyxl.load_workbook(io.BytesIO(data), read_only=True, data_only=True)
    records = []
    for ws in wb.worksheets:
        rows = [list(r) for r in ws.iter_rows(values_only=True)]
        header_i, mapping = _find_header(rows)
        if mapping is None:
            continue
        for row in rows[header_i + 1:]:
            def cell(key):
                idx = mapping.get(key)
                return row[idx] if idx is not None and idx < len(row) else None
            brand = cell("brand")
            if not brand:
                continue
            ada_text = str(cell("ada") or "")
            values = [float(v) for v in _PCT_RE.findall(ada_text)]
            year = cell("year")
            year_match = re.search(r"(19|20)\d{2}", str(year or ""))
            records.append({
                "brand": str(brand).strip(),
                "generic": str(cell("generic") or "").strip() or None,
                "ada_values": [v for v in values if 0 <= v <= 100],
                "ada_text": ada_text[:300],
                "approval_year": int(year_match.group(0)) if year_match else None,
                "application": str(cell("application") or "").strip() or None,
                "modality_raw": str(cell("modality") or "").strip() or None,
                "source": "fda_catalog",
            })
        break
    return records


def fetch_and_parse(url=None):
    url = url or FDA_CATALOG_URL
    if not url:
        raise RuntimeError(
            "未配置 FDA 目录 URL：设置环境变量 ADA_FDA_CATALOG_URL 后重试。"
        )
    resp = _http_get(url)
    resp.raise_for_status()
    return parse_catalog_workbook(resp.content)
