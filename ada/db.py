"""SQLite storage for curated FDA drugs, import logs and prediction cache."""
import json
import re
import sqlite3
import time
from contextlib import contextmanager

from .config import DB_PATH

SCHEMA = """
CREATE TABLE IF NOT EXISTS drugs (
    id              TEXT PRIMARY KEY,
    brand           TEXT NOT NULL,
    generic         TEXT,
    modality        TEXT,
    approval_year   INTEGER,
    indication      TEXT,
    sequence        TEXT,            -- FASTA, may be NULL
    uniprot_acc     TEXT,
    pdb_id          TEXT,
    ada_min_pct     REAL,            -- observed anti-drug antibody incidence (%)
    ada_max_pct     REAL,
    nab_pct         REAL,            -- neutralizing antibody incidence (%)
    ada_category    TEXT,            -- low / moderate / high
    label_url       TEXT,
    label_evidence  TEXT,            -- excerpt justifying the ADA figures
    source          TEXT,            -- curated | fda_label | fda_catalog
    notes           TEXT,
    updated_at      REAL
);

CREATE TABLE IF NOT EXISTS sources_log (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    source      TEXT,
    status      TEXT,
    n_records   INTEGER,
    message     TEXT,
    fetched_at  REAL
);

CREATE TABLE IF NOT EXISTS prediction_cache (
    cache_key   TEXT PRIMARY KEY,
    payload     TEXT,
    created_at  REAL
);

CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT
);
"""


def slugify(text):
    text = text.lower()
    text = re.sub(r"[^a-z0-9]+", "-", text).strip("-")
    return text


@contextmanager
def connect():
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def init_db():
    with connect() as conn:
        conn.executescript(SCHEMA)
        conn.execute(
            "INSERT OR IGNORE INTO meta(key, value) VALUES('schema_version', '1')"
        )


DRUG_FIELDS = [
    "brand", "generic", "modality", "approval_year", "indication",
    "sequence", "uniprot_acc", "pdb_id", "ada_min_pct", "ada_max_pct",
    "nab_pct", "ada_category", "label_url", "label_evidence", "source", "notes",
]


def upsert_drug(record):
    record = dict(record)
    record.setdefault("id", slugify(record["brand"]))
    record["updated_at"] = time.time()
    fields = ["id"] + DRUG_FIELDS + ["updated_at"]
    placeholders = ",".join("?" for _ in fields)
    updates = ",".join("{f}=excluded.{f}".format(f=f) for f in DRUG_FIELDS + ["updated_at"])
    values = [record.get(f) for f in fields]
    with connect() as conn:
        conn.execute(
            "INSERT INTO drugs (%s) VALUES (%s) ON CONFLICT(id) DO UPDATE SET %s"
            % (",".join(fields), placeholders, updates),
            values,
        )
    return record["id"]


def list_drugs(q=None, modality=None, limit=200):
    sql = "SELECT * FROM drugs WHERE 1=1"
    params = []
    if q:
        sql += " AND (lower(brand) LIKE ? OR lower(generic) LIKE ? OR lower(indication) LIKE ?)"
        like = "%" + q.lower() + "%"
        params += [like, like, like]
    if modality:
        sql += " AND modality = ?"
        params.append(modality)
    sql += " ORDER BY COALESCE(ada_max_pct, 0) DESC, brand ASC LIMIT ?"
    params.append(limit)
    with connect() as conn:
        return [dict(r) for r in conn.execute(sql, params).fetchall()]


def get_drug(drug_id):
    with connect() as conn:
        row = conn.execute("SELECT * FROM drugs WHERE id = ?", (drug_id,)).fetchone()
        return dict(row) if row else None


def all_drugs():
    with connect() as conn:
        return [dict(r) for r in conn.execute("SELECT * FROM drugs ORDER BY brand").fetchall()]


def count_drugs():
    with connect() as conn:
        return conn.execute("SELECT COUNT(*) FROM drugs").fetchone()[0]


def modality_stats():
    """Number of drugs and ADA category breakdown per modality."""
    sql = """
        SELECT modality,
               COUNT(*) AS n,
               SUM(CASE WHEN ada_category='low' THEN 1 ELSE 0 END) AS n_low,
               SUM(CASE WHEN ada_category='moderate' THEN 1 ELSE 0 END) AS n_moderate,
               SUM(CASE WHEN ada_category='high' THEN 1 ELSE 0 END) AS n_high,
               AVG(CASE WHEN ada_min_pct IS NOT NULL AND ada_max_pct IS NOT NULL
                        THEN (ada_min_pct + ada_max_pct) / 2.0 END) AS avg_ada
        FROM drugs GROUP BY modality ORDER BY n DESC
    """
    with connect() as conn:
        return [dict(r) for r in conn.execute(sql).fetchall()]


def log_source(source, status, n_records=0, message=""):
    with connect() as conn:
        conn.execute(
            "INSERT INTO sources_log(source,status,n_records,message,fetched_at) VALUES(?,?,?,?,?)",
            (source, status, n_records, message, time.time()),
        )


def recent_logs(limit=10):
    with connect() as conn:
        rows = conn.execute(
            "SELECT * FROM sources_log ORDER BY id DESC LIMIT ?", (limit,)
        ).fetchall()
        return [dict(r) for r in rows]


def cache_get(key):
    with connect() as conn:
        row = conn.execute(
            "SELECT payload FROM prediction_cache WHERE cache_key=?", (key,)
        ).fetchone()
        return json.loads(row["payload"]) if row else None


def cache_set(key, payload):
    with connect() as conn:
        conn.execute(
            "INSERT OR REPLACE INTO prediction_cache(cache_key,payload,created_at) VALUES(?,?,?)",
            (key, json.dumps(payload), time.time()),
        )
