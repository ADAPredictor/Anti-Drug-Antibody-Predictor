"""Global configuration. Everything can be overridden with environment vars."""
import os
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
SEED_DATA_PATH = PACKAGE_DIR / "data" / "seed_drugs.json"

# Runtime data lives in ./data by default (database, cached structures ...)
WORK_DIR = Path(os.environ.get("ADA_WORKDIR", Path.cwd()))
DATA_DIR = Path(os.environ.get("ADA_DATA_DIR", WORK_DIR / "data"))
DB_PATH = Path(os.environ.get("ADA_DB", DATA_DIR / "ada.db"))
STRUCTURE_DIR = DATA_DIR / "structures"

# Network
HTTP_TIMEOUT = float(os.environ.get("ADA_HTTP_TIMEOUT", "20"))
USER_AGENT = "ADA-ImmunogenicityPredictor/0.1 (research; contact: local)"

# Optional FDA "Protein Drug Products" catalog (.xlsx). The FDA periodically
# moves this file; set the URL of the current download to enable the importer:
#   ADA_FDA_CATALOG_URL=https://www.fda.gov/media/xxxx/download
FDA_CATALOG_URL = os.environ.get("ADA_FDA_CATALOG_URL", "")

# openFDA label endpoint (no API key required for low-volume use)
OPENFDA_LABEL_URL = "https://api.fda.gov/drug/label.json"
# accessdata.fda.gov blocks automated clients; use the FDA NCTR label service
# first and DailyMed (NLM mirror of SPL) as fallback.
SPL_XML_URL = "https://nctr-crs.fda.gov/fdalabel/services/spl/set-ids/{set_id}/spl-doc"
SPL_XML_FALLBACK_URL = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{set_id}.xml"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Default HLA panel (common supertype representatives, world-wide coverage)
DEFAULT_CLASS1_ALLELES = [
    "HLA-A*01:01",
    "HLA-A*02:01",
    "HLA-A*03:01",
    "HLA-A*11:01",
    "HLA-A*24:02",
    "HLA-B*07:02",
    "HLA-B*08:01",
    "HLA-B*44:03",
]
DEFAULT_CLASS2_ALLELES = [
    "HLA-DRB1*01:01",
    "HLA-DRB1*03:01",
    "HLA-DRB1*04:01",
    "HLA-DRB1*07:01",
    "HLA-DRB1*11:01",
    "HLA-DRB1*15:01",
]

# Affinity cut-offs (nM), Sette-style
STRONG_AFFINITY_NM = 500.0
WEAK_AFFINITY_NM = 5000.0

MAX_SCAN_LENGTH = 3000  # safety guard for the browser workflow


def ensure_dirs():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)
