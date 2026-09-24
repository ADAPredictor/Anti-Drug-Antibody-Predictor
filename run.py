"""Entry point: python run.py [--host 127.0.0.1] [--port 8000]"""
import argparse

from ada.config import ensure_dirs
from ada.db import init_db
from ada.etl import load_seed_data


def main():
    parser = argparse.ArgumentParser(description="ADA immunogenicity prediction server")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--reload", action="store_true")
    args = parser.parse_args()

    ensure_dirs()
    init_db()
    load_seed_data()

    import uvicorn
    uvicorn.run("ada.web.app:app", host=args.host, port=args.port, reload=args.reload)


if __name__ == "__main__":
    main()
