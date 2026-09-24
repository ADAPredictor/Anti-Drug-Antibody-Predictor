"""Command-line interface.

Examples
--------
python -m ada.cli initdb
python -m ada.cli etl --online --limit 10
python -m ada.cli predict seq.fasta --modality humanized_mab
python -m ada.cli train
python -m ada.cli serve --port 8000
"""
import argparse
import json
import sys

from . import db
from .config import ensure_dirs, DEFAULT_CLASS1_ALLELES
from .etl import load_seed_data
from .pipeline import analyze


def cmd_initdb(_args):
    ensure_dirs()
    db.init_db()
    ids = load_seed_data()
    print("数据库已初始化：%s（内置参考药物 %d 条）" % (db.DB_PATH, len(ids)))


def cmd_etl(args):
    ensure_dirs()
    db.init_db()
    load_seed_data()
    if args.catalog:
        n = __import__("ada.etl", fromlist=["import_fda_catalog"]).import_fda_catalog()
        print("FDA 目录导入：%d 条" % n)
    if args.online:
        summary = __import__("ada.etl", fromlist=["run_online_enrichment"]).run_online_enrichment(
            limit=args.limit)
        print(json.dumps(summary, indent=2, ensure_ascii=False))
    if not args.online and not args.catalog:
        print("内置参考数据已加载。使用 --online 抓取 FDA 标签，--catalog 导入 FDA 目录。")


def cmd_predict(args):
    with open(args.fasta, "r", encoding="utf-8") as fh:
        text = fh.read()
    result = analyze(
        text,
        modality=args.modality,
        class1_alleles=args.allele or DEFAULT_CLASS1_ALLELES,
        run_class2=not args.no_class2,
        uniprot_acc=args.uniprot,
        pdb_id=args.pdb,
    )
    if args.summary:
        s = result["score"]
        print("总分: %s/100  %s" % (s["total"], s["tier_label"]))
        print("MHC-I 强结合表位段: %d | MHC-II: %d | B 细胞表位段: %d"
              % (result["epitopes"]["mhc1"]["n_segments"],
                 result["epitopes"]["mhc2"]["n_segments"] if result["epitopes"]["mhc2"] else 0,
                 result["epitopes"]["bcell"]["n_segments"]))
        for w in result.get("warnings", []):
            print("提示: " + w)
    else:
        json.dump(result, sys.stdout, ensure_ascii=False, indent=2)
        print()


def cmd_train(_args):
    from .models.trainer import train_model
    print(json.dumps(train_model(), indent=2, ensure_ascii=False))


def cmd_serve(args):
    import uvicorn
    ensure_dirs()
    db.init_db()
    load_seed_data()
    uvicorn.run("ada.web.app:app", host=args.host, port=args.port, reload=False)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="ada", description="蛋白药物免疫原性预测")
    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("initdb").set_defaults(func=cmd_initdb)

    p_etl = sub.add_parser("etl", help="数据加载与抓取")
    p_etl.add_argument("--online", action="store_true", help="从 openFDA/SPL 抓取标签")
    p_etl.add_argument("--catalog", action="store_true", help="导入 FDA 蛋白药物目录 xlsx")
    p_etl.add_argument("--limit", type=int, default=None)
    p_etl.set_defaults(func=cmd_etl)

    p_pred = sub.add_parser("predict", help="对序列运行免疫原性预测")
    p_pred.add_argument("fasta")
    p_pred.add_argument("--modality", default="unknown")
    p_pred.add_argument("--allele", action="append")
    p_pred.add_argument("--no-class2", action="store_true")
    p_pred.add_argument("--uniprot", default=None)
    p_pred.add_argument("--pdb", default=None)
    p_pred.add_argument("--summary", action="store_true")
    p_pred.set_defaults(func=cmd_predict)

    sub.add_parser("train").set_defaults(func=cmd_train)

    p_serve = sub.add_parser("serve")
    p_serve.add_argument("--host", default="127.0.0.1")
    p_serve.add_argument("--port", type=int, default=8000)
    p_serve.set_defaults(func=cmd_serve)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
