"""Stdlib unittest smoke tests - run: python -m unittest discover -s tests"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_tmp = tempfile.mkdtemp()
os.environ["ADA_DATA_DIR"] = _tmp
os.environ["ADA_DB"] = os.path.join(_tmp, "test.db")
# keep the unit tests on the fast deterministic heuristic backend
os.environ["ADA_DISABLE_MHCFLURRY"] = "1"

from ada import db
from ada.epitopes import bcell, tcell
from ada.etl import load_seed_data
from ada.features import compute_features
from ada.pipeline import analyze
from ada.scoring import clinical_context
from ada.seqio import parse_fasta, SequenceError
from ada.sources import opfda

INSULIN = """>A
GIVEQCCTSICSLYQLENYCN
>B
FVNQHLCGSHLVEALYLVCGERGFFYTPKT
"""


class TestSequence(unittest.TestCase):
    def test_parse_fasta_two_chains(self):
        chains = parse_fasta(INSULIN)
        self.assertEqual(len(chains), 2)
        self.assertEqual(len(chains[0].seq), 21)
        self.assertEqual(len(chains[1].seq), 30)

    def test_headerless(self):
        chains = parse_fasta("ACDE FGHKL\nMPQR")
        self.assertEqual(chains[0].seq, "ACDEFGHKLMPQR")

    def test_headerless_blankline_blocks_are_chains(self):
        # two blank-line-separated blocks == two chains, identical to
        # parsing the same sequences with ">" headers
        blocks = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTY\n"
                  "\n"
                  "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAV")
        plain = parse_fasta(blocks)
        fasta = parse_fasta(">hc\nEVQLVESGGGLVQPGGSLRLSCAASGFNIKDTY\n"
                            ">lc\nDIQMTQSPSSLSASVGDRVTITCRASQDVNTAV")
        self.assertEqual([c.seq for c in plain], [c.seq for c in fasta])
        self.assertEqual([c.chain_id for c in plain], ["A", "B"])

    def test_sequence_as_header_misformat(self):
        # ">SEQUENCE" (no name) must start a chain, not swallow the sequence
        text = (">EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTY\n"
                ">DIQMTQSPSSLSASVGDRVTITCRASQDVNTAV")
        chains = parse_fasta(text)
        self.assertEqual(len(chains), 2)
        self.assertEqual(chains[0].seq, "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTY")
        self.assertEqual(chains[1].seq, "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAV")

    def test_named_headers_keep_names(self):
        chains = parse_fasta(">trastuzumab heavy chain\nACDEFGHIK\n"
                             ">trastuzumab light chain\nMNPQRSTVWY")
        self.assertEqual(chains[0].chain_id, "trastuzumab")
        self.assertEqual(chains[1].chain_id, "trastuzumab")
        self.assertEqual(chains[0].name, "heavy chain")

    def test_all_seed_sequences_parse(self):
        import json
        seed = os.path.join(os.path.dirname(__file__), "..",
                            "ada", "data", "seed_drugs.json")
        with open(seed, "r", encoding="utf-8") as fh:
            records = json.load(fh)
        n_with_seq = 0
        for rec in records:
            if rec.get("sequence"):
                chains = parse_fasta(rec["sequence"])
                self.assertGreaterEqual(len(chains), 1)
                n_with_seq += 1
        self.assertGreaterEqual(n_with_seq, 25)

    def test_invalid(self):
        with self.assertRaises(SequenceError):
            parse_fasta("ACDEZFG")
        with self.assertRaises(SequenceError):
            parse_fasta("AC")


class TestFeatures(unittest.TestCase):
    def test_insulin_properties(self):
        chains = parse_fasta(INSULIN)
        features, patches = compute_features(chains)
        # human insulin is 5807.8 Da (average mass)
        self.assertAlmostEqual(features["molecular_weight_da"], 5807.8, delta=12.0)
        self.assertTrue(4.8 < features["pI"] < 6.2)
        self.assertEqual(features["length_total"], 51)
        self.assertEqual(features["n_chains"], 2)
        self.assertEqual(features["n_cysteine"], 6)
        self.assertEqual(len(patches), 51)


class TestEpitopes(unittest.TestCase):
    def test_a02_known_binder(self):
        # influenza M1 58-66 GILGFVFTL - canonical strong A*02:01 epitope
        aff = tcell.heuristic_class1_affinity("GILGFVFTL", "HLA-A*02:01")
        self.assertLess(aff, 500.0)
        self.assertEqual(tcell.affinity_level(aff), "strong")

    def test_supertype_parsing(self):
        self.assertEqual(tcell.parse_supertype("HLA-A*02:01"), "A02")
        self.assertEqual(tcell.parse_supertype("HLA-B*07:02"), "B07")
        self.assertEqual(tcell.parse_supertype("HLA-C*07:01"), "generic")

    def test_scan_shapes(self):
        seq = "ACDEFGHIKLMNPQRSTVWYK"  # 21 aa
        r1 = tcell.scan_class1(seq, ["HLA-A*02:01"])
        self.assertEqual(len(r1["track"]), 21)
        self.assertEqual(r1["peptides_scanned"], 13)
        self.assertIn("backend", r1)
        r2 = tcell.scan_class2(seq * 2)
        self.assertEqual(len(r2["track"]), 42)

    def test_bcell(self):
        r = bcell.predict_linear("MMMMMMMMMEEEEEEEEEGGGSSPPPAAAALLVVMMM")
        self.assertEqual(len(r["scores"]), 37)
        for seg in r["segments"]:
            self.assertGreaterEqual(seg["end"] - seg["start"], 5)


class TestLabelExtraction(unittest.TestCase):
    def test_percentages_and_fractions(self):
        sample = (
            "<title>6.2 Immunogenicity</title><text>"
            "Antibodies to X were detected in 9 of 60 patients. "
            "Neutralizing antibodies were found in 3.5% of patients."
            "</text><title>6.3 Postmarketing Experience</title>")
        section = opfda.extract_immunogenicity_section(sample)
        self.assertIn("Antibodies", section)
        vals = opfda.extract_ada_values(section)
        self.assertIn(15.0, vals["percentages"])  # 9/60 fraction parsed
        self.assertTrue(any(3.0 < v < 4.0 for v in vals["neutralizing"]))


class TestPipelineAndDB(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        db.init_db()
        load_seed_data()

    def test_seed_loaded(self):
        drugs = db.all_drugs()
        self.assertGreaterEqual(len(drugs), 28)
        stats = {r["modality"]: r for r in db.modality_stats()}
        self.assertIn("murine_mab", stats)
        # drugs with NULL ADA ranges must not count as zero in the averages:
        # cytokine average is over epoetin (1-2%) and filgrastim (1-3%) only,
        # interferon alfa-2b (null figures) is excluded
        self.assertAlmostEqual(stats["cytokine"]["avg_ada"], 1.75, places=2)

    def test_novelty_similarity_metric(self):
        from ada.pipeline import _chain_similarity
        vh = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTR"
              "YADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
        heavy = vh + "ASTKGPSVFPLAPSSK"  # fragment is a prefix of the full chain
        self.assertEqual(_chain_similarity(vh, vh), 1.0)
        self.assertGreaterEqual(_chain_similarity(vh, heavy), 0.9)
        # an unrelated 20-mer must not look identical through scattered
        # single-residue matches
        self.assertLess(_chain_similarity("WCPRMMMWWWCCCCYYYY", heavy), 0.3)

    def test_novelty_end_to_end(self):
        vh = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTR"
              "YADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
        result = analyze(">trast_VH\n" + vh, modality="humanized_mab",
                         run_class2=False)
        identity = result["score"]["nearest_reference_identity"]
        self.assertIsNotNone(identity)
        self.assertGreaterEqual(identity, 0.99)
        novelty = next(c for c in result["score"]["components"]
                       if c["key"] == "novelty")
        self.assertEqual(novelty["value"], 0.0)

    def test_end_to_end(self):
        result = analyze(INSULIN, modality="hormone_peptide",
                         class1_alleles=["HLA-A*02:01"], run_class2=True)
        score = result["score"]
        self.assertTrue(0 <= score["total"] <= 100)
        self.assertEqual(len(result["epitopes"]["mhc1"]["track"]), 51)
        self.assertIn(result["score"]["tier"],
                      {"low", "low_moderate", "moderate_high", "high"})
        # cache round-trip
        cached = analyze(INSULIN, modality="hormone_peptide",
                         class1_alleles=["HLA-A*02:01"])
        self.assertTrue(cached["cached"])


class TestClinicalContext(unittest.TestCase):
    def test_none_when_nothing_specified(self):
        self.assertIsNone(clinical_context(None))
        self.assertIsNone(clinical_context({}))
        self.assertIsNone(clinical_context({"route": "unknown",
                                            "frequency": "unknown"}))

    def test_sum_and_rows(self):
        ctx = clinical_context({"route": "sc",
                                "frequency": "weekly_or_more",
                                "comedication": "none",
                                "population": "pediatric"})
        self.assertEqual(ctx["adjustment"], 7.0)
        by = dict((r["factor"], r["points"]) for r in ctx["inputs"])
        self.assertEqual((by["route"], by["frequency"], by["comedication"],
                          by["population"]), (2.0, 2.0, 2.0, 1.0))

    def test_invalid_choice_treated_unknown(self):
        ctx = clinical_context({"route": "nasal", "frequency": "daily"})
        self.assertIsNone(ctx)
        ctx = clinical_context({"route": "nasal", "frequency": "weekly_or_more"})
        self.assertEqual(ctx["adjustment"], 2.0)


if __name__ == "__main__":
    unittest.main()
