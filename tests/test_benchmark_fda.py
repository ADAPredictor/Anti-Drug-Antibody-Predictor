# -*- coding: utf-8 -*-
"""FDA-approved protein benchmark (real, traceable sequences).

Ten marketed proteins spanning peptides, hormones and cytokines.  This is a
*regression / behaviour* benchmark, not a validation against ADA rates:
the heuristic epitope backend has no calibrated relationship with clinical
ADA percentages (which are assay-, route- and indication-dependent).

Sequences are mature-chain sequences from UniProt (accessions below),
publicly available; exenatide is exendin-4 (Gila monster).  Scores produced
here are the same outputs the web API returns and can be compared across
backends (heuristic vs mhcflurry) without changing the test.
"""
import unittest

from ada.pipeline import analyze

# (name, uniprot, modality, sequence, total_length)
BENCHMARK = [
    ("human insulin", "P01308", "hormone_peptide",
     ">A\nGIVEQCCTSICSLYQLENYCN\n>B\nFVNQHLCGSHLVEALYLVCGERGFFYTPKT", 51),
    ("glucagon", "P01275", "hormone_peptide",
     "HSQGTFTSDYSKYLDSRRAQDFVQWLMNT", 29),
    ("teriparatide (PTH 1-34)", "P01270", "hormone_peptide",
     "SVSEIQLMHNLGKHLNSMERVEWLRKKLQDVHNF", 34),
    ("exenatide (exendin-4)", None, "hormone_peptide",
     "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS", 39),
    ("pramlintide (amylin analog)", None, "hormone_peptide",
     "KCNTATCATQRLANFLVHSSNNFGPILPPTNVGSNTY", 37),
    ("bivalirudin", None, "other", "FPRPGGGGNGDFEEIPEEYL", 20),
    ("filgrastim (G-CSF)", "P09919", "cytokine",
     "MTPLGPASSLPQSFLLKCLEQVRKIQGDGAALQEKLCATYKLCHPEELVLLGHSLGIPWAPLSSCPSQALQLAGCLSQLHSGLFLYQGLLQALEGISPELGPTLDTLQLDVADFATTIWQQMEELGMAPALQPTQGAMPAFASAFQRRAGGVLVASHLQSFLEVSYRVLRHLAQP", 175),
    ("epoetin alfa (EPO)", "P01588", "cytokine",
     "APPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR", 166),
    ("interferon alfa-2b", "P01563", "cytokine",
     "CDLPQTHSLGSRRTLMLLAQMRRISLFSCLKDRHDFGFPQEEFGNQFQKAETIPVLHEMIQQIFNLFSTKDSSAAWDETLLDKFYTELYQQLNDLEACVIQGVGVTETPLMKEDSILAVRKYFQRITLYLKEKKYSPCAWEVVRAEIMRSFSLSTNLQESLRSKE", 165),
    ("somatropin (hGH)", "P01241", "hormone_peptide",
     "FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTSLCFSESIPTPSNREETQQKSNLELLRISLLLIQSWLEPVQFLRSVFANSLVYGASDSNVYDLLKDLEEGIQTLMGRLEDGSPRTGQIFKQTYSKFDTNSHNDDALLKNYGLLYCFRKDMDKVETFLRIVQCRSVEGSCGF", 191),
]


class TestFDABenchmark(unittest.TestCase):
    def test_all_cases_run_with_sane_outputs(self):
        for name, _acc, modality, seq, length in BENCHMARK:
            r = analyze(seq, modality=modality, run_class2=True,
                        use_cache=False)
            self.assertEqual(r["features"]["length_total"], length, name)
            total = r["score"]["total"]
            self.assertGreaterEqual(total, 0.0, name)
            self.assertLessEqual(total, 100.0, name)
            self.assertIn(r["score"]["tier"],
                          {"low", "low_moderate", "moderate_high", "high"},
                          name)
            for key in ("mhc1", "mhc2", "bcell"):
                ep = r["epitopes"][key] if key != "mhc2" else \
                    r["epitopes"]["mhc2"]
                self.assertEqual(len(ep["track"]), length, "%s:%s" % (name, key))
            cov = r["score"]["coverage"]
            for key in ("mhc1_strong_rate", "mhc2_strong_rate"):
                self.assertGreaterEqual(cov[key], 0.0, name)
                self.assertLessEqual(cov[key], 1.0, name)
            self.assertGreaterEqual(cov["bcell_residues"], 0, name)
            self.assertLessEqual(cov["bcell_residues"], length, name)

    def test_reproducible_with_cache(self):
        name, _acc, modality, seq, _len = BENCHMARK[0]
        r1 = analyze(seq, modality=modality, use_cache=False)
        r2 = analyze(seq, modality=modality, use_cache=True)
        self.assertTrue(r2.get("cached"))
        self.assertEqual(r1["score"]["total"], r2["score"]["total"])

    def test_bcell_scores_in_unit_range(self):
        for name, _acc, modality, seq, length in BENCHMARK:
            r = analyze(seq, modality=modality, use_cache=False)
            for v in r["epitopes"]["bcell"]["track"]:
                self.assertGreaterEqual(v, 0.0, name)
                self.assertLessEqual(v, 1.0, name)


if __name__ == "__main__":
    # manual usage: print the benchmark table
    for name, acc, modality, seq, length in BENCHMARK:
        r = analyze(seq, modality=modality, use_cache=False)
        s, e = r["score"], r["epitopes"]
        print("%-26s L=%-4d %5.1f %-15s MHC1=%d MHC2=%d B=%d" % (
            name, length, s["total"], s["tier"],
            e["mhc1"]["n_segments"],
            e["mhc2"]["n_segments"] if e["mhc2"] else 0,
            e["bcell"]["n_segments"]))
