# ADA Predictor: An Integrated Sequence-Intrinsic Immunogenicity Risk Predictor for Protein Therapeutics with a 289-Drug Clinical Reference Library

**Authors:** Kaitong Li¹, Min Tang¹, Feilin Wu¹, Mei Yang¹, Lijie Jiang¹, Guoyou Liu¹, Zonglei Zhou¹, Peng Guo¹, Jiaona Guang¹, Jingwen Dong², Zhen Zhang², Lin Pang³, Chunming Rao³

**Affiliations:**
¹ JOINN Medical Testing Laboratories (Beijing) Co., Ltd., Beijing, China
² China Pharmaceutical University, Nanjing, China
³ JOINN Pharmaceutical Quality Research and Testing (Beijing) Co., Ltd., Beijing, China

**Corresponding author:** Kaitong Li, likaitong@joinn-lab.com

---

## Abstract

Immunogenicity — the development of anti-drug antibodies (ADA) against therapeutic proteins — drives reduced efficacy, altered pharmacokinetics, and adverse events. Existing tools address single axes (MHC binding or B-cell epitopes) in isolation, and none integrate sequence-derived T/B-cell epitope signals with antibody CDR/germline calibration and a large clinical ADA reference. We present ADA Predictor, a web platform on the *immunopred* engine that combines MHC-I/II epitope prediction, B-cell linear epitope consensus scoring, antibody CDR delimitation and germline V-gene alignment, multi-signal risk flags with capped uplift, and an independent NAb risk score; a 289-biologic reference library compiled from public sources (IDC-DS-V1 [12], FDA/EMA labels, Thera-SAbDab) provides clinical ADA/NAb rates for background comparison without participating in the sequence score. For infliximab (chimeric anti-TNFα), ADA Predictor reported a sequence-intrinsic risk of 58.2% (High) versus a clinical ADA of 34.23%; for adalimumab (fully human), 14.4% (Low) against a clinical ADA of 36.67% — illustrating that clinical ADA incidence is driven by non-sequence factors (dosing, aggregation, assay) and cannot be inferred from the molecule alone, while CDR/germline calibration distinguished chimeric, humanized, and fully human lineages. ADA Predictor offers transparent, reproducible, sequence-intrinsic risk assessment that complements clinical immunogenicity studies; decoupling molecular from clinical determinants enables honest early-stage candidate screening.

**Keywords:** immunogenicity prediction; anti-drug antibody; MHC epitope; B-cell epitope; antibody CDR; germline; neutralizing antibody; protein therapeutics

---

## 1. Introduction

Therapeutic proteins, particularly monoclonal antibodies (mAbs), have transformed the treatment of cancer, autoimmune disease, and metabolic disorders. However, because these molecules are recognized as foreign by the patient's immune system, a substantial fraction of treated patients develop anti-drug antibodies (ADA) that can neutralize drug activity, alter pharmacokinetics, or precipitate infusion/hypersensitivity reactions [1, 14]. The incidence of ADA varies widely — from <1% for some fully human mAbs to >80% for bacterial-derived enzymes — and is influenced by factors including the drug's molecular origin, route and frequency of administration, concomitant immunosuppression, and the patient's HLA genotype [14].

A growing body of work has sought to predict immunogenicity *in silico*, thereby reducing reliance on expensive and time-consuming wet-lab assays during early drug development. The dominant paradigms include:

- **MHC-binding predictors** (e.g., NetMHCpan [9], IEDB MHC [10]) that rank peptides by predicted HLA binding affinity, focusing on the T-cell epitope axis.
- **B-cell epitope predictors** (e.g., Kolaskar–Tongaonkar antigenicity [1], Parker hydrophilicity [2], Emini surface accessibility [3]) that identify linear epitope segments.
- **Antibody-specific analyses** that examine CDR foreignness and germline proximity [5, 6], recognizing that the variable regions of therapeutic antibodies are a major determinant of immunogenicity.

Each single-axis method has a measurable limit. MHC binding is necessary but not sufficient for immunogenicity: POPI, which selected 23 physicochemical properties to classify 428 human MHC class I binders, reached 64.72% leave-one-out accuracy versus 54.91% and 53.23% for two alignment baselines [22]. Linear B-cell prediction remains weak even with deep learning: EpiDope reported a ROC-AUC of 0.67 ± 0.07 [23]. On the antibody side, humanness scores carry a real but partial signal: Hu-mAb's repertoire-based random forest explained 34% of reported ADA variance (R = −0.58), and BioPhi's OASis 9-mer prevalence score explained 28% (R = −0.53); crucially, restricting OASis to non-human 9-mers that NetMHCIIpan also called binders collapsed explained variance to 0.01 (R = 0.44) [20], showing that naively multiplying foreignness by MHC binding erases rather than strengthens the sequence-ADA signal.

However, these tools are fragmented: a drug developer must run separate analyses for each axis and manually integrate the results. More critically, many approaches attempt to predict *clinical ADA incidence* directly from sequence features, conflating the molecule's intrinsic epitope burden with extrinsic clinical factors (assay sensitivity, dosing regimen, patient population HLA distribution) that are not determinable from sequence alone.

This fragmentation is not incidental but has historical roots: T-cell epitope, B-cell epitope, and antibody-humanness methods belong to different methodological communities, each iterating on its own metric with little incentive for cross-axis integration. Meanwhile, clinical ADA incidence was long scattered across drug labels and the literature in non-machine-readable form, and only recently systematized into a computable dataset [12], precluding integration with sequence tools. The temptation to fit sequence features directly to clinical ADA incidence as a regression target is therefore methodologically appealing — small sample sizes, many dimensions, and a complex noise structure — yet readily yields overfit, uninterpretable models that fold variance from disparate sources into a single target.

We present **ADA Predictor**, an integrated web platform that:

1. Combines MHC-I, MHC-II, B-cell epitope, antibody CDR/germline, and physico-chemical aggregation signals into a single 0–100% **sequence-intrinsic immunogenicity risk** score.
2. Explicitly decouples the sequence-derived risk from clinical ADA incidence: the 289-drug reference library provides clinical ADA/NAb rates for *background comparison only* and does not participate in the scoring computation.
3. Provides transparent, reproducible scoring with bilingual (Chinese/English) annotations, enabling both computational and non-computational users to understand the contribution of each signal.

---

## 2. System and Methods

### 2.1 System Architecture

ADA Predictor is implemented as a FastAPI web application with a single-page frontend. The backend (*immunopred*) is a Python package with 21 core modules handling sequence parsing, feature extraction, calibration, and scoring. The frontend renders results in a card-based UI with collapsible technical explanations, score-composition bar charts, and a risk snapshot table.

The reference library (`ada_reference.csv`) contains 289 marketed biologics with fields including: drug name, INN, drug class, ADA rate (%), NAb rate (%), administration route, approval status, sequence (heavy/light for antibodies), and literature citations. Data are compiled from multiple public sources: the IDC-DS-V1 clinical immunogenicity dataset (CC BY 4.0 [12]), FDA/EMA drug labels (public-domain regulatory data), Thera-SAbDab antibody sequences, and UniProt/PDB records. Of all therapeutic molecules in IDC-DS-V1, 96 were included; the others lacked complete information and were excluded.

### 2.2 MHC-I Epitope Prediction

A 9-mer sliding window scans the input sequence against pocket-preference matrices for 7 common HLA-A/B alleles (HLA-A\*02:01, A\*24:02, A\*03:01, B\*07:02, B\*08:01, B\*35:01, B\*44:02). Each matrix encodes position-specific amino-acid preferences at anchor positions (P2 and the C-terminus). Peptides are scored and ranked by intra-protein percentile; the top 10% with a matrix score ≥0.15 are flagged as strong binders.

### 2.3 MHC-II Epitope Prediction

A 15-mer sliding window with a 9-mer core scan is used against position-specific scoring matrices (PSSM) for 6 HLA-DRB1 alleles (DRB1\*01:01, \*03:01, \*04:01, \*07:01, \*11:01, \*15:01). The PSSM encodes preferences at P1, P4, P6, and P9 pocket positions. The strong-binding threshold is set at allele-specific threshold × 0.5 with a percentile ≥90.

### 2.4 B-Cell Linear Epitope Consensus Scoring

Four classical algorithms are combined as a z-score mean:

| Component | Method | Reference |
|---|---|---|
| Antigenicity | Kolaskar–Tongaonkar | [1] |
| Hydrophilicity | Parker et al. | [2] |
| Surface accessibility | Emini et al. | [3] |
| β-turn propensity | Chou–Fasman | [4] |

Consecutive segments of ≥6 amino acids with a consensus z-score ≥0.5 are identified as candidate B-cell epitopes. We note that >90% of B-cell epitopes are conformational; linear scoring serves only as a candidate-region screen.

### 2.5 Antibody CDR Delimitation and Germline Calibration

CDRs are delimited by conserved Cys/Trp residues and WGXG/FGXG anchor motifs (Chothia-like [5]). The variable-region sequence is aligned against a bundled IMGT human V-gene pool [6] using Biopython [8] local alignment to obtain:

- Framework (FR) mismatch fraction
- CDR mismatch fraction
- Best germline V-gene hit (chain, gene, % identity)

Lineage (chimeric/humanized/fully human/murine) is inferred from both the INN suffix (-ximab/-zumab/-umab [11, 13]) and the germline FR identity, with user-alerting on modality mismatches.

### 2.6 Multi-Signal Risk Flags and Capped Uplift

Nine risk flags are evaluated (Table 1). Single MHC/B-cell epitope density flags alone do not trigger uplift; germline-overlap or framework-distance corroboration is required. The uplift is capped at +15% to prevent runaway amplification.

The mechanistic rationale is that single-axis epitope-density signals carry a high false-positive rate — any protein contains some peptides of high antigenicity or binding potential. Requiring that such an epitope also fall on germline-mismatch residues (indicating it resides on a foreign segment) or on a framework distant from the low-ADA human pool materially raises the posterior probability that it is presented and triggers T-cell recognition. The cap prevents any single strong flag from dominating the composite through repeated amplification, keeping the final risk anchored to the raw epitope burden.

| Flag ID | Signal (EN) |
|---|---|
| t_cdr | Elevated CDR non-human-germline epitope burden |
| b_epi | Elevated B-cell linear epitope coverage |
| germline_hetero | Elevated heterologous residues vs best-hit V germline |
| t_on_hetero | MHC-II peptides overlap germline CDR mismatch residues |
| b_on_hetero | B-cell epitopes overlap germline CDR mismatch residues |
| agg | Elevated aggregation-related risk |
| fr_distant | Framework distant from low-ADA human pool |
| charge | Elevated absolute charge density |
| hetero_cdr | Elevated CDR foreignness proxy (fallback) |

### 2.7 Scoring Model

**Step 1 — Raw epitope signal weighting:**

$$\text{FullSeqEpitope} = 0.45 \cdot T_{\text{density}} + 0.20 \cdot T_{\text{score}} + 0.20 \cdot B_{\text{coverage}} + 0.08 \cdot M1_{\text{density}} + 0.07 \cdot A_{\text{agg}}$$

**Step 2 — Path calibration:**

- *Antibody path:* Sequence foreignness (scaffold + FR/CDR mismatch + format exposure) × CDR-weighted epitope burden × lineage tolerance factor (fully human 0.55, humanized 0.75, chimeric 1.15) → baseline; multi-signal uplift added (cap +15%).
- *Protein path:* Class-specific foreignness × epitope × class tolerance factor.
- *Heuristic path:* Near-endogenous ×0.35; human-like Ab ×0.40.

The lineage tolerance factor encodes the a priori expectation of immune recognition from the source of the variable region: fully human antibodies derive from human germline V-genes within the immune-tolerance repertoire; humanized antibodies retain murine CDRs grafted onto human frameworks; chimeric antibodies retain entire murine variable domains. This factor multiplies the foreignness (FR/CDR mismatch) measured in the previous step, so the final risk is jointly constrained by both the "source prior" and the "measured foreignness," rather than by either alone — avoiding the one-sidedness of classifying by lineage only or by residue alignment only.

The antibody path decomposes foreignness and epitope burden into named, weighted components before the lineage factor is applied.

**Foreignness (0–100)** mixes five components: germline scaffold distance (weight 0.35), heterologous residue fraction (0.20), CDR mismatch to the best human V-gene hit (0.20), excess CDR divergence beyond the 0.78 germline cutoff (0.05), and format exposure (0.20).

**Epitope burden (0–100)** mixes CDR-local MHC-II density (0.32), non-germline residue burden (0.33), linear B-cell coverage (0.20), and a framework-distance term (0.15). Its relevance to the final score is gated by foreignness:

$$\text{relevance} = 0.20 + 0.80 \cdot \frac{\text{foreignness}}{100}$$

so burden contributes less when the variable region is close to human germline. The baseline score is:

$$\text{base} = \text{clip}\left(0.45 \cdot \text{foreignness} \cdot \lambda + 0.55 \cdot \text{epitope} \cdot \text{relevance},\ 1,\ 80\right)$$

Soft ceilings apply to highly human scaffolds: fully human scaffold < 0.40 capped at 22; humanized scaffold < 0.45 capped at 16; humanized scaffold < 0.55 capped at 19; chimeric baseline capped at 58. These are structural priors, not coefficients fit to the 289 clinical ADA rates.

**Step 3 — Corroborated uplift and mechanism context:** The nine risk flags (Table 1) contribute a capped uplift: a lone MHC or B-cell density flag adds 0; overlap of a strong MHC-II peptide or B-cell segment with germline-mismatch residues adds a small uplift; two concurrent flags add 7 points; three or more add 12; the total is capped at +15. Fully human and humanized lineages further scale this uplift down unless framework distance is high.

A separate context multiplier is then applied for antibody class: CD19/CD20 B-cell depleting antibodies ×0.12 (cap 6); bispecific, ADC or T-cell engager ×0.48 (cap 28); oncology mAb ×0.40 (cap 24); acute fragment ×0.40 (cap 18). This is the only place where a class-level clinical pattern informs the sequence score; it is a hand rule, not a regression residual from any single drug's ADA rate.

**Classification:** <20% Low; 20–44.9% Moderate; ≥45% High.

### 2.8 NAb Risk Score

An independent three-factor weighted score:

$$\text{NAb} = w_1 \cdot S_{\text{structure}} + w_2 \cdot A_{\text{ADA-incidence}} + w_3 \cdot L_{\text{literature}}$$

- **Structure/class factor:** Baseline 25 + format-exposure bonus + foreignness bonus + class-specific rules (ERT +35, bacterial/toxin +40, chimeric +25, humanized mAb −10, insulin/hormone −15).
- **ADA-incidence factor:** Derived from the nearest neighbor drug's clinical ADA rate.
- **Literature factor:** From IDC-DS-V1 NAb annotations and published reports.

Levels: Low (<25%), Moderate (25–55%), High (≥55%).

### 2.9 Reference Library and kNN Matching

Sequence matching is chain-aware:
- **Fv (heavy + light) ≥95% similarity:** high-confidence same-drug hit.
- **VH-only or VL-only ≥85%:** labeled "single-chain framework only" with a warning that human Ab frameworks are commonly >85% similar.
- **k-mer local fragment:** hint only.

Clinical ADA/NAb values from matched drugs are displayed for background comparison and used in the NAb incidence factor, but are never fed back into the sequence risk score.

### 2.10 Benchmark Against Established Predictors

To position ADA Predictor's built-in pocket/PSSM heuristics against established neural-network predictors, we queried the public IEDB Analysis Resource API (NetMHCpan for MHC-I and NetMHCIIpan for MHC-II [9, 16]) on a panel of nine marketed mAbs with validated VH/VL sequences in the reference library, spanning three lineages: chimeric (infliximab, rituximab, cetuximab), humanized (trastuzumab, bevacizumab, omalizumab), and fully human (adalimumab, panitumumab, golimumab). For every heavy and light chain (18 chains, 9,389 scored windows across 90 chain–allele comparisons; MHC-I alleles HLA-A\*02:01, A\*01:01, B\*07:02; MHC-II alleles HLA-DRB1\*01:01, \*04:01), we compared (i) per-peptide rank agreement (Spearman's ρ between ADA Predictor's internal score and the inverted IEDB percentile rank), (ii) overlap of each method's top-decile peptides (Jaccard index), and (iii) strong-binder counts under each method's threshold (ADA Predictor: intra-protein top 10% plus the engine's absolute threshold; NetMHCpan/NetMHCIIpan: percentile rank ≤2). The fully reproducible script and raw API responses are archived in `paper/benchmark_iedb.py` and `paper/_iedb_cache.json`.

---

## 3. Results

### 3.1 Chimeric Antibody Case: Infliximab

Infliximab (Remicade®) is a chimeric anti-TNFα IgG1 mAb with murine variable regions and human constant regions. Input of its VH and VL sequences yielded:

| Metric | Value |
|---|---|
| Sequence-intrinsic risk | 58.2% (High) |
| Clinical ADA reference | 34.23% |
| Clinical NAb reference | 15.0% |
| Predicted NAb risk | 23.2% (Low) |
| Inferred lineage | Chimeric |
| Scoring path | Antibody CDR/germline calibration |
| Sequence foreignness | 49.5% |
| Base score | 46.2% |
| Multi-signal uplift | +12.0% (capped at +15%) |
| Full-length epitope load | 68.9% (pre-calibration) |

The high risk was driven by elevated CDR non-human-germline epitope burden (foreignness 49.5%) and framework distance from the low-ADA human pool. The multi-signal uplift (+12.0%) was triggered by 5 concurrent flags: CDR non-human-germline epitope burden, germline V-gene divergence, MHC-II peptides overlapping germline-mismatch CDR residues, B-cell epitopes overlapping germline-mismatch CDR residues, and framework distant from the low-ADA human pool.

### 3.2 Fully Human Antibody Case: Adalimumab

Adalimumab (Humira®) is a fully human anti-TNFα mAb. The system predicted:

| Metric | Value |
|---|---|
| Sequence-intrinsic risk | 14.4% (Low) |
| Clinical ADA reference | 36.67% |
| Clinical NAb reference | 5.0% |
| Predicted NAb risk | 17.0% (Low) |
| Inferred lineage | Fully human |
| Lineage tolerance factor | 0.55 |
| Sequence foreignness | 20.7% |
| Multi-signal uplift | 0.0% (no uplift triggered) |
| Full-length epitope load | 62.8% (pre-calibration) |

Despite a pre-calibration epitope load of 62.8% — comparable to infliximab's 68.9% — the fully human lineage tolerance factor (0.55) and low foreignness (20.7%) calibrated the final risk down to 14.4%. Notably, adalimumab's clinical ADA of 36.67% (patient-weighted mean across 159 cohorts) in the IDC-DS-V1 dataset is *higher* than its sequence-intrinsic risk, a divergence attributable to non-sequence factors (subcutaneous route, chronic high-dose administration, aggregation propensity, and assay methodology) rather than molecular immunogenicity. This case exemplifies the core design principle: sequence-intrinsic risk and clinical ADA incidence measure different things and should not be conflated.

### 3.3 Reference Library Coverage

The 289-entry library spans the following drug classes:

| Class | Count |
|---|---|
| Fully human mAb | 26 |
| Humanized mAb | 68 |
| Chimeric mAb | 19 |
| Other antibody / fragment | 3 |
| Enzyme (ERT / thrombolytic) | 37 |
| Hormone / peptide | 58 |
| Cytokine / growth factor | 30 |
| Coagulation factor | 22 |
| Antibody-drug conjugate (ADC) | 7 |
| Fc fusion protein | 6 |
| Other (toxin, plasma protein) | 13 |

*Lineage (fully human / humanized / chimeric) was inferred from the framework-region heterologous fraction against the human germline V-gene pool; non-antibody classes are grouped by the `drug_class` field.*

### 3.4 Benchmark Against NetMHCpan / NetMHCIIpan

Across the nine-mAb panel, the lightweight pocket/PSSM heuristics inside ADA Predictor showed moderate per-peptide ranking agreement with the neural-network predictors: Spearman's ρ averaged 0.485 (median 0.473, SD 0.108) for MHC-I and 0.385 (median 0.399, SD 0.135) for MHC-II over 90 chain–allele comparisons (9,389 peptides); 89 of 90 comparisons were positively correlated, and 96% (MHC-I) / 75% (MHC-II) exceeded ρ = 0.3 (Table 2). The top-decile peptide sets overlapped with Jaccard 0.334 (MHC-I) and 0.250 (MHC-II). Agreement was highest for HLA-A\*02:01 (ρ = 0.602), the best-characterized allele.

**Table 2.** Per-peptide agreement between ADA Predictor internal heuristics and IEDB predictors on 9 marketed mAbs (18 VH/VL chains).

| Axis | Comparisons | Total peptides | ρ mean (median ± SD) | Top-decile Jaccard | ρ > 0.3 |
|---|---|---|---|---|---|
| MHC-I vs NetMHCpan | 54 | ~5,600 | 0.485 (0.473 ± 0.108) | 0.334 | 96% |
| MHC-II vs NetMHCIIpan | 36 | ~3,800 | 0.385 (0.399 ± 0.135) | 0.250 | 75% |

Per-allele mean ρ: HLA-A\*02:01 0.602, HLA-A\*01:01 0.447, HLA-B\*07:02 0.405; HLA-DRB1\*01:01 0.424, DRB1\*04:01 0.347.

ADA Predictor flagged approximately 3.2-fold more strong binders than the IEDB percentile-rank ≤2 threshold (MHC-I 593 vs 187; MHC-II 380 vs 119), indicating higher sensitivity but lower specificity — appropriate for an early-screening heuristic, but it should not be used as a substitute for NetMHCpan for allele-specific deimmunization decisions.

At the composite level, the nine-drug panel is illustrative rather than a validation set (n = 9). ADA Predictor sequence-intrinsic risk retained the expected molecular ordering — chimeric 29.1% > humanized 15.6% > fully human 11.2% group means — whereas clinical ADA incidence did not (14.1%, 2.1%, 16.3%; Spearman ρ = −0.13 between sequence risk and clinical ADA). The two largest discordances are mechanistically interpretable: rituximab (chimeric, but a CD20 B-cell-depleting antibody; sequence score suppressed 12-fold by mechanism-based context tolerance while clinical ADA stays ~5%) and adalimumab (fully human sequence, yet 36.67% clinical ADA in this dataset, reflecting dosing, aggregation, and assay factors outside sequence). This divergence is the expected consequence of the design principle stated in §2.9 and mirrors published observations that even production-grade composite scores explain only a fraction of clinical ADA variance — a 2026 industry review reported that Novo Nordisk's iPREDICT immunogenicity risk score correlated r = 0.35 with clinical ADA and inverted humanness r = 0.50 on 217 mAbs (secondary report from a 2024 conference presentation [19]), and BioPhi OASis humanness explained ~28% variance of clinical ADA in its benchmark [20].

---

## 4. Discussion

### 4.1 Comparison with Existing Tools

The immunogenicity-prediction landscape is fragmented across three layers (Table 3). **MHC/T-cell tools** are the most mature: NetMHCpan/NetMHCIIpan [9, 16] and the IEDB [10] provide accurate per-allele binding and eluted-ligand predictions; PickPocket [25] extends coverage to under-characterized alleles by interpolating from receptor-pocket similarity, complementing the neural-net predictors. Recent deep-learning systems such as T-SCAPE [18] predict T-cell activation and even ADA-inducing potential without explicit MHC inputs. **B-cell tools** remain axis-specific, with linear methods (BepiPred) and conformational structure-based methods such as DiscoTope-3.0 [17]; even the latter achieve only modest F1 (~0.2–0.3). **Antibody-level tools** assess humanness or classify whole molecules: BioPhi OASis/Sapiens [20] provides interpretable 9-mer humanness against the natural OAS repertoire (humanness classification ROC-AUC 0.966; ~28% explained variance of clinical ADA); Hu-mAb and TAP [24] address humanization and developability on top of the ANARCI numbering framework [26]; and AbImmPred [21] classifies antibody immunogenicity from AntiBERTy embeddings (72.7% accuracy on a 22-mAb independent test). Pharmaceutical pipelines (Novo Nordisk iPREDICT; Novartis iSHAPe; EpiVax EpiMatrix) integrate MHC presentation with self-similarity and population weighting but are largely proprietary [19].

**Table 3.** Feature coverage of representative immunogenicity-related tools.

| Tool | MHC-I | MHC-II | B-cell | CDR/germline humanness | Aggregation | NAb | Clinical reference | Integrated 0–100 score | Clinical ADA in score | Access |
|---|---|---|---|---|---|---|---|---|---|---|
| NetMHCpan 4.1 [9] | ✓ | – | – | – | – | – | – | – | – | free (IEDB) |
| NetMHCIIpan 4.x [16] | – | ✓ | – | – | – | – | – | – | – | free (IEDB) |
| IEDB toolkit [10] | ✓ | ✓ | ✓ | – | – | – | – | – | – | free |
| BepiPred / DiscoTope | – | – | ✓ | – | – | – | – | – | – | free web |
| DiscoTope-3.0 [17] | – | – | ✓ (3D) | – | – | – | – | – | – | free server |
| BioPhi (OASis) [20] | – | – | – | ✓ | – | – | partial | – | – | open source |
| TAP / Hu-mAb (OPIG) [24] | – | – | – | ✓ | ✓ | – | – | – | – | free web |
| AbImmPred [21] | – | – | – | ✓ | – | – | yes (labels) | ✓ (class) | yes | web |
| iPREDICT / iSHAPe [19] | ✓ | ✓ | – | ✓ | – | – | yes | ✓ | yes | proprietary |
| T-SCAPE [18] | ✓ | ✓ | – | – | – | – | yes | ✓ | yes | research |
| **ADA Predictor (this work)** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | **289-drug, decoupled** | ✓ | **no (decoupled)** | open web |

ADA Predictor's distinctive position is not per-allele accuracy — the benchmark in §3.4 confirms that the neural predictors remain more accurate per-allele — but *coverage and interpretability*: a single, dependency-light, fully transparent platform that spans T-cell, B-cell, antibody germline, aggregation, and an independent NAb axis, with every composite component shown as an editable weighted bar and every threshold documented to reproducibility. Critically, ADA Predictor is the only platform that presents a large clinical ADA/NAb reference library while *structurally excluding* clinical rates from the sequence score, preventing the optimistic and confounded "sequence-features predict clinical outcome" framing that small-n studies invite.

### 4.2 Sequence Risk vs. Clinical ADA

A key design principle is the strict separation of sequence-intrinsic risk from clinical ADA incidence. Clinical ADA rates are influenced by non-sequence factors (route, frequency, concomitant immunosuppression, assay sensitivity, patient HLA distribution) that are not determinable from the molecule alone. By displaying clinical rates as background reference only — never feeding them into the score — ADA Predictor avoids the methodological error of predicting clinical outcomes from partial molecular information.

A dialectical caveat is warranted: published work does show that sequence features carry *partial* predictive signal for clinical ADA. BioPhi OASis humanness explained ~28% of clinical ADA variance in its benchmark [20], AbImmPred's AntiBERTy-embedding classification reached 72.7% accuracy on 22 mAbs [21], and a 2026 industry review reported iPREDICT's immunogenicity risk score correlating r = 0.35 with clinical ADA on 217 mAbs (secondary report [19]). These signals are real and useful, yet 28% explained variance leaves 72% attributable to non-sequence factors, and a classification on n = 22 is preliminary. ADA Predictor's position is therefore not that "sequence is irrelevant to clinical ADA," but that "sequence is necessary yet insufficient information" — displaying it separately from clinical incidence preserves the explanatory power of each, rather than forcing a weak correlation into a deterministic prediction. The NAb score's ADA-incidence factor is the sole exception, where the nearest neighbor's clinical ADA rate serves as a *proximity-weighted prior* rather than a direct predictor.

### 4.3 Limitations

1. **MHC binding heuristics:** The pocket-preference matrices are not as accurate as NetMHCpan's neural-network models for specific allele-peptide pairs.
2. **Linear B-cell epitopes:** >90% of B-cell epitopes are conformational; linear scoring has limited surface-discrimination power.
3. **CDR delimitation:** Anchor-based heuristics approximate — but do not replace — IMGT/ANARCI formal numbering.
4. **No 3D structure:** Format-exposure is inferred from sequence architecture, not crystal-structure docking.
5. **Clinical data variability:** ADA/NAb rates in the reference library are affected by heterogeneous assay methods and trial designs; cross-drug comparison requires caution.

### 4.4 Future Directions

- Addressing Limitation 1 (§4.3) and the §3.4 benchmark (MHC-I ρ = 0.485, lower specificity than NetMHCpan), integrate NetMHCpan-4.1 as an optional MHC backend to provide higher accuracy for allele-specific deimmunization, while retaining the built-in heuristics as a lightweight default.
- Addressing Limitation 3, replace the anchor heuristic with ANARCI for formal IMGT numbering to meet the traceability requirements of regulatory filings.
- Addressing Limitation 2 (>90% of B-cell epitopes are conformational), introduce conformational B-cell epitope prediction via structure-prediction models to compensate for the surface-discrimination limits of linear scoring.
- Addressing Limitation 5 and the emergence of new modalities, expand the reference library to include biosimilars and bispecifics/CAR-T constructs, both increasing the clinical-background sample size and providing decoupled references for new modalities.

---

## Acknowledgements

The authors thank the ADA Predictor development team for their collaborative effort in designing, implementing, and validating the platform. We gratefully acknowledge the open-source tools and public resources that made this work possible: Biopython [8] for sequence alignment, NumPy and pandas for numerical computation, FastAPI for the web service, and the IEDB [10], IMGT [6], and PDB [7] databases for immunological and structural reference data. MHC binding predictions rely on NetMHCpan/NetMHCIIpan [9, 16] and the optional mhcflurry backend. The 289-drug clinical reference table compiles data from IDC-DS-V1 [12], FDA/EMA labels, and Thera-SAbDab, and is released under CC BY 4.0.

## Funding

No funding to declare.

## Conflict of Interest

The authors declare no competing interests.

## Data and Code Availability

The *immunopred* source code is available at https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor. The 289-drug clinical reference table used in this study compiles clinical ADA/NAb rates from the IDC-DS-V1 dataset (Immunogenicity Database Collaborative, CC BY 4.0 [12], https://immunogenicitydb.org/), FDA/EMA drug labels (public-domain regulatory data), and antibody sequences from Thera-SAbDab and UniProt/PDB. This reference table, together with the benchmark composite values used in the case studies (`paper/benchmark_composite.csv`), is available at the same GitHub repository under CC BY 4.0. IEDB percentile-rank caches are derived from public IEDB data.

---

## References

1. Kolaskar AS, Tongaonkar PC. A semi-empirical method for prediction of antigenic determinants on protein antigens. *FEBS Lett.* 1990;276(1-2):172-176.
2. Parker JM, Guo D, Hodges RS. New hydrophilicity scale derived from high-performance liquid chromatography peptide retention data. *Biochemistry.* 1986;25(19):5425-5432.
3. Emini EA, Hughes JV, Perlow DS, Boger J. Induction of hepatitis A virus-neutralizing antibody by a virus-specific synthetic peptide. *J Virol.* 1985;55(3):836-839.
4. Chou PY, Fasman GD. Prediction of the secondary structure of proteins from their amino acid sequence. *Adv Enzymol.* 1978;47:45-148.
5. Chothia C, Lesk AM. Canonical structures for the hypervariable regions of immunoglobulins. *J Mol Biol.* 1987;196(4):901-917.
6. Lefranc MP, Giudicelli V, Ginestoux C, et al. IMGT, the international ImMunoGeneTics information system. *Nucleic Acids Res.* 2009;37(suppl_1):D1006-D1012.
7. Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic Acids Res.* 2000;28(1):235-242.
8. Cock PJ, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics.* 2009;25(11):1422-1423.
9. Nielsen M, Lundegaard C, Blicher T, et al. NetMHCpan, a method for quantitative prediction of peptide binding to any HLA-A and -B locus. *PLoS One.* 2008;3(3):e1797.
10. Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database (IEDB): 2018 update. *Nucleic Acids Res.* 2019;47(D1):D339-D343.
11. World Health Organization. International Nonproprietary Names (INN) for biological and biotechnological substances. *INN Working Document 05.179.* 2005.
12. Immunogenicity Database Collaborative. IDC DS V1: Clinical immunogenicity dataset of marketed biologics. *CC BY 4.0.* https://immunogenicitydb.org/
13. Mayrhofer P, Kunert R. Nomenclature of humanized mAbs: early concepts, current challenges and future perspectives. *Hum Antibodies.* 2018;27(1):37-51.
14. Jawa V, Cousens LP, Awwad M, Wakshull E, Kropshofer H, De Groot AS. T-cell dependent immunogenicity of protein therapeutics: preclinical assessment and mitigation. *Clin Immunol.* 2013;149(3):534-555.
15. Rosenberg AS, Sauna ZE. Immunogenicity assessment during the development of protein therapeutics. *J Pharm Pharmacol.* 2018;70(5):584-594.
16. Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions across MHC alleles using refined data and analysis. *Nucleic Acids Res.* 2020;48(W1):W449-W454.
17. Sanchez-Trincado JL, Gomez-Perosanz M, Reche PA. Fundamentals and methods for T- and B-cell epitope prediction. *J Immunol Res.* 2017;2017:2680160.
18. Kim J, et al. T-SCAPE: a novel T-cell activation prediction method identifying immunogenic regions in protein therapeutics. *Sci Adv.* 2025;11(49):eadz8759.
19. Hickling TP, Nielsen M, Meysman P, Rose RH, Obrezanova O. Review: application and opportunities for machine learning and artificial intelligence in preclinical immunogenicity risk assessment. *Front Immunol.* 2026;17:1720928.
20. Prihoda D, Maamary J, Waight A, Juan V, Fayadat-Dilman L, Svozil D, Bitton DA. BioPhi: a platform for antibody design, humanization, and humanness evaluation based on natural antibody repertoires and deep learning. *MAbs.* 2022;14(1):2020203.
21. Wang H, Hao X, He Y, Fan L. AbImmPred: an immunogenicity prediction method for therapeutic antibodies using AntiBERTy-based sequence features. *PLoS One.* 2024;19(2):e0296737.
22. Tung CW, Ho SY. POPI: predicting immunogenicity of MHC class I binding peptides by mining informative physicochemical properties. *Bioinformatics.* 2007;23(8):942-949.
23. Collatz M, Mock F, Barth E, Hölzer M, Sachse K, Marz M. EpiDope: a deep neural network for linear B-cell epitope prediction. *Bioinformatics.* 2021;37(4):448-455.
24. Marks C, Hummer AM, Chin M, Deane CM. Humanization of antibodies using a machine learning approach on large-scale repertoire data. *Bioinformatics.* 2021;37(22):4041-4047.
25. Zhang H, Lund O, Nielsen M. The PickPocket method for predicting binding specificities for receptors based on receptor pocket similarities: application to MHC-peptide binding. *Bioinformatics.* 2009;25(10):1293-1299.
26. Dunbar J, Deane CM. ANARCI: antigen receptor numbering and receptor classification. *Bioinformatics.* 2016;32(2):298-300.

---

*Manuscript prepared: 2026-09-22. For submission to a peer-reviewed journal in computational biology / bioinformatics.*
