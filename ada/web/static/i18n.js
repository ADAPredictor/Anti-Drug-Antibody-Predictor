"use strict";

/* ============================================================
 * i18n: bilingual (zh / en) string table + switching logic
 * ============================================================ */

const I18N = {
  zh: {
    /* ---- header / nav ---- */
    doc_title: "ADA Predictor · 蛋白药物免疫原性预测平台",
    subtitle: "蛋白药物免疫原性预测 · Anti-Drug Antibody Predictor",
    tab_predict: "免疫原性预测",
    tab_drugs: "上市药物参照库",
    tab_methods: "方法说明",
    badge_loading: "后端检测中…",
    badge_offline: "后端不可用",
    /* ---- predict form ---- */
    form_title: "输入蛋白序列",
    form_hint: "支持纯序列或 FASTA；多条链（重链/轻链）请使用多个 > 标题。MHC 等位基因内置，无需手动指定。",
    seq_name_label: "序列名称（可选）",
    seq_name_ph: "如 adalimumab（留空则用 FASTA 标题）",
    seq_ph: "粘贴氨基酸序列或 FASTA（支持多链，如 >adalimumab heavy chain 换行 VH序列 换行 >adalimumab light chain 换行 VL序列）",
    drug_class_label: "药物类别（提示，留空自动推断）",
    run_btn: "运行预测",
    predicting: "预测中…",
    predict_fail: "预测失败",
    /* ---- placeholder ---- */
    placeholder_title: "预测结果将在此显示",
    placeholder_desc: "综合 MHC-I/II 表位负荷、B 细胞表位、CDR/种系异源性与聚集倾向，输出 0–100% 的<b>序列内在免疫原性风险</b>；上市药物的临床 ADA 发生率仅作背景对照，不参与评分。",
    /* ---- results: section titles ---- */
    seq_risk: "序列风险",
    composition: "评分构成",
    raw_stage: "① 原始表位信号加权（校准前）",
    weighted_total: "加权合计",
    ab_calib: "② 抗体 CDR/种系校准",
    ctx_tolerance: "③ 用药情境耐受（基于作用机制，非临床 ADA）",
    no_ctx: "（无情境耐受校正）",
    snapshot: "风险速览",
    nab_risk: "中和抗体（NAb）风险",
    nab_level: "中和等级",
    clinical_ref: "上市药物序列参照",
    mhc1_strong: "MHC-I（CD8）强结合表位",
    mhc2_strong: "MHC-II（CD4）强结合表位",
    bcell_linear: "B 细胞线性表位",
    ab_cdr: "抗体 CDR 区域注释",
    context_notes: "情境注释",
    context_notes_desc: "序列分与真实临床 ADA 之间常见偏差来源的自动注释（仅解释，不改变评分）。",
    ada_rules: "ADA 分级规则",
    full_explain: "完整技术解释",
    items_suffix: " 项，展开）",
    disclaimer: "免责声明",
    /* ---- composition labels ---- */
    comp_weight: "权重",
    calib_baseline: "校准基线",
    multi_uplift: "多信号上浮",
    calib_value: "校准值",
    seq_risk_eq: "序列风险",
    base_explain: "基线 = 序列异源度 × CDR 加权表位负担 × 谱系耐受系数；多信号（CDR 非人源种系表位、MHC/B 覆盖错配残基等）交叉印证后有限上浮（封顶 +15%）。",
    ctx_cap: "，封顶 ",
    clinical_note: "临床 ADA 发生率<b>不参与</b>上述任何一步计算，仅在下方参照表中作背景对照。",
    oncology_note: "例如肿瘤向全抗、CD19/CD20 耗竭等机制会压低<b>实际测得</b>的 ADA，但不改变分子的内在表位负荷。",
    /* ---- protein path ---- */
    protein_calib: "② 蛋白药类别校准",
    lbl_protein_class: "蛋白类别",
    lbl_foreignness: "序列异源度",
    lbl_raw_epitope: "原始表位负荷",
    lbl_tol_mult: "类别耐受系数",
    lbl_gate: "门控系数",
    lbl_calib_seq_risk: "校准后序列风险",
    /* ---- heuristic path ---- */
    tolerance: "② 耐受校正",
    lbl_calib_type: "校正类型",
    lbl_seq_risk: "序列风险",
    tol_endogenous: "近内源耐受 ×0.35",
    tol_humanlike: "人源抗体后备耐受 ×0.40",
    tol_none: "无额外耐受校正",
    /* ---- snapshot fields ---- */
    lbl_path: "评分路径",
    lbl_lineage: "推断谱系",
    lbl_germline_hit: "种系最佳命中",
    lbl_fr_mismatch: "框架区(FR)错配率",
    lbl_cdr_mismatch: "CDR 错配率",
    lbl_ng_burden: "非种系负担",
    lbl_mhc2_strong: "MHC-II（CD4）强表位",
    lbl_mhc1_strong: "MHC-I（CD8）强表位",
    lbl_bcell_cov: "B 细胞表位覆盖率",
    lbl_multi_risk: "多信号风险",
    no_uplift: "未触发多信号上浮",
    ep_count: "个",
    pieces: "段",
    per_100aa: "/100aa",
    heavy_chain: "重链 ",
    light_chain: "轻链 ",
    /* ---- NAb fields ---- */
    lbl_structure: "结构/类别因子",
    lbl_ada_incidence: "邻近药 ADA 发生率因子",
    lbl_literature: "文献/注释因子",
    nab_explain: "中和等级解释（展开）",
    /* ---- neighbors ---- */
    high_conf_note: "命中高置信同药/近同药参照，下列 ADA/NAb 为其<b>临床报告值</b>，仅供背景对照，<b>不参与评分</b>。",
    low_conf_warn: "⚠ 无高置信同药命中——下表仅为框架/局部序列近似的药物（人源抗体框架彼此高度相似属常见现象），其临床 ADA 与输入序列无直接对应关系，<b>不参与评分</b>。",
    th_drug: "药物",
    th_inn: "INN",
    th_class: "类别",
    th_similarity: "相似度",
    th_clinical_ada: "临床 ADA(%)",
    th_clinical_nab: "临床 NAb(%)",
    th_trust: "匹配可信度",
    trust_high_fv: "高置信·Fv 同源",
    trust_fv_near: "Fv 近似",
    trust_vh_only: "⚠ 仅重链框架相似",
    trust_vh_title: "仅重链框架相似，轻链/CDR可能完全不同",
    trust_vl_only: "⚠ 仅轻链框架相似",
    trust_vl_title: "仅轻链框架相似，重链/CDR可能完全不同",
    trust_partial: "⚠ 局部片段相似",
    ref_note: "ADA/NAb 数值来自 IDC-DS-V1 等公开来源的患者加权汇总，受检测方法、适应症、给药途径与合并用药影响，跨药物比较需谨慎。",
    /* ---- epitope tables ---- */
    no_mhc1: "未检出强结合表位。",
    no_mhc2: "未检出强结合表位。",
    no_bcell: "未检出达到阈值的 B 细胞表位段。",
    th_start: "起始",
    th_peptide: "肽段",
    th_allele: "等位基因",
    th_score: "分值",
    th_percentile: "百分位",
    th_15mer: "15肽",
    th_core: "9肽核心",
    th_pos: "位置",
    th_seq: "序列",
    th_method: "方法",
    avg_antigenicity: "平均抗原性",
    coverage: "覆盖率",
    /* ---- antibody regions ---- */
    vh_title: "重链可变区 (VH)",
    vl_title: "轻链可变区 (VL)",
    offset: "偏移",
    length: "长度",
    cdr_note: "CDR 由保守 Cys/Trp 与 WGXG/FGXG 锚点启发式划分（类 Chothia），非 IMGT 正式编号；正式申报请用 IMGT/ANARCI 复核。",
    /* ---- drugs page ---- */
    drugs_title: "上市蛋白药物参照库（289 条）",
    search_ph: "搜索名称 / 通用名(INN)",
    all_categories: "全部类别",
    drugs_hint: "ADA / NAb 发生率来自 IDC-DS-V1 等公开临床数据集（CC BY 4.0），受检测方法、给药途径与合并用药影响，跨药物比较需谨慎。序列与临床数据均不参与预测评分。",
    th_name: "名称",
    th_ada: "ADA(%)",
    th_nab: "NAb(%)",
    th_source: "来源",
    lbl_inn: "通用名(INN)",
    lbl_category: "类别",
    lbl_ada_rate: "ADA 发生率",
    lbl_nab_rate: "中和抗体",
    lbl_nab_note: "NAb 注释",
    lbl_route: "给药途径",
    lbl_approved: "批准信息",
    lbl_data_source: "数据来源",
    lbl_ada_citation: "ADA 文献来源",
    lbl_seq_citation: "序列来源",
    lbl_notes: "备注",
    sequence: "序列",
    aa: " aa）",
    /* ---- level badges ---- */
    lvl_low: "低",
    lvl_moderate: "中",
    lvl_high: "高",
    lvl_unknown: "未知",
    /* ---- footer ---- */
    footer: "ADA Predictor v0.1 · 仅供科研与早期药物筛选使用，不构成临床或注册决策依据",
    footer_github: "源代码 (GitHub)",
    /* ---- hero meta ---- */
    hero_seq: "序列：",
    hero_len: "长度 ",
    hero_path: "评分路径：",
    hero_backend: "MHC 后端：",
  },

  en: {
    /* ---- header / nav ---- */
    doc_title: "ADA Predictor · Protein Drug Immunogenicity Predictor",
    subtitle: "Anti-Drug Antibody Predictor",
    tab_predict: "Immunogenicity Prediction",
    tab_drugs: "Marketed Drug Reference",
    tab_methods: "Methods",
    badge_loading: "Checking backend…",
    badge_offline: "Backend unavailable",
    /* ---- predict form ---- */
    form_title: "Input Protein Sequence",
    form_hint: "Plain sequence or FASTA; for multi-chain (heavy/light) use multiple > headers. MHC alleles are built-in — no manual selection needed.",
    seq_name_label: "Sequence name (optional)",
    seq_name_ph: "e.g. adalimumab (defaults to FASTA header)",
    seq_ph: "Paste amino-acid sequence or FASTA (multi-chain supported, e.g. >adalimumab heavy chain\\n VH sequence\\n >adalimumab light chain\\n VL sequence)",
    drug_class_label: "Drug class (hint; auto-inferred if blank)",
    run_btn: "Run Prediction",
    predicting: "Predicting…",
    predict_fail: "Prediction failed",
    /* ---- placeholder ---- */
    placeholder_title: "Prediction results will appear here",
    placeholder_desc: "Combines MHC-I/II epitope load, B-cell epitopes, CDR/germline foreignness and aggregation propensity into a 0–100% <b>sequence-intrinsic immunogenicity risk</b>; marketed-drug clinical ADA rates are shown only as background reference and do not participate in scoring.",
    /* ---- results: section titles ---- */
    seq_risk: "Sequence Risk",
    composition: "Score Composition",
    raw_stage: "① Raw epitope signal weighting (pre-calibration)",
    weighted_total: "Weighted total",
    ab_calib: "② Antibody CDR / germline calibration",
    ctx_tolerance: "③ Use-context tolerance (mechanism-based, not clinical ADA)",
    no_ctx: "(no context-tolerance correction)",
    snapshot: "Risk Snapshot",
    nab_risk: "Neutralizing Antibody (NAb) Risk",
    nab_level: "NAb level",
    clinical_ref: "Marketed Drug Sequence Reference",
    mhc1_strong: "MHC-I (CD8) strong binders",
    mhc2_strong: "MHC-II (CD4) strong binders",
    bcell_linear: "B-cell linear epitopes",
    ab_cdr: "Antibody CDR region annotation",
    context_notes: "Context notes",
    context_notes_desc: "Automated notes on common sources of divergence between sequence score and real-world clinical ADA (explanatory only, no score change).",
    ada_rules: "ADA classification rules",
    full_explain: "Full technical explanation",
    items_suffix: " items, expand)",
    disclaimer: "Disclaimer",
    /* ---- composition labels ---- */
    comp_weight: "weight",
    calib_baseline: "Calibration baseline",
    multi_uplift: "Multi-signal uplift",
    calib_value: "Calibrated value",
    seq_risk_eq: "Sequence risk",
    base_explain: "Baseline = sequence foreignness × CDR-weighted epitope burden × lineage tolerance factor; multi-signal (CDR non-human-germline epitopes, MHC/B coverage of mismatch residues, etc.) cross-confirmation gives a capped uplift (+15% max).",
    ctx_cap: ", cap ",
    clinical_note: "Clinical ADA incidence <b>does not participate</b> in any of the above steps; it appears only as background reference in the table below.",
    oncology_note: "For example, oncology mAbs and CD19/CD20 depleting mechanisms lower <b>actually measured</b> ADA without changing the molecule's intrinsic epitope burden.",
    /* ---- protein path ---- */
    protein_calib: "② Protein-drug class calibration",
    lbl_protein_class: "Protein class",
    lbl_foreignness: "Sequence foreignness",
    lbl_raw_epitope: "Raw epitope load",
    lbl_tol_mult: "Class tolerance factor",
    lbl_gate: "Gating factor",
    lbl_calib_seq_risk: "Calibrated sequence risk",
    /* ---- heuristic path ---- */
    tolerance: "② Tolerance correction",
    lbl_calib_type: "Correction type",
    lbl_seq_risk: "Sequence risk",
    tol_endogenous: "Near-endogenous tolerance ×0.35",
    tol_humanlike: "Human-like Ab backup tolerance ×0.40",
    tol_none: "No extra tolerance correction",
    /* ---- snapshot fields ---- */
    lbl_path: "Scoring path",
    lbl_lineage: "Inferred lineage",
    lbl_germline_hit: "Germline best hit",
    lbl_fr_mismatch: "FR mismatch rate",
    lbl_cdr_mismatch: "CDR mismatch rate",
    lbl_ng_burden: "Non-germline burden",
    lbl_mhc2_strong: "MHC-II (CD4) strong epitopes",
    lbl_mhc1_strong: "MHC-I (CD8) strong epitopes",
    lbl_bcell_cov: "B-cell epitope coverage",
    lbl_multi_risk: "Multi-signal risk",
    no_uplift: "No multi-signal uplift triggered",
    ep_count: " ",
    pieces: " ",
    per_100aa: "/100aa",
    heavy_chain: "Heavy ",
    light_chain: "Light ",
    /* ---- NAb fields ---- */
    lbl_structure: "Structure / class factor",
    lbl_ada_incidence: "Neighbor-drug ADA incidence factor",
    lbl_literature: "Literature / annotation factor",
    nab_explain: "NAb level explanation (expand)",
    /* ---- neighbors ---- */
    high_conf_note: "High-confidence same/near-same drug hit; the ADA/NAb values below are its <b>clinical reported rates</b>, for background reference only — <b>not used in scoring</b>.",
    low_conf_warn: "⚠ No high-confidence same-drug hit — the table below shows framework/local-sequence-similar drugs only (human Ab frameworks are commonly >85% similar); their clinical ADA has no direct correspondence to the input sequence and <b>does not participate in scoring</b>.",
    th_drug: "Drug",
    th_inn: "INN",
    th_class: "Class",
    th_similarity: "Similarity",
    th_clinical_ada: "Clinical ADA(%)",
    th_clinical_nab: "Clinical NAb(%)",
    th_trust: "Match confidence",
    trust_high_fv: "High confidence · Fv homolog",
    trust_fv_near: "Fv near-match",
    trust_vh_only: "⚠ VH-only framework similarity",
    trust_vh_title: "VH framework only; light chain/CDR may differ entirely",
    trust_vl_only: "⚠ VL-only framework similarity",
    trust_vl_title: "VL framework only; heavy chain/CDR may differ entirely",
    trust_partial: "⚠ Local fragment similarity",
    ref_note: "ADA/NAb values are patient-weighted summaries from IDC-DS-V1 and other public sources; they are affected by assay, indication, route and concomitant medication — cross-drug comparison requires caution.",
    /* ---- epitope tables ---- */
    no_mhc1: "No strong-binding epitope detected.",
    no_mhc2: "No strong-binding epitope detected.",
    no_bcell: "No B-cell epitope segment reaching threshold.",
    th_start: "Start",
    th_peptide: "Peptide",
    th_allele: "Allele",
    th_score: "Score",
    th_percentile: "Percentile",
    th_15mer: "15-mer",
    th_core: "9-mer core",
    th_pos: "Position",
    th_seq: "Sequence",
    th_method: "Method",
    avg_antigenicity: "mean antigenicity",
    coverage: "coverage",
    /* ---- antibody regions ---- */
    vh_title: "Heavy-chain variable region (VH)",
    vl_title: "Light-chain variable region (VL)",
    offset: "offset",
    length: "length",
    cdr_note: "CDRs are delimited by conserved Cys/Trp and WGXG/FGXG anchor heuristics (Chothia-like), not IMGT formal numbering; for formal filings verify with IMGT/ANARCI.",
    /* ---- drugs page ---- */
    drugs_title: "Marketed Protein Drug Reference (289 entries)",
    search_ph: "Search name / INN",
    all_categories: "All classes",
    drugs_hint: "ADA / NAb rates are from IDC-DS-V1 and other public clinical datasets (CC BY 4.0); they are affected by assay, route and concomitant medication — cross-drug comparison requires caution. Sequence and clinical data do not participate in the prediction score.",
    th_name: "Name",
    th_ada: "ADA(%)",
    th_nab: "NAb(%)",
    th_source: "Source",
    lbl_inn: "INN",
    lbl_category: "Class",
    lbl_ada_rate: "ADA incidence",
    lbl_nab_rate: "Neutralizing Ab",
    lbl_nab_note: "NAb note",
    lbl_route: "Route",
    lbl_approved: "Approval",
    lbl_data_source: "Data source",
    lbl_ada_citation: "ADA citation",
    lbl_seq_citation: "Sequence citation",
    lbl_notes: "Notes",
    sequence: "Sequence",
    aa: " aa)",
    /* ---- level badges ---- */
    lvl_low: "Low",
    lvl_moderate: "Mod",
    lvl_high: "High",
    lvl_unknown: "Unk",
    /* ---- footer ---- */
    footer: "ADA Predictor v0.1 · For research and early drug screening only — not a basis for clinical or regulatory decisions",
    footer_github: "Source code (GitHub)",
    /* ---- hero meta ---- */
    hero_seq: "Sequence: ",
    hero_len: "length ",
    hero_path: "Scoring path: ",
    hero_backend: "MHC backend: ",
  }
};

/* Dynamic label maps (language-dependent) */
const DRUG_CLASSES_I18N = {
  zh: {
    "": "自动推断",
    antibody: "抗体（mAb / 片段 / 融合）",
    enzyme: "酶（替代疗法 / 溶栓等）",
    hormone_peptide: "激素 / 多肽",
    cytokine: "细胞因子 / 生长因子",
    fc_fusion: "Fc 融合蛋白",
    conjugate: "偶联药物（ADC 等）",
    other: "其他 / 未知"
  },
  en: {
    "": "Auto-infer",
    antibody: "Antibody (mAb / fragment / fusion)",
    enzyme: "Enzyme (replacement / thrombolytic etc.)",
    hormone_peptide: "Hormone / peptide",
    cytokine: "Cytokine / growth factor",
    fc_fusion: "Fc fusion protein",
    conjugate: "Conjugate (ADC etc.)",
    other: "Other / unknown"
  }
};

const LINEAGE_I18N = {
  zh: {
    chimeric: "嵌合（鼠源可变区）",
    humanized: "人源化",
    human: "全人源",
    murine: "鼠源",
    unknown: "未知"
  },
  en: {
    chimeric: "Chimeric (murine variable)",
    humanized: "Humanized",
    human: "Fully human",
    murine: "Murine",
    unknown: "Unknown"
  }
};

const PATH_I18N = {
  zh: {
    antibody: "抗体 CDR/种系校准路径",
    protein: "蛋白药类别校准路径",
    heuristic: "全长启发式路径"
  },
  en: {
    antibody: "Antibody CDR/germline calibration path",
    protein: "Protein-drug class calibration path",
    heuristic: "Full-length heuristic path"
  }
};

const TOL_REASON_I18N = {
  zh: {
    endogenous: "近内源耐受 ×0.35",
    humanlike: "人源抗体后备耐受 ×0.40",
    none: "无额外耐受校正"
  },
  en: {
    endogenous: "Near-endogenous tolerance ×0.35",
    humanlike: "Human-like Ab backup tolerance ×0.40",
    none: "No extra tolerance correction"
  }
};

/* ---- state ---- */
let currentLang = "zh";

/* ---- helpers ---- */
function t(key) {
  return I18N[currentLang][key] || I18N.zh[key] || key;
}

/* Pick zh or en field from a backend bilingual object */
function bi(obj, zhKey, enKey) {
  if (!obj) return "";
  return currentLang === "en" ? (obj[enKey] || obj[zhKey] || "") : (obj[zhKey] || obj[enKey] || "");
}

/* Pick zh or en from a flat dict with "xxx_zh"/"xxx_en" suffixes */
function biFlat(obj, base) {
  if (!obj) return "";
  const en = obj[base + "_en"];
  const zh = obj[base + "_zh"];
  return currentLang === "en" ? (en || zh || "") : (zh || en || "");
}

/* ---- apply to static [data-i18n] elements ---- */
function applyI18n() {
  document.documentElement.lang = currentLang === "en" ? "en" : "zh-CN";
  document.title = t("doc_title");
  document.querySelectorAll("[data-i18n]").forEach(el => {
    const key = el.getAttribute("data-i18n");
    const val = t(key);
    if (val.includes("<")) el.innerHTML = val; else el.textContent = val;
  });
  document.querySelectorAll("[data-i18n-ph]").forEach(el => {
    el.placeholder = t(el.getAttribute("data-i18n-ph"));
  });
  const toggle = document.getElementById("lang-toggle");
  if (toggle) toggle.textContent = currentLang === "zh" ? "EN" : "中文";
  document.dispatchEvent(new CustomEvent("langchange", {detail: currentLang}));
}

function setLang(lang) {
  if (lang === currentLang) return;
  currentLang = lang;
  applyI18n();
}

document.addEventListener("DOMContentLoaded", () => {
  const toggle = document.getElementById("lang-toggle");
  if (toggle) {
    toggle.addEventListener("click", () => {
      setLang(currentLang === "zh" ? "en" : "zh");
    });
  }
});
