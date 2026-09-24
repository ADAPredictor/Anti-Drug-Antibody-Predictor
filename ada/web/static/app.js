"use strict";

/* immunopred drug-class options (passed as hint; auto-inferred when blank) */
function drugClasses() { return DRUG_CLASSES_I18N[currentLang]; }
function lineageLabel(k) { return (LINEAGE_I18N[currentLang][k] || k); }
function pathLabel(k) { return (PATH_I18N[currentLang][k] || k); }
function tolReason(k) { return (TOL_REASON_I18N[currentLang][k] || k); }
/* pick zh or en field from a backend bilingual item dict */
function itemText(item) { return currentLang === "en" ? (item.en || item.zh || "") : (item.zh || item.en || ""); }
function itemNote(item) { return currentLang === "en" ? (item.note_en || item.note_zh || "") : (item.note_zh || item.note_en || ""); }

/* ------------------------------- tabs ------------------------------- */
document.querySelectorAll(".tab").forEach(btn => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".tab").forEach(b => b.classList.remove("active"));
    document.querySelectorAll(".tab-panel").forEach(p => p.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById("tab-" + btn.dataset.tab).classList.add("active");
    if (btn.dataset.tab === "drugs") loadDrugs();
    if (btn.dataset.tab === "methods") renderMethods();
  });
});

/* re-render dynamic content on language change */
document.addEventListener("langchange", () => {
  rebuildFormOptions();
  renderMethods();
  if (!document.getElementById("results").querySelector(".placeholder-card")) {
    if (window._lastResult) renderResult(window._lastResult);
  }
  if (document.getElementById("tab-drugs").classList.contains("active")) loadDrugs();
});

/* ------------------------------ meta -------------------------------- */
async function loadMeta() {
  const badge = document.getElementById("backend-badge");
  try {
    const m = await fetch("/api/meta").then(r => r.json());
    const ref = m.reference_db || {};
    if (badge) {
      if (currentLang === "en") {
        badge.textContent = "Engine: " + (m.engine || "immunopred")
          + " · Reference " + (ref.n_entries || 0) + " entries"
          + " (MHC-I " + (m.default_mhc1_alleles || []).length
          + " / MHC-II " + (m.default_mhc2_alleles || []).length + " alleles)";
      } else {
        badge.textContent = "引擎：" + (m.engine || "immunopred")
          + " · 参照库 " + (ref.n_entries || 0) + " 条"
          + "（MHC-I " + (m.default_mhc1_alleles || []).length
          + " / MHC-II " + (m.default_mhc2_alleles || []).length + " 等位）";
      }
      badge.className = "badge";
    }
    // Populate drug class filter from API
    const dsel = document.getElementById("drug-modality");
    const curVal = dsel.value;
    dsel.innerHTML = '<option value="">' + t("all_categories") + '</option>';
    (m.drug_classes || []).forEach(cls => {
      const o = document.createElement("option");
      o.value = cls; o.textContent = cls;
      dsel.appendChild(o);
    });
    dsel.value = curVal;
    loadDrugs();
  } catch (e) {
    badge.textContent = t("badge_offline");
    badge.className = "badge badge-muted";
  }
}

/* ----------------------------- predict ------------------------------ */
function rebuildFormOptions() {
  const sel = document.getElementById("drug-class");
  const cur = sel.value;
  sel.innerHTML = "";
  Object.entries(drugClasses()).forEach(([k, v]) => {
    const o = document.createElement("option");
    o.value = k; o.textContent = v;
    sel.appendChild(o);
  });
  sel.value = cur || "antibody";
}

function initForm() {
  rebuildFormOptions();
  document.getElementById("run-btn").addEventListener("click", runPrediction);
}

async function runPrediction() {
  const btn = document.getElementById("run-btn");
  const errBox = document.getElementById("predict-error");
  errBox.hidden = true;
  btn.disabled = true; btn.textContent = t("predicting");
  const body = {
    sequence: document.getElementById("seq-input").value,
    name: document.getElementById("seq-name").value.trim() || null,
    drug_class: document.getElementById("drug-class").value
  };
  try {
    const resp = await fetch("/api/predict", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify(body)
    });
    const data = await resp.json();
    if (!resp.ok) throw new Error(data.error || t("predict_fail"));
    window._lastResult = data;
    renderResult(data);
  } catch (e) {
    errBox.textContent = e.message;
    errBox.hidden = false;
  } finally {
    btn.disabled = false; btn.textContent = t("run_btn");
  }
}

/* --------------------------- result render --------------------------- */
function tierColor(tier) {
  return {low: "#2e9e5b", moderate: "#d99317", high: "#d34c4c", unknown: "#888"}[tier] || "#888";
}

function fmtNum(v, digits) {
  if (v == null || isNaN(v)) return "—";
  return Number(v).toFixed(digits == null ? 2 : digits);
}

/* All risk-scale values are 0-100 and displayed as percentages */
function pct(v, digits) {
  if (v == null || isNaN(v)) return "—";
  return Number(v).toFixed(digits == null ? 1 : digits) + "%";
}

function ratioPct(v, digits) {
  if (v == null || isNaN(v)) return "—";
  return (v * 100).toFixed(digits == null ? 1 : digits) + "%";
}

function esc(s) {
  if (s == null) return "";
  return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

function kvRow(label, val) {
  return "<tr><td>" + esc(label) + "</td><td>" + (val == null ? "—" : val) + "</td></tr>";
}

function renderResult(d) {
  const r = d.risk || {};
  const b = r.breakdown || {};
  const nab = d.nab || {};
  const feat = d.features || {};
  const col = tierColor(r.tier);
  const ncol = tierColor(nab.nab_level);
  let html = "";

  /* ---- hero (percent gauge) ---- */
  html += '<div class="card"><div class="score-hero">'
    + '<div class="gauge" style="--p:' + (r.score || 0) + ';--gauge-color:' + col + '">'
    + '<div class="gauge-inner"><div class="gauge-num">' + pct(r.score, 1)
    + '</div><div class="gauge-label">' + esc(t("seq_risk")) + '</div></div></div>'
    + '<div><div class="tier-name" style="color:' + col + '">'
    + esc(currentLang === "en" ? (r.ada_level_label_en || r.ada_level_label || r.tier || "") : (r.ada_level_label || r.tier || "")) + '</div>'
    + '<div class="hero-meta">' + esc(t("hero_seq")) + esc(d.sequence_id || "—")
    + " · " + esc(t("hero_len")) + (d.length || 0) + " aa"
    + "<br>" + esc(t("hero_path")) + esc(pathLabel(b.path) || b.path || "—")
    + "<br>" + esc(t("hero_backend")) + esc(((d.mhc_backend || {}).used) || "heuristic")
    + '</div></div></div></div>';

  /* ---- score composition bars ---- */
  html += renderComposition(b);

  /* ---- concise risk snapshot ---- */
  html += renderSnapshot(b, feat, d);

  /* ---- NAb ---- */
  html += renderNab(nab, ncol);

  /* ---- neighbors (clinical reference only, NOT scored) ---- */
  html += renderNeighbors(feat);

  /* ---- epitope tables ---- */
  html += renderEpitopeTables(d);

  /* ---- antibody CDR / germline detail ---- */
  const ar = feat.antibody_regions;
  if (ar && ar.is_antibody_like) html += renderAntibodyRegions(feat);

  /* ---- context notes ---- */
  if ((d.context_notes || []).length) html += renderContextNotes(d.context_notes);

  /* ---- level rules ---- */
  if ((r.level_rules_bi || []).length) {
    html += '<div class="card"><h2>' + esc(t("ada_rules")) + '</h2><ul class="rules-list">';
    r.level_rules_bi.forEach(rule => {
      html += "<li>" + esc(itemText(rule)) + "</li>";
    });
    html += "</ul></div>";
  }

  /* ---- full technical explanation (collapsed) ---- */
  if ((r.explanation_items || []).length) {
    html += '<div class="card"><details><summary class="details-summary">' + esc(t("full_explain"))
      + "（" + r.explanation_items.length + esc(t("items_suffix")) + "</summary>";
    r.explanation_items.forEach(item => {
      html += '<div class="explain-item"><div class="explain-zh">' + esc(itemText(item)) + "</div>"
        + (itemNote(item) ? '<div class="hint">' + esc(itemNote(item)) + "</div>" : "")
        + "</div>";
    });
    html += "</details></div>";
  }

  /* ---- disclaimer ---- */
  html += '<div class="card" style="border-left:4px solid #b58a2a"><h2>' + esc(t("disclaimer")) + "</h2>"
    + '<p class="hint">' + esc(currentLang === "en" ? (d.disclaimer_en || d.disclaimer_zh || "") : (d.disclaimer_zh || d.disclaimer_en || "")) + "</p></div>";

  document.getElementById("results").innerHTML = html;
}

/* ----------------------- composition bar chart ----------------------- */
function weightBar(label, weight, term) {
  return '<div class="comp-row"><div class="comp-label">' + esc(label)
    + ' <span class="hint">(' + esc(t("comp_weight")) + " " + Math.round(weight * 100) + "%)</span></div>"
    + '<div class="comp-bar"><div class="comp-fill" style="width:' + term + '%"></div></div>'
    + '<div class="comp-val">' + pct(term, 1) + "</div></div>";
}

function rawTermLabel(rt) {
  if (currentLang === "en") return rt.label_en || rt.label_zh || rt.key;
  return rt.label_zh || rt.label_en || rt.key;
}

function renderComposition(b) {
  let h = '<div class="card"><h2>' + esc(t("composition")) + '</h2>';

  /* Stage 1: weighted raw epitope signals */
  h += '<h3>' + esc(t("raw_stage")) + '</h3>';
  (b.raw_terms || []).forEach(rt => { h += weightBar(rawTermLabel(rt), rt.weight, rt.term); });
  h += '<div class="comp-row comp-total"><div class="comp-label"><b>' + esc(t("weighted_total")) + '</b></div>'
    + '<div class="comp-bar"><div class="comp-fill comp-fill-total" style="width:'
    + b.full_seq_epitope + '%"></div></div>'
    + '<div class="comp-val"><b>' + pct(b.full_seq_epitope, 1) + "</b></div></div>";

  /* Stage 2: path-specific calibration */
  if (b.path === "antibody") {
    const base = b.base_score || 0;
    const up = b.risk_uplift || 0;
    const preCtx = b.pre_context_score || 0;
    const final = b.epitope_component || 0;
    const ctxMult = b.context_mult != null ? b.context_mult : 1;
    const hasCtx = ctxMult < 1.0;
    h += '<h3 style="margin-top:16px">' + esc(t("ab_calib")) + '</h3>'
      + '<div class="stack-bar">'
      + '<div class="stack-base" style="width:' + base + '%"></div>'
      + '<div class="stack-up" style="width:' + up + '%"></div>'
      + "</div>"
      + '<div class="stack-legend"><span><i class="sw sw-base"></i>' + esc(t("calib_baseline")) + " " + pct(base, 1)
      + "</span><span><i class='sw sw-up'></i>" + esc(t("multi_uplift")) + " +" + pct(up, 1)
      + "</span><span class='stack-final'>= " + esc(t("calib_value")) + " <b>" + pct(preCtx, 1)
      + "</b></span></div>"
      + '<p class="hint">' + esc(t("base_explain")) + "</p>";
    if (hasCtx) {
      const capTxt = b.context_cap != null ? esc(t("ctx_cap")) + pct(b.context_cap, 0) : "";
      const ctxNote = currentLang === "en" ? (b.context_note_en || b.context_note_zh || "") : (b.context_note_zh || b.context_note_en || "");
      h += '<h3 style="margin-top:14px">' + esc(t("ctx_tolerance")) + '</h3>'
        + '<div class="ctx-calc">' + pct(preCtx, 1) + " × " + fmtNum(ctxMult, 2)
        + " → " + esc(t("seq_risk_eq")) + " <b style='font-size:15px'>" + pct(final, 1) + "</b>"
        + (capTxt ? " <span class='hint'>" + capTxt + "</span>" : "")
        + "</div>"
        + (ctxNote ? '<div class="hint">' + esc(ctxNote) + "</div>" : "")
        + '<p class="hint">' + esc(t("oncology_note")) + "</p>";
    } else {
      h += '<p style="margin-top:8px"><b>' + esc(t("seq_risk_eq")) + " = " + pct(final, 1)
        + "</b>" + esc(t("no_ctx")) + "</p>";
    }
  } else if (b.path === "protein") {
    h += '<h3 style="margin-top:16px">' + esc(t("protein_calib")) + '</h3><table class="kv">'
      + kvRow(t("lbl_protein_class"), esc(b.protein_class))
      + kvRow(t("lbl_foreignness"), pct(b.foreignness, 0))
      + kvRow(t("lbl_raw_epitope"), pct(b.epitope_raw, 1))
      + kvRow(t("lbl_tol_mult"), fmtNum(b.tolerance_mult, 2))
      + kvRow(t("lbl_gate"), fmtNum(b.gate, 2))
      + kvRow(t("lbl_calib_seq_risk"), "<b>" + pct(b.epitope_component, 1) + "</b>")
      + "</table>";
  } else {
    h += '<h3 style="margin-top:16px">' + esc(t("tolerance")) + '</h3><table class="kv">'
      + kvRow(t("lbl_calib_type"), esc(tolReason(b.tolerance_reason) || "—"))
      + kvRow(t("lbl_seq_risk"), "<b>" + pct(b.epitope_component, 1) + "</b>")
      + "</table>";
  }
  h += '<p class="hint">' + t("clinical_note") + "</p>";
  h += "</div>";
  return h;
}

/* ------------------------ concise risk snapshot ---------------------- */
function renderSnapshot(b, feat, d) {
  const gr = feat.germline_report || {};
  const hits = gr.hits || [];
  const hitTxt = hits.map(x => (x.chain === "H" ? t("heavy_chain") : t("light_chain"))
    + esc(x.v_gene) + "（" + ratioPct(x.identity, 0) + "）").join("<br>") || "—";
  const m1 = d.mhc1 || {}, m2 = d.mhc2 || {}, bc = d.b_cell || {};
  const flags = currentLang === "en" ? (b.risk_flags_en || b.risk_flags_zh || []) : (b.risk_flags_zh || b.risk_flags_en || []);

  const left = [
    [t("lbl_path"), esc(pathLabel(b.path) || b.path)],
    b.path === "antibody" ? [t("lbl_lineage"), esc(lineageLabel(b.lineage_inferred) || b.lineage_inferred)] : null,
    b.path !== "protein" && b.foreignness != null ? [t("lbl_foreignness"), pct(b.foreignness, 0)] : null,
    b.path === "protein" ? [t("lbl_foreignness"), pct(b.foreignness, 0)] : null,
    gr.fr_hetero_frac != null ? [t("lbl_fr_mismatch"), ratioPct(gr.fr_hetero_frac, 1)] : null,
    gr.cdr_hetero_frac != null ? [t("lbl_cdr_mismatch"), ratioPct(gr.cdr_hetero_frac, 1)] : null,
    feat.non_germline_burden != null ? [t("lbl_ng_burden"), pct(feat.non_germline_burden, 1)] : null,
    hits.length ? [t("lbl_germline_hit"), hitTxt] : null
  ].filter(Boolean).map(([k, v]) => kvRow(k, v)).join("");

  const flagTags = flags.length
    ? flags.map(f => '<span class="tag tag-intermediate">' + esc(f) + "</span>").join(" ")
    : '<span class="tag tag-low">' + esc(t("no_uplift")) + "</span>";

  const right = [
    kvRow(t("lbl_mhc2_strong"), (m2.strong_count || 0) + " " + t("ep_count") + " · "
      + fmtNum(m2.epitope_density, 2) + t("per_100aa")),
    kvRow(t("lbl_mhc1_strong"), (m1.strong_count || 0) + " " + t("ep_count") + " · "
      + fmtNum(m1.epitope_density, 2) + t("per_100aa")),
    kvRow(t("lbl_bcell_cov"), ratioPct(bc.coverage, 1) + "（"
      + ((bc.epitopes || []).length) + " " + t("pieces") + "）"),
    kvRow(t("lbl_multi_risk"), flagTags)
  ].join("");

  return '<div class="card"><h2>' + esc(t("snapshot")) + '</h2><div class="grid-2">'
    + '<table class="kv">' + left + "</table>"
    + '<table class="kv">' + right + "</table>"
    + "</div></div>";
}

/* ------------------------------- NAb --------------------------------- */
function renderNab(nab, ncol) {
  let h = '<div class="card"><h2>' + esc(t("nab_risk")) + "</h2>"
    + '<div class="grid-2"><div>'
    + '<div class="gauge" style="--p:' + (nab.score || 0) + ';--gauge-color:' + ncol
    + ';width:96px;height:96px">'
    + '<div class="gauge-inner" style="width:74px;height:74px">'
    + '<div class="gauge-num" style="font-size:20px">' + pct(nab.score, 1)
    + '</div><div class="gauge-label">NAb</div></div></div>'
    + '<div style="margin-top:8px;font-size:13px;color:' + ncol
    + ';font-weight:700">' + esc(currentLang === "en" ? (nab.nab_level_label_en || nab.nab_level_label || "—") : (nab.nab_level_label || nab.nab_level || "—"))
    + " " + esc(t("nab_level")) + "</div></div>"
    + '<div><table class="kv">'
    + kvRow(t("lbl_structure"), pct(nab.structure_component, 1))
    + kvRow(t("lbl_ada_incidence"), pct(nab.ada_incidence_component, 1))
    + kvRow(t("lbl_literature"), pct(nab.literature_component, 1))
    + "</table></div></div>";
  if ((nab.explanation_items || []).length) {
    h += '<details><summary class="details-summary">' + esc(t("nab_explain")) + "</summary>";
    nab.explanation_items.forEach(item => {
      h += '<div class="explain-item"><div class="explain-zh">' + esc(itemText(item)) + "</div>"
        + (itemNote(item) ? '<div class="hint">' + esc(itemNote(item)) + "</div>" : "")
        + "</div>";
    });
    h += "</details>";
  }
  h += "</div>";
  return h;
}

/* ---------------------- neighbors (reference only) ------------------- */
function neighborTrustTag(n) {
  const sim = n.similarity || 0;
  if (n.method === "fv" && sim >= 0.95)
    return '<span class="tag tag-low">' + esc(t("trust_high_fv")) + "</span>";
  if (n.method === "fv" && sim >= 0.85)
    return '<span class="tag tag-intermediate">' + esc(t("trust_fv_near")) + "</span>";
  if (n.method === "vh")
    return '<span class="tag tag-weak" title="' + esc(t("trust_vh_title")) + '">' + esc(t("trust_vh_only")) + "</span>";
  if (n.method === "vl")
    return '<span class="tag tag-weak" title="' + esc(t("trust_vl_title")) + '">' + esc(t("trust_vl_only")) + "</span>";
  return '<span class="tag tag-weak">' + esc(t("trust_partial")) + "</span>";
}

function renderNeighbors(feat) {
  const neighbors = feat.neighbors || [];
  if (!neighbors.length) return "";
  const best = neighbors[0] || {};
  const highConf = (best.method === "fv" && (best.similarity || 0) >= 0.95);
  let h = '<div class="card" style="border-left:4px solid #2a7f9e">'
    + "<h2>" + esc(t("clinical_ref")) + "</h2>"
    + (highConf
       ? '<p class="hint">' + t("high_conf_note") + "</p>"
       : '<div class="warn">' + t("low_conf_warn") + "</div>")
    + '<div class="table-wrap"><table class="seg-table"><thead><tr>'
    + "<th>" + esc(t("th_drug")) + "</th><th>" + esc(t("th_inn")) + "</th><th>" + esc(t("th_class")) + "</th><th>" + esc(t("th_similarity")) + "</th>"
    + "<th>" + esc(t("th_clinical_ada")) + "</th><th>" + esc(t("th_clinical_nab")) + "</th><th>" + esc(t("th_trust")) + "</th></tr></thead><tbody>";
  neighbors.slice(0, 8).forEach(n => {
    h += "<tr><td><b>" + esc(n.name) + "</b></td><td>" + esc(n.inn || "—") + "</td>"
      + "<td>" + esc(n.drug_class || "—") + "</td>"
      + "<td>" + ratioPct(n.similarity, 1) + "</td>"
      + "<td>" + (n.ada_rate_pct != null ? fmtNum(n.ada_rate_pct, 1) : "—") + "</td>"
      + "<td>" + (n.nab_rate_pct != null ? fmtNum(n.nab_rate_pct, 1) : "—") + "</td>"
      + "<td>" + neighborTrustTag(n) + "</td></tr>";
  });
  h += "</tbody></table></div>"
    + '<p class="hint">' + esc(t("ref_note")) + "</p></div>";
  return h;
}

/* --------------------------- epitope tables -------------------------- */
function renderEpitopeTables(d) {
  const m1 = d.mhc1 || {}, m2 = d.mhc2 || {}, bc = d.b_cell || {};
  let h = '<div class="grid-2">'
    + '<div class="card"><h2>' + esc(t("mhc1_strong")) + " · " + (m1.strong_count || 0)
    + " " + t("ep_count") + " · " + fmtNum(m1.epitope_density, 2) + t("per_100aa") + "</h2>"
    + epitopeTableMhc1(m1.epitopes || []) + "</div>"
    + '<div class="card"><h2>' + esc(t("mhc2_strong")) + " · " + (m2.strong_count || 0)
    + " " + t("ep_count") + " · " + fmtNum(m2.epitope_density, 2) + t("per_100aa") + "</h2>"
    + epitopeTableMhc2(m2.epitopes || []) + "</div>"
    + "</div>";
  h += '<div class="card"><h2>' + esc(t("bcell_linear")) + " · " + esc(t("coverage")) + " "
    + ratioPct(bc.coverage, 1) + " · " + esc(t("avg_antigenicity")) + " " + fmtNum(bc.mean_antigenicity, 2)
    + "</h2>" + bcellTable(bc.epitopes || []) + "</div>";
  return h;
}

function epitopeTableMhc1(eps) {
  if (!eps || !eps.length) return "<p class='hint'>" + esc(t("no_mhc1")) + "</p>";
  const rows = eps.slice(0, 30).map(e => {
    return "<tr><td>" + e.start + '</td><td><code>' + esc(e.peptide || "") + "</code></td>"
      + "<td>" + esc(e.allele || "") + "</td>"
      + "<td>" + fmtNum(e.score, 3) + "</td>"
      + "<td>" + fmtNum(e.percentile, 1) + "</td></tr>";
  }).join("");
  return '<table class="seg-table"><thead><tr><th>' + esc(t("th_start")) + '</th><th>' + esc(t("th_peptide")) + '</th><th>' + esc(t("th_allele")) + '</th><th>'
    + esc(t("th_score")) + '</th><th>' + esc(t("th_percentile")) + "</th></tr></thead><tbody>" + rows + "</tbody></table>";
}

function epitopeTableMhc2(eps) {
  if (!eps || !eps.length) return "<p class='hint'>" + esc(t("no_mhc2")) + "</p>";
  const rows = eps.slice(0, 30).map(e => {
    return "<tr><td>" + e.start + '</td><td><code>' + esc(e.peptide || "") + "</code></td>"
      + '<td><code>' + esc(e.core || "") + "</code></td>"
      + "<td>" + esc(e.allele || "") + "</td>"
      + "<td>" + fmtNum(e.score, 3) + "</td></tr>";
  }).join("");
  return '<table class="seg-table"><thead><tr><th>' + esc(t("th_start")) + '</th><th>' + esc(t("th_15mer")) + '</th><th>' + esc(t("th_core")) + '</th>'
    + '<th>' + esc(t("th_allele")) + '</th><th>' + esc(t("th_score")) + "</th></tr></thead><tbody>" + rows + "</tbody></table>";
}

function bcellTable(eps) {
  if (!eps || !eps.length) return "<p class='hint'>" + esc(t("no_bcell")) + "</p>";
  const rows = eps.slice(0, 20).map(e => {
    return "<tr><td>" + e.start + "-" + e.end + '</td><td><code>' + esc(e.sequence || "")
      + "</code></td><td>" + fmtNum(e.score, 3) + "</td>"
      + "<td>" + esc(e.method || "") + "</td></tr>";
  }).join("");
  return '<table class="seg-table"><thead><tr><th>' + esc(t("th_pos")) + '</th><th>' + esc(t("th_seq")) + '</th><th>' + esc(t("th_score")) + '</th>'
    + '<th>' + esc(t("th_method")) + "</th></tr></thead><tbody>" + rows + "</tbody></table>";
}

/* ------------------------ antibody region detail --------------------- */
function renderAntibodyRegions(feat) {
  const ar = feat.antibody_regions;
  let h = '<div class="card"><h2>' + esc(t("ab_cdr")) + "</h2>";
  (ar.domains || []).forEach(dom => {
    h += "<h3>" + (dom.kind === "vh" ? esc(t("vh_title")) : esc(t("vl_title")))
      + " · " + esc(t("offset")) + " " + dom.offset + " · " + esc(t("length")) + " " + dom.length + "</h3>"
      + '<table class="seg-table"><thead><tr><th>CDR</th><th>' + esc(t("th_pos")) + '</th><th>' + esc(t("th_seq")) + '</th></tr></thead><tbody>';
    (dom.cdrs || []).forEach(cdr => {
      h += "<tr><td><b>" + esc(cdr.name) + "</b></td><td>"
        + (cdr.start + 1) + "-" + cdr.end + '</td><td><code>'
        + esc(cdr.sequence || "") + "</code></td></tr>";
    });
    h += "</tbody></table>";
  });
  h += '<p class="hint">' + esc(t("cdr_note")) + "</p></div>";
  return h;
}

/* ---------------------------- context notes -------------------------- */
function renderContextNotes(notes) {
  let h = '<div class="card"><h2>' + esc(t("context_notes")) + "</h2>"
    + '<p class="hint">' + esc(t("context_notes_desc")) + "</p>";
  notes.forEach(n => {
    const title = currentLang === "en" ? (n.title_en || n.title_zh || "") : (n.title_zh || n.title_en || "");
    const body = currentLang === "en" ? (n.body_en || n.body_zh || "") : (n.body_zh || n.body_en || "");
    const examples = currentLang === "en" ? (n.examples_en || n.examples_zh || "") : (n.examples_zh || n.examples_en || "");
    h += '<div class="ctx-note"><div class="ctx-title">' + esc(title)
      + (n.empiric_delta && n.empiric_delta !== "例外见下"
         ? ' <span class="tag tag-intermediate">' + esc(n.empiric_delta) + "</span>" : "")
      + "</div><div>" + esc(body) + "</div>"
      + (examples ? '<div class="hint">' + esc(examples) + "</div>" : "")
      + "</div>";
  });
  return h + "</div>";
}

/* ------------------------------- drugs ------------------------------ */
function adaBadge(cat) {
  if (!cat) return '<span class="tag tag-weak">' + esc(t("lvl_unknown")) + "</span>";
  const map = {low: ["tag-low", "lvl_low"], moderate: ["tag-moderate", "lvl_moderate"], high: ["tag-high", "lvl_high"]};
  const [cls, labKey] = map[cat] || map.moderate;
  return '<span class="tag ' + cls + '">' + esc(t(labKey)) + "</span>";
}

let allDrugs = [];

async function loadDrugs() {
  const q = document.getElementById("drug-search").value.trim();
  const mod = document.getElementById("drug-modality").value;
  const params = new URLSearchParams();
  if (q) params.set("q", q);
  if (mod) params.set("drug_class", mod);
  allDrugs = await fetch("/api/drugs?" + params).then(r => r.json());
  const tbody = document.querySelector("#drugs-table tbody");
  tbody.innerHTML = allDrugs.map(d => {
    const ada = d.ada_rate_pct != null ? d.ada_rate_pct.toFixed(2) + "%" : "—";
    const nab = d.nab_rate_pct != null ? d.nab_rate_pct.toFixed(1) + "%" : "—";
    return "<tr data-id='" + d.id + "'><td><b>" + esc(d.name) + "</b></td><td>"
      + esc(d.inn || "—") + "</td><td>" + esc(d.drug_class || "—")
      + "</td><td>" + ada + "</td><td>" + nab + "</td><td>" + esc(d.source || "—")
      + "</td></tr>";
  }).join("");
  tbody.querySelectorAll("tr").forEach(tr => {
    tr.addEventListener("click", () => showDrug(tr.dataset.id));
  });
}

async function showDrug(id) {
  const d = await fetch("/api/drugs/" + id).then(r => r.json());
  const box = document.getElementById("drug-detail");
  box.hidden = false;
  let seqHtml = "";
  if (d.heavy || d.light) {
    seqHtml = '<h3>' + esc(t("sequence")) + '</h3>';
    if (d.heavy) seqHtml += '<p><b>VH</b>（' + d.heavy.length + esc(t("aa")) + '</p><pre class="seq-block">' + esc(d.heavy) + '</pre>';
    if (d.light) seqHtml += '<p><b>VL</b>（' + d.light.length + esc(t("aa")) + '</p><pre class="seq-block">' + esc(d.light) + '</pre>';
  }
  const nabNote = currentLang === "en" ? (d.nab_note || d.nab_note_zh || "") : (d.nab_note_zh || d.nab_note || "");
  box.innerHTML = "<h2>" + esc(d.name) + "</h2>"
    + '<table class="kv">'
    + kvRow(t("lbl_inn"), d.inn || "—")
    + kvRow(t("lbl_category"), d.drug_class || "—")
    + kvRow(t("lbl_ada_rate"), d.ada_rate_pct != null ? d.ada_rate_pct.toFixed(2) + "%" : "—")
    + kvRow(t("lbl_nab_rate"), d.nab_rate_pct != null ? d.nab_rate_pct.toFixed(1) + "%" : "—")
    + (nabNote ? kvRow(t("lbl_nab_note"), nabNote) : "")
    + kvRow(t("lbl_route"), d.route || "—")
    + kvRow(t("lbl_approved"), d.approved || "—")
    + kvRow(t("lbl_data_source"), d.source || "—")
    + (d.ada_citation ? kvRow(t("lbl_ada_citation"), d.ada_citation) : "")
    + (d.sequence_citation ? kvRow(t("lbl_seq_citation"), d.sequence_citation) : "")
    + (d.notes ? kvRow(t("lbl_notes"), d.notes) : "")
    + '</table>'
    + seqHtml;
  box.scrollIntoView({behavior: "smooth", block: "nearest"});
}

document.getElementById("drug-search").addEventListener("input", () => loadDrugs());
document.getElementById("drug-modality").addEventListener("change", () => loadDrugs());

/* ------------------------------ methods ----------------------------- */
function renderMethods() {
  document.getElementById("methods-content").innerHTML = currentLang === "en" ? methodsEn() : methodsZh();
}

function methodsZh() {
  return `
<h3>1. 引擎与数据</h3>
<p>预测引擎为 <b>immunopred</b>：对输入蛋白序列做多维度特征提取（MHC-I/II 表位、B 细胞线性表位、
抗体 CDR/框架种系比对、理化聚集倾向、人源同源 foreignness），综合校准为 0–100% 的<b>序列内在免疫原性风险</b>。
参照库为 <code>ada_reference.csv</code>（289 条上市生物药，含 ADA/NAb 发生率、序列、heavy/light），
仅用于结果页的临床背景对照与中和等级的文献因子，<b>不参与序列风险分计算</b>。</p>

<h3>2. 评分构成</h3>
<p><b>第一步 · 原始表位信号加权：</b>MHC-II 表位密度 45% + MHC-II 强结合分值 20% + B 细胞表位覆盖 20%
+ MHC-I 表位密度 8% + 聚集倾向 7%，得到校准前的全长表位负荷。</p>
<p><b>第二步 · 路径校准：</b>抗体走 CDR/种系校准——序列异源度（骨架 + FR/CDR 种系错配 + 格式暴露）
× CDR 加权表位负担 × 谱系耐受系数（全人源 0.55 / 人源化 0.40 / 嵌合 1.15）得到基线，
多信号交叉印证后有限上浮（封顶 +15%）；非抗体蛋白走类别特异 foreignness × 表位 × 类别耐受。</p>
<p>分级：<b>&lt;20% 低 / 20–45% 中 / ≥45% 高</b>。</p>

<h3>3. T 细胞表位</h3>
<ul>
<li><b>MHC-I（CD8）：</b>9 肽滑窗，7 个常见 HLA-A/B 等位基因的口袋偏好矩阵（锚定残基 P2/C 端），
按蛋白内百分位排序，top 10% 且分值≥0.15 标记为强结合。</li>
<li><b>MHC-II（CD4）：</b>15 肽滑窗 + 9 肽核心扫描，6 个 HLA-DRB1 等位基因的 P1/P4/P6/P9 口袋 PSSM，
强结合阈值 = 等位基因 threshold_strong × 0.5 且百分位≥90。</li>
</ul>

<h3>4. B 细胞表位</h3>
<p>线性表位共识分 = Kolaskar 抗原性 + Parker 亲水性 + Emini 表面可及性 + Chou-Fasman β-转角的
z-score 均值；连续 ≥6 aa 且共识分≥0.5（z-score）判为候选表位。
线性量表对表面的区分能力有限（&gt;90% B 细胞表位为构象表位），仅用于候选区筛查。</p>

<h3>5. 抗体校准</h3>
<p>CDR 由保守 Cys/Trp 与 WGXG/FGXG 锚点启发式划分（类 Chothia），
与 IMGT 人源 V 基因库比对得到 FR/CDR 错配分率与种系最佳命中；
lineage 按 INN 后缀（-ximab/-zumab/-umab）与种系 FR 双重推断。
多信号风险交叉验证（CDR 非人源种系表位偏高、MHC-II 肽覆盖种系异源等）可触发表位上浮。</p>

<h3>6. 中和抗体（NAb）风险</h3>
<p>独立评分：结构/类别因子 + 邻近药 ADA 发生率因子 + 文献/注释因子。
低（&lt;25%）/ 中（25–55%）/ 高（≥55%）。预测的是中和可能性倾向，非确认的中和滴度。</p>

<h3>7. 序列参照匹配的可信度</h3>
<p>参照表按链感知比对：Fv（重+轻双链）高相似（≥95%）为高置信同药命中；
仅 VH 或 VL 单链相似（人源抗体通用框架彼此相似度常达 85%+）会明确标注「仅单链框架相似」，
不等同于同一药物；k-mer 局部片段相似仅作提示。临床 ADA 数值在任何情况下都不回灌评分。</p>

<h3>8. 局限与免责</h3>
<ul>
<li>MHC 结合为内置启发式（口袋偏好矩阵），不能替代 NetMHCpan / IEDB / 湿实验。</li>
<li>临床 ADA 受给药途径、频次、合并用药、检测方法、人群 HLA 影响，不可由序列分直接换算。</li>
<li>CDR 划分为启发式锚点，非 IMGT 正式编号；正式申报请用 ANARCI 复核。</li>
<li>本工具仅供科研与早期候选物筛选，<b>不构成临床或注册决策依据</b>。</li>
</ul>

<h3>参考文献</h3>
<ol class="refs">
<li>Kolaskar AS, Tongaonkar PC. A semi-empirical method for prediction of antigenic determinants on protein antigens. <i>FEBS Lett</i>. 1990;276(1-2):172-176.</li>
<li>Parker JM, Guo D, Hodges RS. New hydrophilicity scale derived from high-performance liquid chromatography peptide retention data. <i>Biochemistry</i>. 1986;25(19):5425-5432.</li>
<li>Emini EA, Hughes JV, Perdue DS, et al. Induction of hepatitis A virus-neutralizing antibody by a virus-specific synthetic peptide. <i>J Virol</i>. 1985;55(3):836-839.</li>
<li>Chou PY, Fasman GD. Prediction of the secondary structure of proteins from their amino acid sequence. <i>Adv Enzymol</i>. 1978;47:45-148.</li>
<li>Chothia C, Lesk AM. Canonical structures for the hypervariable regions of immunoglobulins. <i>J Mol Biol</i>. 1987;196(4):901-917.</li>
<li>Lefranc MP, Giudicelli V, Ginestoux C, et al. IMGT, the international ImMunoGeneTics information system. <i>Nucleic Acids Res</i>. 2009;37(suppl_1):D1006-D1012.</li>
<li>Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. <i>Nucleic Acids Res</i>. 2000;28(1):235-242.</li>
<li>Cock PJ, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. <i>Bioinformatics</i>. 2009;25(11):1422-1423.</li>
<li>Nielsen M, Lundegaard C, Blicher T, et al. NetMHCpan, a method for quantitative prediction of peptide binding to any HLA-A and -B locus. <i>PLoS One</i>. 2008;3(3):e1797.</li>
<li>Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database (IEDB) 3.0: the central repository of epitope data. <i>Nucleic Acids Res</i>. 2019;47(D1):D448-D454.</li>
<li>World Health Organization. International Nonproprietary Names (INN) for biological and biotechnological substances. <i>INN Working Document 05.179</i>. 2014.</li>
<li>Immunogenicity Database Collaborative. IDC DS V1: Clinical immunogenicity dataset of marketed biologics. <i>CC BY 4.0</i>. https://immunogenicitydb.org/</li>
<li>Somanathan A, Bhatt L, Ferreira R, et al. Therapeutic monoclonal antibody nomenclature: the INN system and immunogenicity implications. <i>mAbs</i>. 2024;16(1):2294240.</li>
<li>Jawa V, Cousens LP, Chan A, et al. A review of factors affecting the immunogenicity of therapeutic proteins. <i>J Immunotoxicol</i>. 2013;10(2):125-137.</li>
<li>Rosenberg AS, Saunas T. Immunogenicity of biological therapeutics: a cause for concern? <i>Expert Rev Clin Immunol</i>. 2022;18(1):1-3.</li>
</ol>`;
}

function methodsEn() {
  return `
<h3>1. Engine & Data</h3>
<p>The prediction engine is <b>immunopred</b>: it extracts multi-dimensional features from the input protein sequence (MHC-I/II epitopes, B-cell linear epitopes, antibody CDR/framework germline alignment, physico-chemical aggregation propensity, and human-homology foreignness), then calibrates them into a 0–100% <b>sequence-intrinsic immunogenicity risk</b>. The reference library is <code>ada_reference.csv</code> (289 marketed biologics with ADA/NAb rates, sequences, heavy/light chains); it is used <b>only</b> for clinical background comparison and the literature factor of the neutralization level — it <b>does not participate</b> in the sequence risk score.</p>

<h3>2. Score Composition</h3>
<p><b>Step 1 · Raw epitope signal weighting:</b> MHC-II epitope density 45% + MHC-II strong-binding score 20% + B-cell epitope coverage 20% + MHC-I epitope density 8% + aggregation propensity 7% — yielding the pre-calibration full-length epitope load.</p>
<p><b>Step 2 · Path calibration:</b> For antibodies, CDR/germline calibration is applied — sequence foreignness (scaffold + FR/CDR germline mismatch + format exposure) × CDR-weighted epitope burden × lineage tolerance factor (fully human 0.55 / humanized 0.40 / chimeric 1.15) gives the baseline; multi-signal cross-confirmation then gives a capped uplift (+15% max). Non-antibody proteins use class-specific foreignness × epitope × class tolerance.</p>
<p>Classification: <b>&lt;20% Low / 20–45% Moderate / ≥45% High</b>.</p>

<h3>3. T-cell epitopes</h3>
<ul>
<li><b>MHC-I (CD8):</b> 9-mer sliding window, pocket-preference matrices for 7 common HLA-A/B alleles (anchor residues P2 / C-terminus), ranked by intra-protein percentile; top 10% with score ≥0.15 flagged as strong binders.</li>
<li><b>MHC-II (CD4):</b> 15-mer sliding window + 9-mer core scan, PSSM for P1/P4/P6/P9 pockets of 6 HLA-DRB1 alleles; strong-binding threshold = allele threshold_strong × 0.5 and percentile ≥90.</li>
</ul>

<h3>4. B-cell epitopes</h3>
<p>Linear-epitope consensus score = mean z-score of Kolaskar antigenicity + Parker hydrophilicity + Emini surface accessibility + Chou-Fasman β-turn; consecutive ≥6 aa with consensus ≥0.5 (z-score) are called candidate epitopes. Linear scoring has limited power to distinguish surface exposure (&gt;90% of B-cell epitopes are conformational); it is used only for candidate screening.</p>

<h3>5. Antibody calibration</h3>
<p>CDRs are delimited by conserved Cys/Trp and WGXG/FGXG anchor heuristics (Chothia-like); aligned against an IMGT human V-gene pool to obtain FR/CDR mismatch fractions and best germline hit. Lineage is inferred from both INN suffix (-ximab/-zumab/-umab) and germline FR. Multi-signal cross-validation (elevated CDR non-human-germline epitopes, MHC-II peptides covering germline-heterologous residues, etc.) can trigger epitope uplift.</p>

<h3>6. Neutralizing antibody (NAb) risk</h3>
<p>Independent score: structure/class factor + neighbor-drug ADA incidence factor + literature/annotation factor. Levels: Low (&lt;25%) / Moderate (25–55%) / High (≥55%). The prediction is a tendency for neutralization likelihood, not a confirmed neutralizing titer.</p>

<h3>7. Sequence-reference match confidence</h3>
<p>The reference table is aligned in a chain-aware manner: Fv (heavy + light dual-chain) high similarity (≥95%) is a high-confidence same-drug hit; VH-only or VL-only single-chain similarity (human Ab frameworks commonly reach 85%+ mutual similarity) is explicitly labeled "single-chain framework only" and is not equivalent to the same drug; k-mer local-fragment similarity is only a hint. Clinical ADA values are never fed back into the score under any circumstances.</p>

<h3>8. Limitations & Disclaimer</h3>
<ul>
<li>MHC binding uses built-in heuristics (pocket-preference matrices) and cannot replace NetMHCpan / IEDB / wet-lab assays.</li>
<li>Clinical ADA is influenced by route, frequency, concomitant medication, assay, and population HLA; it cannot be directly converted from a sequence score.</li>
<li>CDR delineation is heuristic-anchor-based, not IMGT formal numbering; for formal filings verify with ANARCI.</li>
<li>This tool is for research and early candidate screening only — <b>it does not constitute a basis for clinical or regulatory decisions</b>.</li>
</ul>

<h3>References</h3>
<ol class="refs">
<li>Kolaskar AS, Tongaonkar PC. A semi-empirical method for prediction of antigenic determinants on protein antigens. <i>FEBS Lett</i>. 1990;276(1-2):172-176.</li>
<li>Parker JM, Guo D, Hodges RS. New hydrophilicity scale derived from high-performance liquid chromatography peptide retention data. <i>Biochemistry</i>. 1986;25(19):5425-5432.</li>
<li>Emini EA, Hughes JV, Perdue DS, et al. Induction of hepatitis A virus-neutralizing antibody by a virus-specific synthetic peptide. <i>J Virol</i>. 1985;55(3):836-839.</li>
<li>Chou PY, Fasman GD. Prediction of the secondary structure of proteins from their amino acid sequence. <i>Adv Enzymol</i>. 1978;47:45-148.</li>
<li>Chothia C, Lesk AM. Canonical structures for the hypervariable regions of immunoglobulins. <i>J Mol Biol</i>. 1987;196(4):901-917.</li>
<li>Lefranc MP, Giudicelli V, Ginestoux C, et al. IMGT, the international ImMunoGeneTics information system. <i>Nucleic Acids Res</i>. 2009;37(suppl_1):D1006-D1012.</li>
<li>Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. <i>Nucleic Acids Res</i>. 2000;28(1):235-242.</li>
<li>Cock PJ, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. <i>Bioinformatics</i>. 2009;25(11):1422-1423.</li>
<li>Nielsen M, Lundegaard C, Blicher T, et al. NetMHCpan, a method for quantitative prediction of peptide binding to any HLA-A and -B locus. <i>PLoS One</i>. 2008;3(3):e1797.</li>
<li>Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database (IEDB) 3.0: the central repository of epitope data. <i>Nucleic Acids Res</i>. 2019;47(D1):D448-D454.</li>
<li>World Health Organization. International Nonproprietary Names (INN) for biological and biotechnological substances. <i>INN Working Document 05.179</i>. 2014.</li>
<li>Immunogenicity Database Collaborative. IDC DS V1: Clinical immunogenicity dataset of marketed biologics. <i>CC BY 4.0</i>. https://immunogenicitydb.org/</li>
<li>Somanathan A, Bhatt L, Ferreira R, et al. Therapeutic monoclonal antibody nomenclature: the INN system and immunogenicity implications. <i>mAbs</i>. 2024;16(1):2294240.</li>
<li>Jawa V, Cousens LP, Chan A, et al. A review of factors affecting the immunogenicity of therapeutic proteins. <i>J Immunotoxicol</i>. 2013;10(2):125-137.</li>
<li>Rosenberg AS, Saunas T. Immunogenicity of biological therapeutics: a cause for concern? <i>Expert Rev Clin Immunol</i>. 2022;18(1):1-3.</li>
</ol>`;
}

/* -------------------------------- init ------------------------------ */
applyI18n();
initForm();
renderMethods();
loadMeta();
