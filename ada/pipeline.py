"""End-to-end analysis pipeline: sequence -> features, epitopes, score."""
import difflib
import hashlib
import json

from . import db
from .config import (DEFAULT_CLASS1_ALLELES, DEFAULT_CLASS2_ALLELES,
                     MAX_SCAN_LENGTH)
from .epitopes import bcell, tcell
from .features import compute_features
from .scoring import clinical_context, compute_score, tier, TIER_LABELS
from .seqio import parse_fasta, total_length


def _segments_to_coverage(segments, length, offset=0):
    flags = [0] * length
    for seg in segments:
        for i in range(seg["start"], seg["end"]):
            if 0 <= offset + i < length:
                flags[offset + i] = 1
    return flags


def _shift_segments(segments, offset, chain_id):
    out = []
    for seg in segments:
        s = dict(seg)
        s["start"] += offset
        s["end"] += offset
        s["peptide_start"] = s.get("peptide_start", seg["start"]) + offset
        s["chain"] = chain_id
        out.append(s)
    return out


# 新颖性局部比对只统计至少 8 aa 的连续匹配段（与最小 T 细胞表位长度一致），
# 避免短查询在长蛋白上由零散单残基匹配造成的假高一致性。
_NOVELTY_MIN_BLOCK = 8


def _chain_similarity(seq_a, seq_b):
    """两条链的序列一致性估计，取全局比对与局部包含度的较大值。

    - 全局：difflib 的 2M/(na+nb)，长度相近时敏感；
    - 局部：长度 >=8 的连续匹配块占较短链的比例，使 VH/VL 片段或多肽
      能与完整抗体重链/融合蛋白比较（旧逻辑的长度门控会直接跳过这些对）。
    """
    matcher = difflib.SequenceMatcher(None, seq_a, seq_b)
    global_r = matcher.ratio()
    matched = sum(block.size for block in matcher.get_matching_blocks()
                  if block.size >= _NOVELTY_MIN_BLOCK)
    contained = matched / float(min(len(seq_a), len(seq_b)))
    return max(global_r, contained)


def _novelty(chains, reference_drugs):
    """(novelty, identity, best_drug) against reference drug sequences."""
    best = 0.0
    best_drug = None
    has_ref = False
    for drug in reference_drugs:
        if not drug.get("sequence"):
            continue
        try:
            ref_chains = parse_fasta(drug["sequence"])
        except Exception:
            continue
        has_ref = True
        for chain in chains:
            for ref in ref_chains:
                sim = _chain_similarity(chain.seq, ref.seq)
                if sim > best:
                    best = sim
                    best_drug = drug
    if not has_ref:
        return None
    return round(1.0 - best, 3), round(best, 3), best_drug


def analyze(sequence_text, modality="unknown",
            class1_alleles=None, run_class2=True, class2_alleles=None,
            uniprot_acc=None, pdb_id=None, use_cache=True, clinical=None):
    """Run the full prediction. Returns a JSON-serializable result dict."""
    chains = parse_fasta(sequence_text)
    length = total_length(chains)
    if length > MAX_SCAN_LENGTH:
        raise ValueError("序列总长 %d 超过在线分析上限 %d；请使用 CLI 批量模式。"
                         % (length, MAX_SCAN_LENGTH))

    class1_alleles = class1_alleles or DEFAULT_CLASS1_ALLELES
    class2_alleles = class2_alleles or DEFAULT_CLASS2_ALLELES
    cache_key = hashlib.sha256(json.dumps({
        "seq": sequence_text, "m": modality, "a1": class1_alleles,
        "a2": run_class2, "u": uniprot_acc, "p": pdb_id,
        "c": clinical,
    }, sort_keys=True).encode()).hexdigest()
    if use_cache:
        cached = db.cache_get(cache_key)
        if cached is not None:
            cached["cached"] = True
            return cached

    warnings = []
    if any("X" in c.seq for c in chains):
        warnings.append("序列含未知残基 X，相关区域的预测结果可信度降低。")

    features, patch_flags = compute_features(chains)

    t1_track = [0.0] * length
    t2_track = [0.0] * length
    b_track = [0.0] * length
    t1_segments, t2_segments, b_segments = [], [], []
    hits1 = windows1 = hits2 = windows2 = 0
    offset = 0
    t1_backend = None
    t2_backend = None
    best1_affinity = None
    best2_affinity = None
    scanned1 = scanned2 = 0
    for chain in chains:
        r1 = tcell.scan_class1(chain.seq, class1_alleles)
        t1_backend = r1["backend"]
        t1_track[offset:offset + len(chain.seq)] = r1["track"]
        t1_segments += _shift_segments(r1["segments"], offset, chain.chain_id)
        hits1 += r1["strong_hit_count"]
        windows1 += r1["windows"]
        scanned1 += r1["peptides_scanned"]
        b1 = r1["summary"].get("best_affinity_nm")
        if b1 is not None:
            best1_affinity = b1 if best1_affinity is None else min(best1_affinity, b1)

        if run_class2:
            r2 = tcell.scan_class2(chain.seq, class2_alleles)
            t2_backend = r2["backend"]
            t2_track[offset:offset + len(chain.seq)] = r2["track"]
            t2_segments += _shift_segments(r2["segments"], offset, chain.chain_id)
            hits2 += r2["strong_hit_count"]
            windows2 += r2["windows"]
            scanned2 += r2["peptides_scanned"]
            b2 = r2["summary"].get("best_affinity_nm")
            if b2 is not None:
                best2_affinity = b2 if best2_affinity is None else min(best2_affinity, b2)

        rb = bcell.predict_linear(chain.seq)
        b_track[offset:offset + len(chain.seq)] = rb["scores"]
        b_segments += _shift_segments(rb["segments"], offset, chain.chain_id)
        offset += len(chain.seq)

    mhc1_rate = hits1 / windows1 if windows1 else 0.0
    mhc2_rate = hits2 / windows2 if windows2 else 0.0
    bcell_cov = sum(_segments_to_coverage(b_segments, length))

    # optional experimental structure
    exposure = None
    if uniprot_acc or pdb_id:
        try:
            from . import structure
            pdb_text, source = structure.fetch_pdb_text(uniprot_acc, pdb_id)
            exposure = structure.compute_exposure(pdb_text, chains)
            exposure["source"] = source
            if exposure["mapped_fraction"] < 0.5:
                warnings.append("结构残基与查询序列匹配率仅 %.0f%%，表面暴露结果仅供参考。"
                                % (100 * exposure["mapped_fraction"]))
        except Exception as exc:
            warnings.append("结构获取/解析失败，已回退到纯序列模式：%s" % exc)

    novelty = _novelty(chains, db.all_drugs())
    clinical_ctx = clinical_context(clinical)
    novelty_value = novelty[0] if novelty else None
    exposed_fraction = exposure["exposed_fraction"] if exposure else None
    total, breakdown = compute_score(
        features, mhc1_rate, mhc2_rate, bcell_cov, exposed_fraction,
        modality, novelty_value,
    )

    # concatenate chains for the sequence track display
    display_seq = "".join(c.seq for c in chains)
    chain_ranges = []
    off = 0
    for c in chains:
        chain_ranges.append({"chain_id": c.chain_id, "name": c.name,
                             "start": off, "end": off + len(c.seq)})
        off += len(c.seq)

    result = {
        "input": {
            "modality": modality,
            "class1_alleles": class1_alleles,
            "class2_alleles": class2_alleles if run_class2 else [],
            "sequence": display_seq,
            "chains": chain_ranges,
            "structure_requested": bool(uniprot_acc or pdb_id),
        },
        "features": features,
        "epitopes": {
            "mhc1": {
                "backend": t1_backend,
                "track": t1_track,
                "segments": _serialize_segments(t1_segments),
                "n_segments": len(t1_segments),
                "peptides_scanned": scanned1,
                "best_affinity_nm": best1_affinity,
            },
            "mhc2": {
                "backend": t2_backend,
                "track": t2_track,
                "segments": _serialize_segments(t2_segments),
                "n_segments": len(t2_segments),
                "peptides_scanned": scanned2,
                "best_affinity_nm": best2_affinity,
            } if run_class2 else None,
            "bcell": {
                "track": b_track,
                "segments": _serialize_segments(b_segments),
                "n_segments": len(b_segments),
                "method": "parker+cf-turn-v2",
            },
            "exposure": exposure,
        },
        "score": {
            "total": total,
            "tier": tier(total),
            "tier_label": TIER_LABELS[tier(total)],
            "components": breakdown,
            "coverage": {"mhc1_strong_rate": round(mhc1_rate, 4),
                         "mhc2_strong_rate": round(mhc2_rate, 4),
                         "mhc1_strong_peptides": hits1,
                         "mhc2_strong_peptides": hits2,
                         "bcell_residues": bcell_cov},
            "nearest_reference_identity": novelty[1] if novelty else None,
            "nearest_reference": {
                "generic": novelty[2]["generic"],
                "identity": novelty[1],
                "modality": novelty[2].get("modality"),
                "ada_min_pct": novelty[2].get("ada_min_pct"),
                "ada_max_pct": novelty[2].get("ada_max_pct"),
                "ada_category": novelty[2].get("ada_category"),
                "label_evidence": novelty[2].get("label_evidence"),
                "source": novelty[2].get("source"),
            } if novelty and novelty[2] else None,
            "clinical_context": clinical_ctx,
            "adjusted_total": (round(min(100.0, max(0.0,
                                total + clinical_ctx["adjustment"])), 1)
                               if clinical_ctx else None),
        },
        "reference": db.modality_stats(),
        "warnings": warnings,
        "cached": False,
    }
    db.cache_set(cache_key, result)
    return result


def _serialize_segments(segments):
    out = []
    for s in segments:
        out.append({
            "start": s["start"],
            "end": s["end"],
            "chain": s.get("chain"),
            "peptide": s.get("peptide"),
            "peptide_start": s.get("peptide_start"),
            "affinity_nm": round(s["affinity"], 1) if s.get("affinity") else None,
            "level": s.get("level"),
            "strong_alleles": s.get("strong_alleles", []),
            "promiscuity": s.get("promiscuity"),
            "mean_parker": (round(s["mean_parker"], 2)
                            if s.get("mean_parker") is not None else None),
            "mean_turn": (round(s["mean_turn"], 2)
                          if s.get("mean_turn") is not None else None),
            "seq": s.get("seq"),
        })
    return out
