"""
Antibody format / architecture from sequence (not INN labels).

Detects IgG / Fab / scFv / VHH-like architectures and estimates
*format exposure* — surfaces that are buried in full IgG but solvent-
accessible in fragments (scFv neoepitopes → pre-existing ADA risk).

Important: many reference rows store VH+VL only (no Fc). Those look like
scFv by length alone — we treat dual-domain, linker-free, Fv-length
sequences as truncated IgG V-regions (low exposure) unless drug_class /
an engineered linker indicates a true fragment.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from .cdr import AntibodyRegionMap, annotate_antibody_regions

# Engineered (G4S)n linkers — definitive scFv evidence
_ENGINEERED_LINKER = re.compile(r"(?:GGGGS){1,4}|(?:GGGS){2,}|GGSGG|GSSGG")
# Short VH–VL junctions only trusted when product class says scFv / single-chain
_SHORT_JUNCTION = re.compile(
    r"(?:SS|AA|NA|GS){1,2}(?=[DEQARIKLMNPVSTWYG]{0,2}(?:EIV|DIQ|DIV|AIQ|SSE|QSA|QSV|DIQM))"
)

_FC_MOTIFS = (
    "EPKSCDK",
    "DKTHTCP",
    "PPCPAP",
    "ELLGGP",
    "GQPREPQ",
    "YNSTYR",
    "KCKVSNK",
)
_CH1_MOTIFS = ("ASTKGPS", "GPSVFPL", "APSVFPL")
_CL_MOTIFS = ("RTVAAP", "GQPKAAP", "KRTVAAP")

# Exposure 0–1: contribution to structural foreignness (format-driven)
_EXPOSURE = {
    "igg": 0.08,
    "fab": 0.28,
    "scfv": 0.72,
    "vhh": 0.55,
    "fv": 0.12,  # V-only reference for a full mAb — not a marketed scFv
    "unknown": 0.15,
}


@dataclass
class AbFormatReport:
    format: str  # igg | fab | scfv | vhh | fv | unknown
    exposure: float
    length: int
    n_domains: int
    has_linker: bool
    has_fc: bool
    has_constant: bool
    method_note_zh: str
    method_note_en: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "format": self.format,
            "exposure": round(float(self.exposure), 4),
            "length": int(self.length),
            "n_domains": int(self.n_domains),
            "has_linker": bool(self.has_linker),
            "has_fc": bool(self.has_fc),
            "has_constant": bool(self.has_constant),
            "method_note_zh": self.method_note_zh,
            "method_note_en": self.method_note_en,
        }


def _has_any(seq: str, motifs: tuple[str, ...]) -> bool:
    return any(m in seq for m in motifs)


def detect_ab_format(
    sequence: str,
    regions: AntibodyRegionMap | None = None,
    *,
    inn: str = "",
    drug_class: str = "",
) -> AbFormatReport:
    """Infer Ab architecture and format-exposure score from sequence (+ optional class)."""
    seq = "".join(c for c in (sequence or "").upper() if c.isalpha())
    L = len(seq)
    blob = f"{inn} {drug_class}".lower()
    reg = regions if regions is not None else annotate_antibody_regions(seq)
    n_dom = len(reg.domains) if reg and reg.domains else 0
    has_fc = _has_any(seq, _FC_MOTIFS)
    has_ch1 = _has_any(seq, _CH1_MOTIFS)
    has_cl = _has_any(seq, _CL_MOTIFS)
    has_const = has_ch1 or has_cl or has_fc
    has_engineered_linker = bool(_ENGINEERED_LINKER.search(seq))
    class_scfv = bool(
        re.search(r"\bscfv\b|single[\s-]?chain|single[\s-]?domain antibody", blob)
    )
    class_fab = bool(re.search(r"\bfab\b|antigen[\s-]?binding fragment", blob))
    class_vhh = bool(re.search(r"\bvhh\b|nanobody", blob))
    class_mab = bool(
        re.search(
            r"monoclonal|whole mab|\bigg\b|chimeric|humanized|fully human|"
            r"(?<!single-domain )(?<!single domain )\bantibody\b|\bmab\b",
            blob,
        )
    ) and not class_scfv and not class_fab and not class_vhh
    # Short junction counts as linker only for known scFv-class products (e.g. Bro)
    has_linker = has_engineered_linker or (
        class_scfv and bool(_SHORT_JUNCTION.search(seq))
    )

    fmt = "unknown"
    if has_fc or L >= 600:
        fmt = "igg"
    elif class_fab or (has_const and L >= 400 and not has_fc):
        fmt = "fab"
    elif class_vhh or (n_dom == 1 and 100 <= L <= 160 and not class_mab):
        fmt = "vhh"
    elif class_scfv or (has_engineered_linker and n_dom >= 2 and L <= 360):
        fmt = "scfv"
    elif n_dom >= 2 and L <= 360 and not has_const:
        fmt = "fv"
    elif class_mab:
        fmt = "igg"
    elif reg and reg.is_antibody_like and L < 400:
        fmt = "fv"
    elif L >= 400:
        fmt = "igg"

    # Full mAb product with V-only reference sequence → IgG exposure, not scFv
    if class_mab and fmt in {"scfv", "fv", "vhh", "unknown"} and not has_engineered_linker:
        fmt = "igg"
    if class_scfv and fmt in {"fv", "unknown", "vhh", "igg"} and L <= 400 and not has_fc:
        fmt = "scfv"

    exposure = float(_EXPOSURE.get(fmt, _EXPOSURE["unknown"]))
    if fmt == "scfv" and has_linker:
        exposure = min(1.0, exposure + 0.05)

    note_zh = (
        f"架构检测：格式={fmt}，长度={L}，可变域数={n_dom}，"
        f"linker={'有' if has_linker else '无'}，Fc={'有' if has_fc else '无'}，"
        f"恒定区指纹={'有' if has_const else '无'}；"
        f"格式暴露面={exposure:.2f}（真 scFv 抬升；库内仅 VH+VL 的全抗按低暴露处理）。"
    )
    note_en = (
        f"Architecture: format={fmt}, len={L}, V-domains={n_dom}, "
        f"linker={'yes' if has_linker else 'no'}, Fc={'yes' if has_fc else 'no'}, "
        f"constant motifs={'yes' if has_const else 'no'}; "
        f"format exposure={exposure:.2f} (true scFv elevated; Fv-only full-mAb refs stay low)."
    )
    return AbFormatReport(
        format=fmt,
        exposure=exposure,
        length=L,
        n_domains=n_dom,
        has_linker=has_linker,
        has_fc=has_fc,
        has_constant=has_const,
        method_note_zh=note_zh,
        method_note_en=note_en,
    )


def format_exposure_score(
    sequence: str,
    regions: AntibodyRegionMap | None = None,
    *,
    inn: str = "",
    drug_class: str = "",
) -> float:
    return float(
        detect_ab_format(sequence, regions, inn=inn, drug_class=drug_class).exposure
    )
