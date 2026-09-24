"""Fetch immunogenicity (ADA) evidence from FDA labels via openFDA + SPL XML.

Both endpoints are public and require no API key for light use:
* https://api.fda.gov/drug/label.json   - locate the current label
* https://www.accessdata.fda.gov/spl/data/{set_id}/{set_id}.xml - full SPL
"""
import re

from ..config import (OPENFDA_LABEL_URL, SPL_XML_URL, SPL_XML_FALLBACK_URL,
                      HTTP_TIMEOUT, USER_AGENT)

_PCT_RE = re.compile(r"(\d+(?:\.\d+)?)\s*%")
_FRAC_RE = re.compile(r"(\d{1,4})\s*(?:of|/)\s*(\d{1,4})")
_TAG_RE = re.compile(r"<[^>]+>")
_WS_RE = re.compile(r"[ \t]+")
_NEXT_HEADING_RE = re.compile(
    r"(?:^|\n)\s*(?:\d{1,2}\.\d{1,2}|\d{1,2})\s+[A-Z][A-Za-z][A-Za-z &/-]{4,}"
)


def _http_get(url, params=None):
    import requests
    resp = requests.get(url, params=params,
                        headers={"User-Agent": USER_AGENT,
                                 "Accept": "application/json"},
                        timeout=HTTP_TIMEOUT)
    return resp


def find_label(brand):
    """Return (set_id, brand_name) or (None, None)."""
    for query in ('openfda.brand_name.exact:"%s"' % brand.upper(),
                  'openfda.brand_name:"%s"' % brand.lower()):
        try:
            resp = _http_get(OPENFDA_LABEL_URL,
                             params={"search": query, "limit": 1})
            if resp.status_code != 200:
                continue
            results = resp.json().get("results", [])
            if results:
                names = results[0].get("openfda", {}).get("brand_name", [brand])
                return results[0].get("set_id"), names[0]
        except Exception:
            continue
    return None, None


def fetch_spl_text(set_id):
    headers = {"User-Agent": "Mozilla/5.0 (compatible; ADA-research/0.1)",
               "Accept": "application/xml,text/xml,*/*"}
    import requests
    last_exc = None
    for template in (SPL_XML_URL, SPL_XML_FALLBACK_URL):
        try:
            resp = requests.get(template.format(set_id=set_id),
                                headers=headers, timeout=HTTP_TIMEOUT)
            resp.raise_for_status()
            if "<title" in resp.text[:5000].lower() or "<section" in resp.text[:5000].lower():
                return resp.text
        except Exception as exc:
            last_exc = exc
    if last_exc:
        raise last_exc
    raise IOError("SPL XML 下载失败：%s" % set_id)


def extract_immunogenicity_section(xml_text):
    """Locate the Immunogenicity section in SPL XML and return plain text."""
    lower = xml_text.lower()
    title_match = None
    for m in re.finditer(r"<title[^>]*>(.*?)</title>", xml_text, re.S | re.I):
        if "immunogenicity" in m.group(1).lower():
            title_match = m
            break
    if title_match is None:
        return ""
    tail = xml_text[title_match.end():title_match.end() + 20000]
    stop = _NEXT_HEADING_RE.search(_TAG_RE.sub("\n", tail[:8000]))
    if stop:
        tail = tail[:max(0, _safe_xml_cut(tail, stop.start()))]
    text = _TAG_RE.sub(" ", tail)
    text = re.sub(r"\s*\n\s*", "\n", text)
    text = _WS_RE.sub(" ", text)
    return text.strip()


def _safe_xml_cut(tail, text_pos):
    """Map a position in tag-stripped text back to approximate raw position."""
    out = 0
    stripped = 0
    in_tag = False
    for ch in tail:
        if ch == "<":
            in_tag = True
            continue
        if ch == ">":
            in_tag = False
            continue
        if not in_tag:
            if stripped >= text_pos:
                return out
            stripped += 1
        out += 1
    return out


def _sentences(text):
    return [s.strip() for s in re.split(r"(?<=[.;])\s+", text) if s.strip()]


def extract_ada_values(section_text):
    """Parse percentages / fractions from antibody-related sentences.

    Returns dict(percentages=[...], neutralizing=[...], evidence=str)
    """
    antibody_sentences = [s for s in _sentences(section_text)
                          if re.search(r"antibod|immunogenic|ada\b", s, re.I)]
    percentages, neutralizing = [], []
    for sent in antibody_sentences:
        vals = [float(x) for x in _PCT_RE.findall(sent)]
        for m in _FRAC_RE.finditer(sent):
            a, b = int(m.group(1)), int(m.group(2))
            if 0 < a <= b and 5 <= b <= 100000:
                vals.append(round(100.0 * a / b, 2))
        vals = [v for v in vals if 0.0 <= v <= 100.0]
        if not vals:
            continue
        if re.search(r"neutraliz", sent, re.I):
            neutralizing.extend(vals)
        else:
            percentages.extend(vals)
    evidence = " ".join(antibody_sentences[:4])[:600]
    return {
        "percentages": sorted(set(percentages)),
        "neutralizing": sorted(set(neutralizing)),
        "evidence": evidence,
    }


def enrich_from_label(brand):
    """Full pipeline for one brand. Returns patch dict or None."""
    set_id, canonical_brand = find_label(brand)
    if not set_id:
        return None
    xml_text = fetch_spl_text(set_id)
    section = extract_immunogenicity_section(xml_text)
    extracted = extract_ada_values(section)
    return {
        "label_url": SPL_XML_URL.format(set_id=set_id),
        "label_evidence": extracted["evidence"],
        "ada_values": extracted["percentages"],
        "neutralizing_values": extracted["neutralizing"],
    }
