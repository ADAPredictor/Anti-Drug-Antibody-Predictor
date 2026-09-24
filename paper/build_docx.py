# -*- coding: utf-8 -*-
"""Build a Bioinformatics (Oxford University Press) submission-styled DOCX
from ADA_manuscript_en.md. Times New Roman, A4 margins, single-paragraph
abstract, OUP numbered references with hanging indent, and the required
back-matter sections (Acknowledgements / Funding / Conflict of Interest)."""
import re, sys
from pathlib import Path
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

PAPER = Path(__file__).resolve().parent
LANG = sys.argv[1].lower() if len(sys.argv) > 1 else "en"
if LANG == "zh":
    SRC = str(PAPER / "ADA_manuscript_zh.md")
    OUT = str(PAPER / "ADA_manuscript_Bioinformatics_zh.docx")
    ABSTRACT_H = "摘要"
    REFS_H = "参考文献"
    KW_H = "关键词："
    CJK_FONT = "宋体"
    BACK = [
        ("致谢", "[投稿前填写——按 OUP 要求致谢贡献者、技术协助与资助机构。]"),
        ("资助", "[投稿前填写——列出资助来源、项目编号及各资助方作用；如无，注明\"未获资助\"。]"),
        ("利益冲突", "[投稿前填写——声明是否存在竞争性财务或非财务利益；如无，注明\"无\"。]"),
    ]
else:
    SRC = str(PAPER / "ADA_manuscript_en.md")
    OUT = str(PAPER / "ADA_manuscript_Bioinformatics.docx")
    ABSTRACT_H = "Abstract"
    REFS_H = "References"
    KW_H = "Keywords: "
    CJK_FONT = "Times New Roman"
    BACK = [
        ("Acknowledgements", "[To be completed before submission — acknowledge contributors, technical assistance, and funding bodies per OUP requirements.]"),
        ("Funding", "[To be completed — list funding sources, grant numbers, and the role of each funder. If none, state 'No funding to declare'.]"),
        ("Conflict of Interest", "[To be completed — declare any competing financial or non-financial interests. If none, state 'None declared'.]"),
    ]
BODY_PT = 11
LINE = 1.5
FONT = "Times New Roman"

# ---- inline markup -> runs ----
INLINE_RE = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)")

def add_runs(p, text, base_pt=BODY_PT, base_bold=False, base_italic=False):
    """Add runs to paragraph p, honoring **bold**, *italic*, `code`."""
    if text is None:
        text = ""
    pos = 0
    for m in INLINE_RE.finditer(text):
        if m.start() > pos:
            _run(p, text[pos:m.start()], base_pt, base_bold, base_italic)
        seg = m.group(0)
        if seg.startswith("**"):
            _run(p, seg[2:-2], base_pt, True, base_italic)
        elif seg.startswith("*"):
            _run(p, seg[1:-1], base_pt, base_bold, True)
        elif seg.startswith("`"):
            r = _run(p, seg[1:-1], base_pt, base_bold, base_italic)
            r.font.name = "Consolas"
        pos = m.end()
    if pos < len(text):
        _run(p, text[pos:], base_pt, base_bold, base_italic)

def _run(p, text, pt, bold, italic):
    r = p.add_run(text)
    r.font.size = Pt(pt)
    r.font.name = FONT
    r.bold = bold
    r.italic = italic
    r.font.element.rPr.rFonts.set(qn("w:eastAsia"), CJK_FONT)
    return r

def clean_formula(tex):
    """Faithful readability pass on LaTeX: keep meaning, strip display clutter."""
    s = tex
    s = s.replace(r"\text{", "").replace(r"\cdot", "·")
    # remove leftover braces that wrapped \text content
    s = s.replace("{", "").replace("}", "")
    s = s.replace(r"\%", "%")
    s = s.strip()
    return s

# ---- document setup ----
def setup_doc():
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = FONT
    style.font.size = Pt(BODY_PT)
    style.element.rPr.rFonts.set(qn("w:eastAsia"), CJK_FONT)
    pf = style.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    pf.line_spacing = LINE
    pf.space_after = Pt(6)
    sec = doc.sections[0]
    sec.page_height = Cm(29.7)   # A4
    sec.page_width = Cm(21.0)
    sec.top_margin = Cm(2.5)
    sec.bottom_margin = Cm(2.5)
    sec.left_margin = Cm(2.5)
    sec.right_margin = Cm(2.5)
    return doc

def set_hanging(p, cm=0.75):
    pPr = p._p.get_or_add_pPr()
    ind = pPr.find(qn("w:ind"))
    if ind is None:
        ind = OxmlElement("w:ind")
        pPr.append(ind)
    ind.set(qn("w:left"), str(int(cm * 567)))
    ind.set(qn("w:hanging"), str(int(cm * 567)))

def shade_cell(cell, fill="D9D9D9"):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:fill"), fill)
    tcPr.append(shd)

def add_table(doc, rows):
    if not rows:
        return
    ncols = len(rows[0])
    t = doc.add_table(rows=len(rows), cols=ncols)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.style = "Table Grid"
    for i, row in enumerate(rows):
        for j, cell_txt in enumerate(row):
            cell = t.rows[i].cells[j]
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.line_spacing = 1.0
            add_runs(p, cell_txt, base_pt=10, base_bold=(i == 0))
            if i == 0:
                shade_cell(cell)
    # spacing after table
    doc.add_paragraph().paragraph_format.space_after = Pt(2)

def parse_table(lines, i):
    """Parse a pipe table starting at line i; return (rows, next_i)."""
    rows = []
    while i < len(lines) and lines[i].strip().startswith("|"):
        line = lines[i].strip()
        if re.match(r"^\|[\s:|-]+\|$", line) and set(line) <= set("|: -"):
            i += 1
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        rows.append(cells)
        i += 1
    return rows, i

# ---- main build ----
def main():
    with open(SRC, encoding="utf-8") as f:
        text = f.read()
    lines = text.split("\n")

    doc = setup_doc()
    i = 0
    in_abstract = False
    abstract_paras = []      # collected abstract body paragraphs (minus keywords)
    keywords_line = None
    refs_started = False
    ref_items = []
    pending_blocks = []      # normal content blocks
    back_matter_in_md = False  # set True if markdown already contains back-matter sections

    # first pass: collect abstract specially and detect refs
    # We'll do a single linear pass with state.
    section_title = None
    title_done = False

    def flush_normal(block):
        pending_blocks.append(block)

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        s = line.strip()

        if s == "":
            i += 1
            continue

        # horizontal rule -> skip
        if s == "---":
            i += 1
            continue

        # title
        if s.startswith("# ") and not title_done:
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(12)
            r = p.add_run(s[2:].strip())
            r.font.size = Pt(15)
            r.bold = True
            r.font.name = FONT
            title_done = True
            i += 1
            continue

        # author meta lines (start with **Word:**) — but not inside the abstract
        if re.match(r"^\*\*[^*]+?[:：]\*\*", s) and not s.startswith("|") and not in_abstract:
            label_m = re.match(r"^\*\*([^*]+?)[:：]\*\*\s*(.*)$", s)
            if label_m:
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_after = Pt(2)
                _run(p, label_m.group(1) + ": ", BODY_PT, True, False)
                add_runs(p, label_m.group(2), base_pt=BODY_PT)
                i += 1
                continue

        # abstract section
        if s.startswith("## Abstract") or s.startswith("## 摘要"):
            p = doc.add_paragraph()
            r = p.add_run(ABSTRACT_H)
            r.bold = True
            r.font.size = Pt(BODY_PT)
            r.font.name = FONT
            p.paragraph_format.space_after = Pt(4)
            in_abstract = True
            i += 1
            continue

        if in_abstract:
            # end of abstract when next ## or end
            if s.startswith("## "):
                in_abstract = False
                # don't increment; let the ## handler below catch it
            else:
                # collect abstract content
                m = re.match(r"^\*\*([^*]+?)[:：]\*\*\s*(.*)$", s)
                if m:
                    label = m.group(1)
                    content = m.group(2)
                    if "keyword" in label.lower() or "关键词" in label:
                        keywords_line = content
                    else:
                        abstract_paras.append(content)
                else:
                    abstract_paras.append(s)
                i += 1
                continue

        # references section
        if s.startswith("## References") or s.startswith("## 参考文献"):
            # first, render back-matter placeholders before references (OUP order)
            # skip if markdown already contains back-matter sections
            if not back_matter_in_md:
                add_back_matter(doc)
            p = doc.add_paragraph()
            r = p.add_run(REFS_H)
            r.bold = True
            r.font.size = Pt(BODY_PT)
            r.font.name = FONT
            p.paragraph_format.space_before = Pt(10)
            refs_started = True
            i += 1
            continue

        if refs_started:
            # numbered reference "N. Author..."
            mref = re.match(r"^(\d+)\.\s+(.*)$", s)
            if mref:
                num = mref.group(1)
                body_ref = mref.group(2)
                p = doc.add_paragraph()
                p.paragraph_format.space_after = Pt(4)
                p.paragraph_format.line_spacing = 1.0
                set_hanging(p, 0.8)
                _run(p, "[" + num + "] ", 10, False, False)
                add_runs(p, body_ref, base_pt=10)
            else:
                # continuation line of a reference -> append to last? just add as small para
                p = doc.add_paragraph()
                p.paragraph_format.space_after = Pt(2)
                p.paragraph_format.line_spacing = 1.0
                add_runs(p, s, base_pt=10)
            i += 1
            continue

        # major heading ##
        if s.startswith("## "):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(6)
            p.paragraph_format.keep_with_next = True
            r = p.add_run(s[3:].strip())
            r.bold = True
            r.font.size = Pt(12)
            r.font.name = FONT
            # detect back-matter headings present in markdown
            if any(bm in s for bm in ["Acknowledgements", "Funding", "Conflict of Interest", "Data", "致谢", "资助", "利益冲突", "数据与代码"]):
                back_matter_in_md = True
            i += 1
            continue

        # sub heading ###
        if s.startswith("### "):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(8)
            p.paragraph_format.space_after = Pt(4)
            p.paragraph_format.keep_with_next = True
            r = p.add_run(s[4:].strip())
            r.bold = True
            r.italic = True
            r.font.size = Pt(BODY_PT)
            r.font.name = FONT
            i += 1
            continue

        # pipe table
        if s.startswith("|"):
            rows, i = parse_table(lines, i)
            add_table(doc, rows)
            continue

        # formula block $$...$$
        if s.startswith("$$"):
            # could be inline single or multi-line
            formula = s.strip("$").strip()
            j = i + 1
            if not s.endswith("$$") or s == "$$":
                while j < n and not lines[j].strip().endswith("$$"):
                    formula += " " + lines[j].strip()
                    j += 1
                if j < n:
                    formula += " " + lines[j].strip().rstrip("$").strip()
                    j += 1
            else:
                j = i + 1
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(4)
            p.paragraph_format.space_after = Pt(4)
            r = p.add_run(clean_formula(formula))
            r.italic = True
            r.font.size = Pt(BODY_PT)
            r.font.name = "Cambria Math"
            i = j
            continue

        # bullet list
        if s.startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            p.paragraph_format.space_after = Pt(3)
            add_runs(p, s[2:].strip())
            i += 1
            continue

        # numbered list "1. "
        mnl = re.match(r"^\d+\.\s+(.*)$", s)
        if mnl:
            p = doc.add_paragraph(style="List Number")
            p.paragraph_format.space_after = Pt(3)
            add_runs(p, mnl.group(1))
            i += 1
            continue

        # normal paragraph (may span until blank / next markup)
        para_buf = s
        j = i + 1
        while j < n:
            nxt = lines[j].strip()
            if nxt == "" or nxt.startswith("#") or nxt.startswith("|") or \
               nxt.startswith("- ") or nxt.startswith("$$") or \
               re.match(r"^\d+\.\s+", nxt) or nxt == "---":
                break
            para_buf += " " + nxt
            j += 1
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        add_runs(p, para_buf)
        i = j

    # render abstract (single unstructured paragraph) BEFORE the rest?
    # The abstract was collected; render it now inserted right after title block.
    # python-docx appended normal content after abstract heading already;
    # but we deferred abstract body. Insert it as one paragraph at the doc top area.
    # Simpler: we already added "Abstract" heading inline; now add the merged body
    # as a paragraph. However order: heading was added when '## Abstract' seen,
    # which was before sections. So we can just append the merged paragraph now
    # — but that would put it at the END. Need to insert right after the Abstract heading.
    # To keep correct order, we instead render abstract body immediately by reordering.
    # Workaround: build abstract paragraph and move it to follow the Abstract heading.

    # Build the merged abstract paragraph and move it right after the 'Abstract' heading.
    abs_para = doc.add_paragraph()
    abs_para.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    abs_para.paragraph_format.space_after = Pt(8)
    merged = " ".join(abstract_paras)
    add_runs(abs_para, merged)
    # move abs_para to just after the 'Abstract' heading paragraph
    # find index of the paragraph containing run text 'Abstract' that is a heading
    body_el = abs_para._p
    head = None
    for p in doc.paragraphs:
        if p.text.strip() == ABSTRACT_H:
            head = p._p
            break
    if head is not None:
        head.addnext(body_el)

    # keywords paragraph (after abstract body)
    if keywords_line:
        kp = doc.add_paragraph()
        kp.paragraph_format.space_after = Pt(8)
        _run(kp, KW_H, BODY_PT, True, False)
        add_runs(kp, keywords_line, base_pt=BODY_PT)
        body_el.addnext(kp._p)

    doc.save(OUT)
    print("WROTE", OUT)
    print("abstract words:", len(merged.split()))

def add_back_matter(doc):
    """Insert OUP-required back-matter sections with clearly-marked placeholders."""
    for heading, body in BACK:
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(10)
        p.paragraph_format.space_after = Pt(2)
        r = p.add_run(heading)
        r.bold = True
        r.font.size = Pt(BODY_PT)
        r.font.name = FONT
        bp = doc.add_paragraph()
        bp.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        rr = bp.add_run(body)
        rr.italic = True
        rr.font.size = Pt(BODY_PT)
        rr.font.name = FONT

if __name__ == "__main__":
    main()
