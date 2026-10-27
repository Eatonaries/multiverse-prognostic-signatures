# -*- coding: utf-8 -*-
"""Build the JBI submission package.

Produces a self-contained Word manuscript from manuscript_v0.1.md plus the
figure legends, and assembles the upload folder C:\\emupload_JBI.

JBI required parts (from the live Guide for Authors):
  Cover Letter, Graphical Abstract, Statement of Significance (in the
  Introduction), Declaration of Competing Interest, Title page,
  Structured abstract + Keywords, Body (Introduction, Related Work, Methods,
  Results, Discussion, Conclusion; <=6000 words; figures+tables <=8)
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

import docx
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(r"D:\deepseek\jieru\_bio4")
MS_MD = ROOT / "manuscript" / "manuscript_v0.1.md"
LEGENDS = ROOT / "manuscript" / "FigureLegends_JBI.md"
COVER = ROOT / "manuscript" / "CoverLetter_JBI.md"
FIG = ROOT / "results" / "figures"
PKG = Path(r"C:\emupload_JBI")

# ---------------------------------------------------------------- markdown ---
INLINE = re.compile(r"(\*\*.+?\*\*|\*.+?\*|`[^`]+`)")


def add_runs(par, text, base_size=11, base_bold=False):
    """Render markdown inline emphasis into runs."""
    for tok in INLINE.split(text):
        if not tok:
            continue
        bold, italic, mono = base_bold, False, False
        t = tok
        if t.startswith("**") and t.endswith("**") and len(t) > 4:
            bold, t = True, t[2:-2]
        elif t.startswith("*") and t.endswith("*") and len(t) > 2:
            italic, t = True, t[1:-1]
        elif t.startswith("`") and t.endswith("`") and len(t) > 2:
            mono, t = True, t[1:-1]
        # strip any remaining markdown emphasis markers
        t = t.replace("**", "").replace("*", "")
        r = par.add_run(t)
        r.bold = bold
        r.italic = italic
        r.font.size = Pt(base_size)
        if mono:
            r.font.name = "Consolas"
    return par


def md_table_to_doc(doc, lines, size=9):
    """Render a markdown pipe table."""
    rows = []
    for ln in lines:
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if all(re.fullmatch(r":?-{2,}:?", c) for c in cells if c):
            continue
        rows.append(cells)
    if not rows:
        return
    ncol = max(len(r) for r in rows)
    t = doc.add_table(rows=0, cols=ncol)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, r in enumerate(rows):
        cells = t.add_row().cells
        for j in range(ncol):
            txt = r[j] if j < len(r) else ""
            p = cells[j].paragraphs[0]
            p.paragraph_format.space_after = Pt(1)
            add_runs(p, txt, base_size=size, base_bold=(i == 0))
    doc.add_paragraph()


def render_md(doc, md: str, skip_h1=True, size=11):
    lines = md.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        s = ln.strip()
        if not s:
            i += 1
            continue
        if s.startswith("|") and i + 1 < len(lines) and lines[i + 1].strip().startswith("|"):
            block = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i])
                i += 1
            md_table_to_doc(doc, block, size=size - 2)
            continue
        if s.startswith("### "):
            p = doc.add_heading(s[4:].strip(), level=3)
            i += 1
            continue
        if s.startswith("## "):
            p = doc.add_heading(s[3:].strip(), level=2)
            i += 1
            continue
        if s.startswith("# ") and not skip_h1:
            doc.add_heading(s[2:].strip(), level=1)
            i += 1
            continue
        if s.startswith("---"):
            i += 1
            continue
        if s.startswith("> "):
            # blockquote -> indented italic paragraph
            buf = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                buf.append(lines[i].strip().lstrip(">").strip())
                i += 1
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.3)
            add_runs(p, " ".join(b for b in buf if b), base_size=size - 1)
            for r in p.runs:
                r.italic = True
                r.font.color.rgb = RGBColor(0x44, 0x44, 0x44)
            continue
        if re.match(r"^[-*]\s+", s) or re.match(r"^\d+\.\s+", s):
            m = re.match(r"^([-*]|\d+\.)\s+", s)
            ordered = m.group(1) != "-"
            p = doc.add_paragraph(style="List Number" if ordered else "List Bullet")
            add_runs(p, s[m.end():], base_size=size)
            i += 1
            continue
        # normal paragraph: join wrapped lines
        buf = [s]
        i += 1
        while i < len(lines):
            nxt = lines[i].strip()
            if (not nxt or nxt.startswith(("#", "|", ">", "---"))
                    or re.match(r"^[-*]\s+", nxt) or re.match(r"^\d+\.\s+", nxt)):
                break
            buf.append(nxt)
            i += 1
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(6)
        add_runs(p, " ".join(buf), base_size=size)


def strip_front_matter(md: str) -> str:
    """Drop the internal header (title-of-draft, target-journal note)."""
    lines = md.splitlines()
    out = []
    for ln in lines:
        if ln.startswith("# Manuscript draft"):
            continue
        if ln.startswith("**Target journal**"):
            continue
        if ln.startswith("Verified hard limits"):
            continue
        out.append(ln)
    return "\n".join(out)


def main():
    PKG.mkdir(parents=True, exist_ok=True)
    for old in PKG.glob("*"):
        if old.is_file():
            old.unlink()

    # ---------------------------------------------------------- main docx ---
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(11)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(1)
        s.left_margin = s.right_margin = Inches(1)

    md = strip_front_matter(MS_MD.read_text(encoding="utf-8"))
    render_md(doc, md, size=11)

    # page break before the references
    for p in doc.paragraphs:
        if p.text.strip() == "References" and p.style.name.startswith("Heading"):
            p.paragraph_format.page_break_before = True
            break

    out_main = PKG / "1_Manuscript_JBI.docx"
    doc.save(out_main)
    print("wrote", out_main, out_main.stat().st_size, "bytes")

    # --------------------------------------------------------- cover letter ---
    cd = Document()
    cd.styles["Normal"].font.name = "Times New Roman"
    cd.styles["Normal"].font.size = Pt(11)
    cover_md = COVER.read_text(encoding="utf-8")
    cover_md = re.sub(r"^# .*\n+", "", cover_md)
    cover_md = re.sub(r"^---$", "", cover_md, flags=re.M)
    render_md(cd, cover_md, size=11)
    out_cover = PKG / "2_Cover_Letter.docx"
    cd.save(out_cover)
    print("wrote", out_cover, out_cover.stat().st_size, "bytes")

    # ------------------------------------------------------------- figures ---
    figs = {
        "3_Figure_1_specification_curve.pdf": "fig1_specification_curve.pdf",
        "4_Figure_2_variance_decomposition.pdf": "fig2_variance_decomposition.pdf",
        "5_Figure_3_external_attenuation.pdf": "fig3_external_attenuation.pdf",
        "6_Figure_4_cross_cancer_swing.pdf": "fig4_cross_cancer_swing.pdf",
        "7_Figure_5_variance_two_panels.pdf": "fig5_variance_two_panels.pdf",
        "8_Figure_6_random_control.pdf": "fig6_random_control.pdf",
    }
    for dst, src in figs.items():
        s = FIG / src
        if s.exists():
            shutil.copy2(s, PKG / dst)
            print("copied", src, "->", dst)
        else:
            print("MISSING FIGURE", src)

    # ---------------------------------------------------- graphical abstract ---
    ga = FIG / "graphical_abstract.pdf"
    if ga.exists():
        shutil.copy2(ga, PKG / "9_Graphical_Abstract.pdf")
        print("copied graphical_abstract.pdf")
    else:
        print("MISSING graphical_abstract.pdf")

    # ---------------------------------------------------------- supplements ---
    supp = sorted(FIG.glob("supplementary_figS*.pdf"))
    for i, s in enumerate(supp, 1):
        shutil.copy2(s, PKG / f"10_Supplementary_Figure_S{i}.pdf")
        print("copied", s.name)

    # ------------------------------------------------------------ legends ----
    leg = LEGENDS.read_text(encoding="utf-8")
    (PKG / "Figure_Legends.md").write_text(leg, encoding="utf-8")
    print("wrote Figure_Legends.md")

    print("\npackage contents:")
    for f in sorted(PKG.iterdir()):
        print(f"   {f.name:44} {f.stat().st_size:>9,}")


if __name__ == "__main__":
    main()
