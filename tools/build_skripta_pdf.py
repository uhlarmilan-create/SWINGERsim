"""Build the Czech textbook (skripta) docs/skripta/cs/*.md into docs/pdf/SWINGERsim-skripta-cs.pdf with pandoc and
XeLaTeX (fonts: DejaVu, Latin Modern Math, FreeSerif).

    python tools/build_skripta_pdf.py

The chapters are written for the Claude Doc of the skripta; here their headings move one level up (## -> chapter),
```latex blocks become display math, wide displays are scaled to the line and long inline code may break.
"""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = str(ROOT / "docs/skripta/cs")
OUT = str(ROOT / "build/skripta")
PDF = ROOT / "docs/pdf/SWINGERsim-skripta-cs.pdf"
os.makedirs(OUT, exist_ok=True)
ORDER = ["ch01", "ch02", "ch03", "ch04", "ch05", "ch06", "ch07", "ch08", "ch09", "ch10", "ch11", "ch12", "prilohy"]


def fix(md: str) -> str:
    out, in_code, lang = [], False, None
    for line in md.split("\n"):
        m = re.match(r"^(\s*)```(\w*)\s*$", line)
        if m and not in_code:
            in_code, lang = True, m.group(2)
            out.append("$$" if lang == "latex" else line)
            continue
        if m and in_code and m.group(2) == "":
            in_code = False
            out.append("$$" if lang == "latex" else line)
            continue
        if in_code:
            out.append(line)
            continue
        # headings one level up (## Kapitola -> chapter)
        h = re.match(r"^(#{2,6}) ", line)
        if h:
            line = line[1:]
        out.append(line)
    md = "\n".join(out)

    md = re.sub(r"!\[([^\]]*)\]\((figures/[^)]+)\)[ \t]*",
                lambda m: f"![]({SRC}/{m.group(2)}){{width=96%}}\n\n", md)
    return md


parts = []
for name in ORDER:
    parts.append(fix(open(f"{SRC}/{name}.md", encoding="utf-8").read()))
body = "\n\n".join(parts)
# Úvod and Přílohy unnumbered in the TOC sense: chapters carry their own numbers in the title
meta = """---
title: "SWINGERsim"
subtitle: "Fyzika a numerické modelování vícenásobných hvězd"
date: "listopad 2026"
lang: cs
---

"""
open(f"{OUT}/skripta.md", "w", encoding="utf-8").write(meta + body)

header = r"""
\usepackage{newunicodechar}
\newfontfamily\fbmath{DejaVu Math TeX Gyre}
\newfontfamily\fbsym{FreeSerif}[ItalicFont=FreeSerif,BoldItalicFont=FreeSerif Bold]
\newfontfamily\fbsans{DejaVu Sans}
\newunicodechar{ᐟ}{\ifmmode\text{{\fbsans\symbol{"141F}}}\else{\fbsans\symbol{"141F}}\fi}
\newunicodechar{𝒜}{\ifmmode\text{{\fbmath\symbol{"1D49C}}}\else{\fbmath\symbol{"1D49C}}\fi}
\newunicodechar{𝒩}{\ifmmode\text{{\fbmath\symbol{"1D4A9}}}\else{\fbmath\symbol{"1D4A9}}\fi}
\newunicodechar{𝒯}{\ifmmode\text{{\fbmath\symbol{"1D4AF}}}\else{\fbmath\symbol{"1D4AF}}\fi}
\newunicodechar{𝓗}{\ifmmode\text{{\fbmath\symbol{"1D4D7}}}\else{\fbmath\symbol{"1D4D7}}\fi}
\newunicodechar{𝓘}{\ifmmode\text{{\fbmath\symbol{"1D4D8}}}\else{\fbmath\symbol{"1D4D8}}\fi}
\newunicodechar{𝓛}{\ifmmode\text{{\fbmath\symbol{"1D4DB}}}\else{\fbmath\symbol{"1D4DB}}\fi}
\newunicodechar{𝓜}{\ifmmode\text{{\fbmath\symbol{"1D4DC}}}\else{\fbmath\symbol{"1D4DC}}\fi}
\newunicodechar{𝓡}{\ifmmode\text{{\fbmath\symbol{"1D4E1}}}\else{\fbmath\symbol{"1D4E1}}\fi}
\newunicodechar{𝖨}{\ifmmode\text{{\fbmath\symbol{"1D5A8}}}\else{\fbmath\symbol{"1D5A8}}\fi}
\newunicodechar{ℓ}{\ifmmode\text{{\fbmath\symbol{"2113}}}\else{\fbmath\symbol{"2113}}\fi}
\newunicodechar{ℛ}{\ifmmode\text{{\fbmath\symbol{"211B}}}\else{\fbmath\symbol{"211B}}\fi}
\newunicodechar{∮}{\ifmmode\text{{\fbmath\symbol{"222E}}}\else{\fbmath\symbol{"222E}}\fi}
\newunicodechar{≪}{\ifmmode\text{{\fbmath\symbol{"226A}}}\else{\fbmath\symbol{"226A}}\fi}
\newunicodechar{≫}{\ifmmode\text{{\fbmath\symbol{"226B}}}\else{\fbmath\symbol{"226B}}\fi}
\newunicodechar{≲}{\ifmmode\text{{\fbmath\symbol{"2272}}}\else{\fbmath\symbol{"2272}}\fi}
\newunicodechar{≳}{\ifmmode\text{{\fbmath\symbol{"2273}}}\else{\fbmath\symbol{"2273}}\fi}
\newunicodechar{⋆}{\ifmmode\text{{\fbmath\symbol{"22C6}}}\else{\fbmath\symbol{"22C6}}\fi}
\newunicodechar{☉}{\ifmmode\text{{\fbsym\symbol{"2609}}}\else{\fbsym\symbol{"2609}}\fi}
\newunicodechar{☾}{\ifmmode\text{{\fbsym\symbol{"263E}}}\else{\fbsym\symbol{"263E}}\fi}
\usepackage{fvextra}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,breakanywhere,commandchars=\\\{\}}
\fvset{breaklines,breakanywhere}
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\small}
\setcounter{secnumdepth}{-2}
\usepackage{titlesec}
\titleformat{\chapter}[block]{\normalfont\huge\bfseries}{}{0pt}{}
\titlespacing*{\chapter}{0pt}{0pt}{20pt}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\nouppercase{\leftmark}}
\fancyfoot[C]{\thepage}
\renewcommand{\chaptermark}[1]{\markboth{#1}{}}
\renewcommand{\headrulewidth}{0.3pt}
\fancypagestyle{plain}{\fancyhf{}\fancyfoot[C]{\thepage}\renewcommand{\headrulewidth}{0pt}}
\emergencystretch=3em
\setlength{\parskip}{4pt plus 1pt}
\setlength{\parindent}{0pt}
\newsavebox{\fitbox}
\newcommand{\fitdisplay}[1]{\sbox{\fitbox}{$#1$}\ifdim\wd\fitbox>\linewidth\resizebox{\linewidth}{!}{\usebox{\fitbox}}\else\usebox{\fitbox}\fi}
\makeatletter
\def\maxwidth{\ifdim\Gin@nat@width>\linewidth\linewidth\else\Gin@nat@width\fi}
\makeatother
"""
open(f"{OUT}/header.tex", "w").write(header)

cmd = ["pandoc", f"{OUT}/skripta.md", "-f", "markdown-implicit_figures-raw_html+tex_math_dollars",
       "-o", f"{OUT}/skripta.pdf",
       "--pdf-engine=xelatex", "--toc", "--toc-depth=2", "--top-level-division=chapter",
       "--include-in-header", f"{OUT}/header.tex",
       "-V", "documentclass=report", "-V", "papersize=a4", "-V", "geometry:margin=2.2cm",
       "-V", "mainfont=DejaVu Serif", "-V", "sansfont=DejaVu Sans", "-V", "monofont=DejaVu Sans Mono",
       "-V", "monofontoptions=Scale=0.82", "-V", "mathfont=Latin Modern Math", "-V", "fontsize=10pt",
       "-V", "colorlinks=true", "-V", "linkcolor=black", "-V", "urlcolor=blue", "-V", "toc-title=Obsah"]
cmd[cmd.index("-o") + 1] = f"{OUT}/skripta.tex"
cmd = [c for c in cmd if c != "--pdf-engine=xelatex"] + ["-s"]
r = subprocess.run(cmd, capture_output=True, text=True)
print("pandoc", r.returncode, r.stderr[-3000:])

tex = open(f"{OUT}/skripta.tex", encoding="utf-8").read()
# display math: shrink to the line width when too wide
tex = re.sub(r"\\\[(.*?)\\\]",
             lambda m: "\\[\\fitdisplay{\\displaystyle " + m.group(1).strip() + "}\\]", tex, flags=re.S)


# inline code: allow breaks after / _ . and -
def brk(m):
    s = m.group(1)
    s = re.sub(r"(\\_|/|\.|-|=|,)", r"\1\\allowbreak{}", s)
    return "\\texttt{" + s + "}"


tex = re.sub(r"\\texttt\{((?:[^{}]|\{[^{}]*\})*)\}", brk, tex)
open(f"{OUT}/skripta.tex", "w", encoding="utf-8").write(tex)
for i in range(3):
    r = subprocess.run(["xelatex", "-interaction=nonstopmode", "skripta.tex"], cwd=OUT,
                       capture_output=True, text=True)
    print("xelatex", i, r.returncode)
shutil.copy(f"{OUT}/skripta.pdf", PDF)
print("built", PDF.relative_to(ROOT))
