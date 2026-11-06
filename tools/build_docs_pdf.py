"""Build the SWINGERsim documentation (PDF) from docs/<lang>/<document>/*.md.

Requires pandoc and XeLaTeX:
    sudo apt install pandoc texlive-xetex texlive-fonts-recommended texlive-lang-czechslovak

    python tools/build_docs_pdf.py                       # every document, both languages
    python tools/build_docs_pdf.py guide part3           # some documents
    python tools/build_docs_pdf.py --lang cs methods     # one language

Documents (one folder each, chapters in file-name order):
    guide    User Guide              methods  Physics and Methods
    part1    Testing, Part I         part2    Testing, Part II
    part3    Testing, Part III       part4    Testing, Part IV
    outlook  Open questions and future work

Each document starts with a YAML block in its first chapter (title, subtitle). Inside a document
links are [text](#anchor) or [text](chapter.md#anchor); figures are ../../figures/<path>.
Output: docs/pdf/SWINGERsim-<document>-<lang>.pdf
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
OUT = DOCS / "pdf"

# document: (book layout with chapters?, title of the table of contents per language)
DOCUMENTS = {
    "guide": True,
    "methods": True,
    "part1": True,
    "part2": False,
    "part3": False,
    "part4": False,
    "outlook": False,
}
TOC = {"en": "Contents", "cs": "Obsah"}
LANG = {"en": "en-GB", "cs": "cs-CZ"}

FIGURE = re.compile(r"\]\((?:\.\./)+figures/([^)\s]+)\)")
CHAPTER_LINK = re.compile(r"\]\([0-9a-z_]+\.md(#[^)]+)?\)")


def heading_slug(title: str) -> str:
    """Anchor of a heading: lower case, accents dropped, dashes (also en dashes and minus signs)
    kept as hyphens, other punctuation removed, spaces as hyphens."""
    import unicodedata

    t = re.sub("[‐-―−]", "-", title.lower())
    t = unicodedata.normalize("NFKD", t).encode("ascii", "ignore").decode()
    t = re.sub(r"[^a-z0-9 _-]", " ", t)
    return re.sub(r"-+", "-", re.sub(r"\s+", "-", t.strip())).strip("-")


def _cells(line: str) -> list[str]:
    return [c.strip() for c in re.split(r"(?<!\\)\|", line.strip())[1:-1]]


def proportional_tables(text: str) -> str:
    """Give pipe-table columns widths proportional to their content (pandoc reads the dashes)."""
    lines = text.split("\n")
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        nxt = lines[i + 1] if i + 1 < len(lines) else ""
        if line.startswith("|") and re.match(r"^\|(\s*:?-{3,}:?\s*\|)+\s*$", nxt):
            j = i + 2
            while j < len(lines) and lines[j].startswith("|"):
                j += 1
            rows = [_cells(r) for r in [line, *lines[i + 2 : j]]]
            ncol = len(rows[0])
            weight = [
                max(
                    4.0,
                    sum(min(len(r[k]), 80) for r in rows if k < len(r)) / len(rows),
                    1.2
                    * max((len(w) for r in rows if k < len(r) for w in r[k].split()), default=3),
                )
                for k in range(ncol)
            ]
            aligns = [
                (":" if c.strip().startswith(":") else "", ":" if c.strip().endswith(":") else "")
                for c in _cells(nxt)
            ]
            dashes = [max(3, round(100 * w / sum(weight))) for w in weight]
            sep = [
                a + "-" * (d - len(a) - len(b)) + b
                for d, (a, b) in zip(dashes, aligns + [("", "")] * ncol, strict=False)
            ]
            out += [line, "|" + "|".join(sep) + "|", *lines[i + 2 : j]]
            i = j
            continue
        out.append(line)
        i += 1
    return "\n".join(out)


def chapters(lang: str, doc: str) -> list[Path]:
    return sorted((DOCS / lang / doc).glob("*.md"))


def assemble(lang: str, doc: str) -> str:
    """The chapters of a document joined; links between its chapters point at anchors, figure
    paths become absolute; the YAML block of the first chapter stays at the top."""
    parts = []
    for path in chapters(lang, doc):
        text = path.read_text(encoding="utf-8")
        text = CHAPTER_LINK.sub(lambda m: f"]({m.group(1) or '#'})", text)
        text = FIGURE.sub(lambda m: f"]({(DOCS / 'figures' / m.group(1)).as_posix()})", text)
        parts.append(proportional_tables(text).strip())
    return "\n\n".join(parts) + "\n"


HEADER = r"""
\usepackage{fvextra}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
\usepackage{float}
\floatplacement{figure}{H}
\renewcommand{\arraystretch}{1.25}
\usepackage{newunicodechar}
\newunicodechar{☉}{\ensuremath{_\odot}}
\newunicodechar{⊙}{\ensuremath{_\odot}}
\newunicodechar{≳}{\ensuremath{\gtrsim}}
\newunicodechar{≲}{\ensuremath{\lesssim}}
\newunicodechar{≪}{\ensuremath{\ll}}
\newunicodechar{≫}{\ensuremath{\gg}}
\newunicodechar{ℓ}{\ensuremath{\ell}}
\usepackage{amssymb}
\usepackage{caption}
\captionsetup[figure]{labelformat=empty,font=small}
\captionsetup[table]{labelformat=empty,font=small}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\textsc{SWINGERsim}}
\fancyhead[R]{\small\nouppercase{\leftmark}}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.3pt}
"""


def build(lang: str, doc: str) -> Path:
    if not chapters(lang, doc):
        raise FileNotFoundError(f"no chapters in docs/{lang}/{doc}/")
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"SWINGERsim-{doc}-{lang}.pdf"
    book = DOCUMENTS[doc]
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / f"{doc}-{lang}.md"
        src.write_text(assemble(lang, doc), encoding="utf-8")
        hdr = Path(tmp) / "header.tex"
        hdr.write_text(HEADER, encoding="utf-8")
        cmd = [
            "pandoc",
            str(src),
            "-o",
            str(out),
            "--pdf-engine=xelatex",
            "--toc",
            "--toc-depth=2",
            "--number-sections",
            "--include-in-header",
            str(hdr),
            "-V",
            f"documentclass={'report' if book else 'article'}",
            *(["--top-level-division=chapter"] if book else []),
            "-V",
            "papersize=a4",
            "-V",
            "geometry:margin=2.3cm",
            "-V",
            "mainfont=DejaVu Serif",
            "-V",
            "sansfont=DejaVu Sans",
            "-V",
            "monofont=DejaVu Sans Mono",
            "-V",
            "monofontoptions=Scale=0.85",
            "-V",
            "fontsize=10pt",
            "-V",
            "colorlinks=true",
            "-V",
            "linkcolor=blue",
            "-V",
            "urlcolor=blue",
            "-V",
            f"toc-title={TOC[lang]}",
            "-M",
            f"date={dt.datetime.now(dt.UTC).date().isoformat()}",
            "-M",
            f"lang={LANG[lang]}",
        ]
        subprocess.run(cmd, check=True)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("documents", nargs="*", help="documents to build (default: all)")
    p.add_argument("--lang", choices=sorted(TOC), action="append")
    args = p.parse_args(argv)
    if shutil.which("pandoc") is None or shutil.which("xelatex") is None:
        print("pandoc and xelatex are required, see the docstring of this script", file=sys.stderr)
        return 1
    bad = [d for d in args.documents if d not in DOCUMENTS]
    if bad:
        p.error(f"unknown document(s) {bad}; choose from {list(DOCUMENTS)}")
    for lang in args.lang or sorted(TOC):
        for doc in args.documents or DOCUMENTS:
            if not chapters(lang, doc):
                print(f"skipped {doc} ({lang}): no chapters yet")
                continue
            print(f"built {build(lang, doc).relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
