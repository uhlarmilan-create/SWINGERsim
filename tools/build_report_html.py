"""Interactive HTML edition of a report in docs/<lang>/ (one self-contained file).

    python tools/build_report_html.py docs/en/squadra_test3.md [docs/cs/squadra_test3.md ...]

The Markdown is converted with pandoc. Every figure whose PNG has a plotly JSON next to it (made by the report
scripts, e.g. tools/squadra3_report.py) is replaced by the interactive plotly figure: zoom, pan, hover values,
click the legend to hide a curve. Other images are embedded as data URIs. plotly.js is inlined from the plotly
Python package, so the file works offline. Output: docs/html/<name>-<lang>.html.
"""

from __future__ import annotations

import argparse
import base64
import html
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "html"

CSS = """
:root { --bg:#ffffff; --fg:#1d2126; --muted:#5b6470; --rule:#e3e6ea; --accent:#1F4FB0; --code:#f4f5f7; }
html { background: var(--bg); }
body { background: var(--bg); color: var(--fg); font: 16px/1.6 "DejaVu Serif", Georgia, serif; margin: 0; }
main { max-width: 1060px; margin: 0 auto; padding: 24px 16px 80px; }
h1 { font: 700 30px/1.25 "DejaVu Sans", Arial, sans-serif; margin: 12px 0 18px; }
h2 { font: 700 23px/1.3 "DejaVu Sans", Arial, sans-serif; margin: 42px 0 10px; padding-top: 10px;
     border-top: 1px solid var(--rule); }
h3 { font: 700 18px/1.35 "DejaVu Sans", Arial, sans-serif; margin: 28px 0 8px; }
a { color: var(--accent); }
table { border-collapse: collapse; margin: 14px 0 20px; font: 14px/1.45 "DejaVu Sans", Arial, sans-serif;
        display: block; overflow-x: auto; max-width: 100%; }
th, td { border-bottom: 1px solid var(--rule); padding: 6px 10px; text-align: left; vertical-align: top; }
th { background: var(--code); }
code { background: var(--code); padding: 1px 4px; border-radius: 3px; font-size: 0.9em; }
figure { margin: 22px 0 6px; }
figcaption { font: 14px/1.45 "DejaVu Sans", Arial, sans-serif; color: var(--muted); margin-top: 6px; }
.plot { width: 100%; }
.hint { font: 13px/1.4 "DejaVu Sans", Arial, sans-serif; color: var(--muted); background: var(--code);
        padding: 8px 12px; border-radius: 6px; margin: 14px 0; }
img { max-width: 100%; }
"""


def _pandoc(md: str) -> str:
    r = subprocess.run(
        ["pandoc", "-f", "markdown", "-t", "html5", "--mathjax"],
        input=md,
        capture_output=True,
        text=True,
        check=True,
    )
    return r.stdout


def build(md_path: Path) -> Path:
    md_path = md_path.resolve()
    lang = md_path.parent.name
    text = md_path.read_text()
    title = re.search(r"^# (.+)$", text, re.MULTILINE).group(1)
    figs: list[tuple[str, str, Path]] = []

    def fig(m):
        cap, src = m.group(1), m.group(2)
        p = (md_path.parent / src).resolve()
        figs.append((cap, src, p))
        return f"\n\nFIGURE_SLOT_{len(figs) - 1}\n\n"

    body = _pandoc(re.sub(r"!\[(.*?)\]\((.*?)\)", fig, text))
    plots = []
    for k, (cap, _src, p) in enumerate(figs):
        cap_html = _pandoc(cap).strip()
        cap_html = re.sub(r"^<p>|</p>$", "", cap_html)
        j = p.with_suffix(".json")
        if j.exists():
            plots.append((f"plot{k}", j.read_text()))
            block = f'<figure><div id="plot{k}" class="plot"></div><figcaption>{cap_html}</figcaption></figure>'
        elif p.exists():
            b64 = base64.b64encode(p.read_bytes()).decode()
            block = f'<figure><img alt="" src="data:image/png;base64,{b64}"><figcaption>{cap_html}</figcaption></figure>'
        else:
            block = f"<figure><figcaption>{cap_html} (figure missing)</figcaption></figure>"
        body = re.sub(rf"<p>FIGURE_SLOT_{k}</p>", lambda _m, b=block: b, body)
    # links between the language editions point to the HTML editions
    body = re.sub(
        r'href="\.\./(en|cs)/([\w_]+)\.md"',
        lambda m: f'href="{m.group(2)}-{m.group(1)}.html"',
        body,
    )
    import plotly

    js = (Path(plotly.__file__).parent / "package_data" / "plotly.min.js").read_text()
    hint = "The figures are interactive: drag to zoom, double-click to reset, click a legend entry to hide a curve, hover to read values."
    if lang == "cs":
        hint = "Grafy jsou interaktivní: tažením myší zvětšíte výřez, dvojklik vrátí původní pohled, kliknutím na položku legendy křivku skryjete, po najetí myší se zobrazí hodnoty."
    init = "\n".join(
        f"(function(){{var f={fj};f.layout.width=undefined;f.layout.autosize=true;"
        f"Plotly.newPlot('{pid}',f.data,f.layout,{{responsive:true,displaylogo:false}});}})();"
        for pid, fj in plots
    )
    page = f"""<!doctype html>
<html lang="{lang}"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title><style>{CSS}</style><script>{js}</script></head>
<body><main><div class="hint">{hint}</div>
{body}
</main><script>{init}</script></body></html>
"""
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{md_path.stem}-{lang}.html"
    out.write_text(page)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", type=Path)
    a = ap.parse_args(argv)
    for f in a.files:
        print(build(f))
    return 0


if __name__ == "__main__":
    sys.exit(main())
