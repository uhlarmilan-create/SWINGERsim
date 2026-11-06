# SWINGERsim documentation – style guide (for every author of a chapter)

The documentation is being rewritten from scratch into seven documents, each in English (docs/en/)
and Czech (docs/cs/). The English text is written first; the Czech one is a faithful translation.

## The seven documents (one folder each, chapters in file-name order: 01_xxx.md, 02_xxx.md …)

| folder   | title (EN)                                                         | title (CS)                                                          |
|----------|--------------------------------------------------------------------|---------------------------------------------------------------------|
| guide    | SWINGERsim – User Guide                                            | SWINGERsim – Uživatelská příručka                                   |
| methods  | SWINGERsim – Physics and Methods                                   | SWINGERsim – Fyzika a metody                                        |
| part1    | Testing SWINGERsim, Part I: Validation against the literature      | Testování SWINGERsimu, část I: Ověření proti literatuře             |
| part2    | Testing SWINGERsim, Part II: Three published multiple stars        | Testování SWINGERsimu, část II: Tři publikované vícenásobné hvězdy  |
| part3    | Testing SWINGERsim, Part III: Eclipsing binaries with a changing inclination | Testování SWINGERsimu, část III: Zákrytové dvojhvězdy s proměnným sklonem |
| part4    | Testing SWINGERsim, Part IV: Special cases                         | Testování SWINGERsimu, část IV: Zvláštní případy                    |
| outlook  | SWINGERsim – Open questions and future work                        | SWINGERsim – Otevřené otázky a další práce                          |

The first chapter of each document starts with a YAML block:

    ---
    title: "SWINGERsim"
    subtitle: "User Guide"
    ---

Build: `python tools/build_docs_pdf.py <folder> --lang en` → docs/pdf/SWINGERsim-<folder>-en.pdf.
In guide, methods and part1 a level-1 heading (`#`) is a chapter; in the other documents it is a
section. Every chapter file starts with one `#` heading.

## Hard rules

0. **Own work vs external tools.** SWINGERsim is not a wrapper: much of it was programmed for it.
   Wherever a method is described, say what is SWINGERsim's own implementation and what comes from
   an external library or data set (REBOUND, REBOUNDx, emcee, MIST/SSE tracks, Claret tables …).

1. **Name.** The software is **SWINGERsim** (the name is final; the code rename quadsim → swingersim
   follows right after the documentation). Write `swingersim` for the command, the Python package
   (`import swingersim`) and paths (`src/swingersim/...`). Never write "quadsim".
2. **No versions, no history.** SWINGERsim is version 1.0 and development is closed. Never mention
   version numbers (v0.27, "since version 0.31", "0.36.0"), job numbers ("job 062", "jobs 069–085"),
   campaign group numbers ("group 6"), dates of development steps, or "we added / fixed / changed".
   Describe what the software **is** and **does**. Problems met during testing may be described as
   findings ("a member that is unstable at the start can hang the integrator; therefore each member
   has a wall-clock limit"), never as a changelog.
3. **Audience.** Astronomers and astrophysics students, not programmers. Tone: an astrophysics
   student presenting results to professors. Every technical term is defined at first use (in the
   text or in the glossary); every claim has its basis (equation, paper, test) and its limits.
   Anticipate the obvious question a referee would ask and answer it.
4. **Physics is validated against the literature.** Cite papers as (Author Year) or (Author et al.
   Year) in the text, with the full reference in the document's reference list (ADS-style: authors,
   year, journal, volume, page; arXiv id if useful). Do not invent references; use those already in
   docs/en/references.md and the source pages, or the papers in the code docstrings.
5. **Numbers** come only from the sources named in your brief (existing docs, output files, project
   notes). Do not invent or round away numbers; if two sources disagree, use the later/final one
   and say nothing about the earlier one.
6. **Units and notation.** M☉, R☉, L☉, au, d, yr, Myr, Gyr, K, deg (°). Minus sign "−" in O−C.
   Inclination i, mutual inclination i_mut (or $i_\mathrm{mut}$), eccentricity e, period P.
   Math in LaTeX `$…$`. Use "ZLK" (von Zeipel–Lidov–Kozai) after defining it once.
7. **Figures.** `![**Figure N.** Caption: what is shown, what to look at, what it means.](../../figures/<path>.png)`
   — numbered by hand per document (Figure 1, 2, …), every figure referred to in the text.
   Use existing figures under docs/figures/; new figures go under docs/figures/<document>/.
   Tables: pipe tables with a caption line above, `**Table N.** …`.
   **Every figure, graph and table is explained in detail**: the caption says what is plotted
   (axes, units, colours, symbols, panels) and what to look at; the text has a paragraph that
   walks the reader through it and says what it shows (and, in the test parts, how it compares
   with the published result). For tables: what each column means and how to read a row.
8. **Links.** Inside a document: `[text](#anchor)` with explicit anchors `# Heading {#anchor}`.
   Between documents: plain text ("see *Physics and Methods*, chapter Tides and rotation").
9. **Czech version.** Czech (not Slovak), decimal comma in the text and tables (0,64; 2,4 M☉), the
   same figures, the same anchors (anchors stay in English), the same structure. Terms: zákrytová
   dvojhvězda, sklon, vzájemný sklon, precese uzlů, apsidální stáčení, slapy, Rocheův lalok,
   přenos hmoty, ansámbl (ensemble), věrohodnost, apriorní rozdělení (prior), aposteriorní
   rozdělení (posterior), okamžik minima, diagram O−C, třetí světlo, efekt světelného času (LTTE).
10. **Style.** Short paragraphs, active voice, present tense. Lists only when the reader scans or
    compares. Tables for parameters and results. No marketing words ("powerful", "cutting-edge").
    Commands and configuration keys in `code`; complete YAML/shell examples in fenced blocks
    (```yaml, ```bash) that actually work with the current code.

## The command line after the rename

`swingersim run`, `swingersim ensemble`, `swingersim gui`, … (the same sub-commands as
`quadsim --help` lists today; check `src/quadsim/cli.py`). Configuration files are YAML as today.
