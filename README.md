# SWINGERsim

**SWINGERsim** is a simulator of hierarchical multiple stars: triples, 2+2 and 3+1 quadruples
and nested hierarchies of up to ten stars (2+(2+1), 2+2+2, 3+3, …). It follows a system from
the pre-main sequence to the end of the asymptotic giant branch and computes what an observer
of eclipsing multiples measures.

- **Dynamics.** Direct N-body integration (REBOUND IAS15 with REBOUNDx: general relativity,
  tides and spins, mass change) and its own secular (orbit-averaged) equations for 2+2 and
  for any nested hierarchy, with Kozai–Lidov cycles, precession, resonances, stability criteria
  and the fate of the system (stable, merger, ejection, break-up).
- **Stars.** MIST and SSE evolutionary tracks, apsidal constants and gyration radii (Claret),
  equilibrium and dynamical tides (Hurley et al. 2002), inertial waves, winds, stable and
  tertiary mass transfer, common envelopes, magnetic braking, mergers, pulsations, discs.
- **Observables.** Light curves with Roche geometry, times of minima and O−C (light-time and
  dynamical effects), radial velocities, eclipse windows, changes of the inclination and of
  eclipse depths.
- **Inference.** Fits of times of minima, radial velocities, light curves and inclination
  histories with one N-body model of the whole system (least squares and MCMC), Monte Carlo
  ensembles conditioned on observations (importance weights, Bayesian evidence), emulators and
  simulation-based inference.
- **Web application** (`swingersim gui`), a job **worker** for a separate computer and
  self-contained HTML reports of every result.

Much of the physics, all of the secular and fitting machinery, the observables and the
application are SWINGERsim's own code; the external libraries it builds on (REBOUND,
REBOUNDx, SciPy, emcee, optionally PHOEBE, MESA) are listed in the User Guide, chapter
*How SWINGERsim works*.

## Installation

Linux, macOS, or Windows through WSL (Windows Subsystem for Linux), with Python 3.11 or newer:

```bash
git clone https://github.com/uhlarmilan-create/SWINGERsim.git
cd SWINGERsim
./install.sh
source .venv/bin/activate
```

`install.sh` creates a virtual environment `.venv`, installs SWINGERsim with its common optional
packages, makes a test run and runs the fast benchmarks. Options: `--minimal` (core only),
`--with-sbi` (simulation-based inference, about 1 GB more), `--with-phoebe` (PHOEBE 2 light
curves), `--venv DIR`, `--python PY`, `--no-check`; `./install.sh --help` lists them.

**Windows.** Install WSL once in PowerShell (as administrator) with `wsl --install`, open the
Ubuntu terminal and follow the steps above inside it. `install.ps1` does the same from
PowerShell: it checks WSL and runs `install.sh` inside it.

Manual installation: `python3 -m venv .venv && source .venv/bin/activate && pip install -e ".[all,dev]"`.

## First steps

```bash
swingersim gui                                                  # web application in the browser
swingersim run configs/examples/example_2p2.yaml -o output/example_2p2.h5
swingersim ensemble configs/examples/example_v994_ensemble.yaml -o output/v994_ens -j 4
swingersim ensemble-report output/v994_ens
swingersim fit configs/examples/example_fit.yaml
swingersim validate --all                                       # all benchmarks against published results
```

A run writes next to its HDF5 file an HTML report, a summary, CSV time series, an interactive 3D
viewer and PNG figures. `swingersim --help` lists all commands and `swingersim <command> --help`
their options.

## Documentation

All documents exist in English (`-en`) and Czech (`-cs`) as PDF in [`docs/pdf/`](docs/pdf) and
as Markdown in [`docs/en/`](docs/en) and [`docs/cs/`](docs/cs).

| Document | Content |
| --- | --- |
| [User Guide](docs/pdf/SWINGERsim-guide-en.pdf) | purpose, how it works, installation, quick start, web application, command line, workflows, examples, configuration reference, outputs, troubleshooting, glossary |
| [Methods](docs/pdf/SWINGERsim-methods-en.pdf) | the physics and numerics: equations, approximations, validity and references |
| [Testing, Part I](docs/pdf/SWINGERsim-part1-en.pdf) | validation against published results and independent codes |
| [Testing, Part II](docs/pdf/SWINGERsim-part2-en.pdf) | fits of observed eclipsing multiples |
| [Testing, Part III](docs/pdf/SWINGERsim-part3-en.pdf) | changing inclinations and disappearing eclipses (V907 Sco, QX Cas, AH Cep) |
| [Testing, Part IV](docs/pdf/SWINGERsim-part4-en.pdf) | special cases: the future of CzeV343, the hierarchy of CzeV4315, MCMC convergence |
| [Open questions and outlook](docs/pdf/SWINGERsim-outlook-en.pdf) | what the tests left open and possible extensions |

**Textbook (Czech).** *SWINGERsim: fyzika a numerické modelování vícenásobných hvězd* – the physics and numerical
modelling of multiple stars, from the two-body problem to the full workflow on V907 Sco, with solved examples and
exercises: [PDF](docs/pdf/SWINGERsim-skripta-cs.pdf), sources in [`docs/skripta/cs/`](docs/skripta/cs)
(`python tools/build_skripta_pdf.py`).

## Tests

```bash
pytest                                  # unit tests (fast; slow ones are deselected)
pytest -m slow                          # long integrations
swingersim validate --all               # benchmarks against published results
```

## Repository layout

| Folder | Content |
| --- | --- |
| `src/swingersim/` | the package, with its stellar tracks and tables in `src/swingersim/data/` |
| `configs/` | example configurations and the configurations of the tests |
| `data/observed/` | observed data of the test systems (times of minima, radial velocities, light-curve extracts) |
| `docs/` | the documentation (Markdown sources, figures, PDF) |
| `tools/` | scripts that make the figures and numbers of the test reports and build the PDFs |
| `tests/` | the test suite |

Results (`output/`) and job folders (`jobs/`) are created locally and are not part of the
repository.

## Citing

If you use SWINGERsim in a publication, please cite it as given in [`CITATION.cff`](CITATION.cff)
and cite the external codes and data it used (REBOUND, REBOUNDx, MIST, Claret tables, …) as
listed in the Methods document.

## License

MIT – see [`LICENSE`](LICENSE). The stellar tracks and tables in `src/swingersim/data/` come from
published data sets; cite their original papers (Methods document, References).
