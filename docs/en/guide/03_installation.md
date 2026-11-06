# Installation {#installation}

SWINGERsim is a Python package that runs on Linux and macOS (on Windows through the Windows
Subsystem for Linux, WSL). This chapter lists what it needs, how to install it and check the
installation, the optional parts, the worker for long computations and how to build this
documentation.

## Requirements {#requirements}

| Item | Requirement | Why |
| --- | --- | --- |
| Python | 3.11 or newer | the language of the package |
| C compiler | `gcc` or `clang` | REBOUND and REBOUNDx, the N-body libraries, are compiled during installation where pip finds no ready-made package |
| operating system | Linux or macOS; Windows through WSL | the worker uses process groups and signals of POSIX systems |
| disk | about 100 MB for the stellar tracks and tables of the package, plus the optional libraries (the `sbi` group alone about 1 GB) | the data are part of the package |

**Table 7.** Requirements. The left column names the item, the middle the requirement, the
right column why it is needed. Memory depends on the task: a parallel ensemble needs memory for
every process, which is why the worker watches the free memory (below).

## Steps {#install-steps}

```bash
git clone https://github.com/uhlarmilan-create/SWINGERsim.git
cd SWINGERsim
./install.sh
source .venv/bin/activate
```

The installation script `install.sh` finds Python 3.11 or newer, checks that the stellar
tracks and tables of the package are present, creates an isolated Python environment in the
folder `.venv` (so that the packages of SWINGERsim do not mix with the rest of the system),
installs SWINGERsim with the groups `all` and `dev` of Table 8, makes a test run and runs the
fast benchmarks. It never uses `sudo`: when a system package is missing (Python, its `venv`
module, a C compiler), it stops and prints the command that installs it. Its options:

| Option | Effect |
| --- | --- |
| `--minimal` | only the core and the test suite (groups `dev`) |
| `--with-sbi` | also the group `sbi` (PyTorch, about 1 GB) |
| `--with-phoebe` | also the group `phoebe` |
| `--venv DIR` | the environment in `DIR` instead of `.venv` |
| `--python PY` | this Python interpreter |
| `--no-check` | no test run and no benchmarks at the end |

**Table 7a.** Options of `install.sh` (`./install.sh --help` lists them). They can be combined,
e.g. `./install.sh --with-sbi --no-check`.

Activate the environment with `source .venv/bin/activate` in every new terminal; the command
`swingersim` is then available. **Windows:** install WSL once (PowerShell as administrator:
`wsl --install`, restart, open *Ubuntu* once to create a user) and run the steps above in the
Ubuntu terminal; alternatively `powershell -ExecutionPolicy Bypass -File install.ps1` checks
WSL and runs `install.sh` inside it. The web application started in WSL opens in the Windows
browser at the printed address.

The same by hand:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[all,dev]"
swingersim --help
```

`pip install -e` installs the package in editable mode: the source in `src/swingersim/` is
used directly, so the configurations, data and outputs in the folder of the repository are
found. The names in square brackets are optional groups of dependencies (Table 8). The last
command lists the sub-commands; if it prints them, the installation works.

| Group | Packages | What it enables |
| --- | --- | --- |
| (always) | numpy, scipy, pydantic, pyyaml, h5py, rebound, reboundx | runs, all physics, HDF5 output |
| `dev` | pytest, ruff | the test suite |
| `all` | the groups `viz`, `ensemble`, `fast`, `emulator`, `fit` | everything except `sbi` and `phoebe`; installed by `install.sh` |
| `viz` | plotly, kaleido, matplotlib | the interactive 3D viewer, the plots inside the HTML reports, PNG figures (matplotlib; kaleido with a Chrome browser for PNG files of plotly figures); the web application uses plotly if present, otherwise it loads plotly.js from the internet |
| `fit` | emcee, corner, astropy | MCMC errors of fits, corner plots, conversion of HJD/JD (UTC) to BJD_TDB |
| `fast` | numba | the compiled secular equations, about a hundred times faster secular runs and ensembles |
| `ensemble` | pandas, pyarrow | ensemble summaries also as Parquet files |
| `emulator` | scikit-learn | neural emulators (`emulator-train`, `emulator-report`) and the stability classifiers below |
| `sbi` | sbi, torch | simulation-based inference (`swingersim sbi`); about 1 GB |
| `phoebe` | phoebe | PHOEBE 2 as the light-curve engine of the inner binaries (`--engine phoebe`) |

**Table 8.** Optional groups of dependencies. Install a group with
`pip install -e ".[<group>]"`; several groups are separated by commas. The right column says
which commands or options need the group; without it, they stop with a message naming the
missing package.

## Check the installation {#install-check}

```bash
pytest -m "not slow and not validation"     # unit tests, a few minutes
swingersim validate                          # fast benchmarks against published results
swingersim run configs/examples/example_2p2.yaml -o output/check/example.h5
```

The unit tests check the conversions between elements and positions, the hierarchical
initialisation, energy conservation, the relativistic periastron advance, the configuration
schema, the web application and the documentation. `swingersim validate` reproduces a set of
published results and prints PASS, EXPLAINED (a difference whose reason is documented) or FAIL
for each; the benchmarks and their verdicts are discussed in *Testing, Part I*. The last command
is the first run of the [Quick start](#quick-start).

## Optional: stability classifiers of Vynatheya et al. (2023) {#install-vynatheya}

The machine-learning stability classifiers for 2+2 and 3+1 quadruples (Vynatheya et al. 2023)
are published by their authors and are not part of SWINGERsim. To use them, download their
repository into `data/`:

```bash
git clone https://github.com/pavanvyn/quadruple-stability data/quadruple-stability
pip install -e ".[emulator]"
```

`swingersim stability config.yaml` then reports their verdict for every sub-quadruple. The
networks were saved with scikit-learn 1.2.2; newer versions read them, and the version warning
is suppressed. Another location can be given with the environment variable
`SWINGERSIM_VYNATHEYA_DIR`.

## Optional: MSE for independent comparisons {#install-mse}

The Multiple Stellar Evolution code MSE (Hamers et al. 2021) is an independent code for
hierarchical multiples that was used to cross-check SWINGERsim (*Testing, Part I*). It is not
needed to run SWINGERsim. To make own comparisons, install it from its repository; it needs C,
C++ and Fortran compilers:

```bash
sudo apt install gfortran
git clone https://github.com/hamers/mse ~/mse
cd ~/mse && make
```

## Optional: MESA for own stellar-model grids {#install-mesa}

SWINGERsim can compute its own grids of stellar models with MESA (Modules for Experiments in
Stellar Astrophysics; Paxton et al. 2011; Jermyn et al. 2023): pre-main-sequence and envelope
structure, $k_2$ and the dissipation of inertial waves from full stellar profiles. This is
needed only when the built-in MIST tracks and the Claret tables do not cover a star (an unusual
metallicity, a particular physics set-up). MESA is installed once by hand; SWINGERsim then runs
it (also through the worker).

Requirements: Linux x86_64, at least 8 GB memory, about 30 GB of free disk (MESA about 20 GB,
its software development kit about 2 GB, a grid a few GB) and one to two hours.

1. System packages (the MESA software development kit, SDK, brings its own compilers):

   ```bash
   sudo apt install binutils make perl libx11-6 libx11-dev zlib1g zlib1g-dev tcsh unzip python3
   ```

2. Download the MESA SDK for Linux x86_64 from
   <http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk> and unpack it into `~/mesa`.
3. Download a MESA release (the Zenodo record <https://doi.org/10.5281/zenodo.2602941> or
   <https://github.com/MESAHub/mesa/releases>) and unpack it into `~/mesa`. The SDK must be
   the version the MESA release asks for.
4. Add to `~/.bashrc` and open a new terminal:

   ```bash
   export MESASDK_ROOT=$HOME/mesa/mesasdk
   source $MESASDK_ROOT/bin/mesasdk_init.sh
   export MESA_DIR=$HOME/mesa/mesa-<release>
   export OMP_NUM_THREADS=2
   ```

5. Compile MESA (30–60 minutes): `cd $MESA_DIR && ./install`.
6. Test with SWINGERsim – one 1 M☉ track from the pre-main sequence to the end of the main
   sequence:

   ```bash
   swingersim mesa-grid configs/mesa/lowmass_pms_ms.yaml --test
   ```

7. Compute the whole grid and reduce it (hours; both steps can run as worker jobs):

   ```bash
   swingersim mesa-grid configs/mesa/lowmass_pms_ms.yaml -j 6
   swingersim mesa-reduce configs/mesa/lowmass_pms_ms.yaml
   ```

The reduced grid is `mesa/lowmass_pms_ms/lowmass_pms_ms.npz`; a configuration uses it with
`stellar_model.structure: mesa:lowmass_pms_ms`.

## The worker for long computations {#install-worker}

Large ensembles (tens of thousands of members, N-body ensembles), fits of many free parameters
and MESA grids take hours to days. They are best computed on a computer with many cores that
stays switched on. `swingersim worker` turns such a computer into a job runner: it watches a
folder, runs every job file that appears there and writes the results into `output/`.

```bash
cd swingersim
source .venv/bin/activate
swingersim worker                  # watches jobs/, Ctrl+C stops it
swingersim worker --cores 8        # use at most 8 cores for all jobs together
```

The worker creates the folders of Table 9 and writes a heartbeat file `jobs/worker.json` (name
of the computer, number of cores, the versions it was started with, the current job, the time of
its last check) every few seconds. Jobs that need one core (a run, a light curve, eclipse
times) run side by side up to the number of cores; a job that needs all cores (an ensemble, a
fit or a chaos map with `processes: 0`, `validate --all`) waits for the others and runs alone.

| Folder or file | Content |
| --- | --- |
| `jobs/drafts/` | configurations written in the web application (not run by the worker) |
| `jobs/submitted/` | frozen copies of the configurations of submitted jobs |
| `jobs/inbox/` | `<id>.yaml`: queued jobs, taken in the order of their names; `<id>.yaml.hold`: jobs waiting until they are queued |
| `jobs/running/` | the jobs being computed, each with `<id>.yaml`, its screen output `<id>.log` and its progress `<id>.status.json` |
| `jobs/done/`, `jobs/failed/` | finished and failed jobs with the same three files |
| `jobs/cancel/` | an empty file `<id>` stops that job; its output inside `output/` is deleted and its files are kept here |
| `jobs/pause/` | an empty file `<id>` pauses the running job; removing the file lets it continue exactly where it stopped |
| `jobs/stop` | an empty file stops the worker; running jobs go back to the inbox |
| `jobs/worker.json`, `jobs/worker.lock` | heartbeat; lock that prevents a second worker on the same folder |

**Table 9.** The job folder of the worker. The left column is the folder or file inside the
SWINGERsim folder, the right column what it holds. The web application uses the same layout.

To keep the worker running after logging out, start it in `tmux` or as a user service, e.g.
`systemd-run --user --unit=swingersim-worker swingersim worker --root ~/swingersim`. When the
free memory of the computer falls below 8 % of the total (at least 1.5 GB), the worker stops
the job that uses the most memory and marks it failed with an explanation, so that the computer
stays usable. After an update of SWINGERsim, restart the worker; `worker.json` shows the version
it runs. Job files and the rules for them are described in [Workflows](#wf-worker).

## Building this documentation {#install-docs}

The documentation is written in Markdown in `docs/en/` (English) and `docs/cs/` (Czech), one
folder per document. The PDF files are built with pandoc and XeLaTeX:

```bash
sudo apt install pandoc texlive-xetex texlive-fonts-recommended texlive-lang-czechslovak
python tools/build_docs_pdf.py guide --lang en      # this guide: docs/pdf/SWINGERsim-guide-en.pdf
python tools/build_docs_pdf.py                      # every document in both languages
```
