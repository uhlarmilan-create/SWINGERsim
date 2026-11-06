#!/usr/bin/env bash
# SWINGERsim installer for Linux, macOS and Windows (inside WSL).
#
#   ./install.sh                 virtual environment .venv + SWINGERsim with all common extras
#   ./install.sh --minimal       only the core (runs, physics, HDF5 output)
#   ./install.sh --with-sbi      also simulation-based inference (PyTorch, about 1 GB)
#   ./install.sh --with-phoebe   also PHOEBE 2 as a light-curve engine
#   ./install.sh --venv DIR      put the virtual environment into DIR (default: .venv)
#   ./install.sh --python PY     use this Python interpreter (default: the newest python3.x >= 3.11)
#   ./install.sh --no-check      skip the check after the installation
#
# The script never uses sudo: missing system packages are reported with the command to install them.
set -euo pipefail
cd "$(dirname "$0")"

VENV=.venv
PY=""
EXTRAS="all,dev"
CHECK=1
while [ $# -gt 0 ]; do
  case "$1" in
    --minimal) EXTRAS="dev" ;;
    --with-sbi) EXTRAS="$EXTRAS,sbi" ;;
    --with-phoebe) EXTRAS="$EXTRAS,phoebe" ;;
    --venv) VENV="$2"; shift ;;
    --python) PY="$2"; shift ;;
    --no-check) CHECK=0 ;;
    -h|--help) sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1 (see ./install.sh --help)"; exit 2 ;;
  esac
  shift
done

say() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
fail() { printf '\n\033[31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }

# --- 1. Python >= 3.11 ---------------------------------------------------------------------
say "Looking for Python 3.11 or newer"
if [ -z "$PY" ]; then
  for c in python3.13 python3.12 python3.11 python3; do
    if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import sys; sys.exit(sys.version_info < (3, 11))'; then
      PY="$c"; break
    fi
  done
fi
[ -n "$PY" ] || fail "Python 3.11 or newer was not found.
  Debian/Ubuntu: sudo apt install python3 python3-venv python3-dev
  Fedora:        sudo dnf install python3 python3-devel
  macOS:         brew install python@3.12"
"$PY" -c 'import sys; sys.exit(sys.version_info < (3, 11))' || fail "$PY is older than 3.11"
echo "using $PY ($("$PY" --version 2>&1))"
"$PY" -c 'import venv, ensurepip' 2>/dev/null || fail "the Python module venv/ensurepip is missing.
  Debian/Ubuntu: sudo apt install python3-venv   (or python3.X-venv for your version)"

# --- 2. C compiler (REBOUND/REBOUNDx are built from source where no wheel exists) ------------
if ! command -v cc >/dev/null 2>&1 && ! command -v gcc >/dev/null 2>&1 && ! command -v clang >/dev/null 2>&1; then
  echo "warning: no C compiler found. If pip has to build REBOUND or REBOUNDx from source, install one:"
  echo "  Debian/Ubuntu: sudo apt install build-essential    macOS: xcode-select --install"
fi

# --- 3. Package data (stellar tracks and tables) ---------------------------------------------
say "Checking the data of the package"
DATA=src/swingersim/data
for f in mist sse_tracks.npz claret2023_structure.npz claret2004_structure.npz claret2011_ld.npz passbands; do
  [ -e "$DATA/$f" ] || fail "$DATA/$f is missing. Clone the whole repository (git clone, not a partial download)."
done
echo "stellar tracks and tables: $(du -sh "$DATA" | cut -f1)"

# --- 4. Virtual environment and installation -------------------------------------------------
say "Creating the virtual environment $VENV"
[ -x "$VENV/bin/python" ] || "$PY" -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"
python -m pip install --upgrade pip setuptools wheel >/dev/null

say "Installing SWINGERsim with the extras [$EXTRAS] (a few minutes)"
python -m pip install -e ".[$EXTRAS]"

# Where no REBOUNDx wheel exists (e.g. Python 3.14), pip builds it in a temporary environment and
# the library keeps a path to that environment's REBOUND, which no longer exists. Rebuild it
# against the REBOUND installed here.
if ! python -c "import reboundx" 2>/dev/null; then
  say "Rebuilding REBOUNDx against the installed REBOUND"
  python -m pip install --no-build-isolation --no-cache-dir --force-reinstall --no-deps \
    "reboundx==$(python -c 'import importlib.metadata as m; print(m.version("reboundx"))')"
  python -c "import reboundx" || fail "REBOUNDx does not load (see the message above)"
fi

# --- 5. Check --------------------------------------------------------------------------------
say "Checking the installation"
swingersim --help >/dev/null || fail "the command swingersim does not start"
python - <<'EOF'
import importlib, swingersim
print(f"SWINGERsim {swingersim.__version__} from {swingersim.__file__.rsplit('/', 2)[0]}")
for m in ("rebound", "reboundx", "numpy", "scipy", "h5py", "plotly", "numba", "emcee", "matplotlib"):
    try:
        mod = importlib.import_module(m)
        print(f"  {m:11s} {getattr(mod, '__version__', 'ok')}")
    except ImportError:
        print(f"  {m:11s} not installed (optional)")
EOF
if [ "$CHECK" = 1 ]; then
  mkdir -p output/check
  swingersim run configs/examples/example_2p2.yaml -o output/check/example.h5 >/dev/null \
    && echo "test run: output/check/example.h5" || fail "the test run failed"
  echo "fast benchmarks against published results (about a minute):"
  swingersim validate | tail -n 3
fi

say "Done"
cat <<EOF
Activate the environment in every new terminal:
    source $VENV/bin/activate
Then, for example:
    swingersim gui                                   # web application (opens in the browser)
    swingersim run configs/examples/example_2p2.yaml -o output/example.h5
    swingersim validate --all                        # all benchmarks (about an hour, all cores)
User Guide: docs/pdf/SWINGERsim-guide-en.pdf (Czech: SWINGERsim-guide-cs.pdf)
EOF
