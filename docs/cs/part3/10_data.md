# Data a reprodukovatelnost

Konfigurace jsou v `configs/squadra_test3/`, pozorovaná data v `data/observed/<system>/` (minima s jejich
zdroji, sklony a meze, radiální rychlosti, světelné křivky TESS) a výstupy fitů a ansámblů
v `output/squadra_test3/<system>/`. Minima z TESS byla změřena nástrojem `tools/tess_minima.py`, hloubky
QX Cas nástrojem `tools/qx_cas_extra.py`, data Gaia dotazuje `tools/gaia_check.py`. Obrázky a všechna
uváděná čísla vytváří `python tools/squadra3_report.py` a jsou uložena
v `docs/figures/squadra_test3/numbers.json`.
