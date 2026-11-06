# Data and reproducibility

The configurations are in `configs/squadra_test3/`, the observed data in `data/observed/<system>/`
(minima with their sources, inclinations and limits, radial velocities, TESS light curves), and the
outputs of the fits and ensembles in `output/squadra_test3/<system>/`. TESS minima were measured with
`tools/tess_minima.py`, the QX Cas depths with `tools/qx_cas_extra.py`, the Gaia data are queried by
`tools/gaia_check.py`. The figures and all
numbers quoted are produced by `python tools/squadra3_report.py` and stored in
`docs/figures/squadra_test3/numbers.json`.
