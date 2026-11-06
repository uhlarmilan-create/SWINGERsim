# Data and reproducibility

The configurations are in `configs/squadra_test/` (published, anchored, fit and ensemble files for each
system), the measured minima, radial velocities and TESS light curves in `data/observed/<system>/`, and the
outputs of the fits in `output/squadra_test/<system>/`. Figures and numbers are made by
`python tools/squadra_report.py` and `python tools/squadra_fits.py`; the oscillation of V994 Her is removed by
`tools/remove_oscillation.py`, the SPOC crowding values are read by `tools/tess_crowding.py`. Parts I, III and IV of
the test series are separate documents (*Testing SWINGERsim*, Parts I, III, IV).
