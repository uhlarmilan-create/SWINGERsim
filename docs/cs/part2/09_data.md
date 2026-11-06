# Data a reprodukovatelnost

Konfigurace jsou v `configs/squadra_test/` (publikované, ukotvené, fitovací a ansámblové soubory pro každou
soustavu), naměřená minima, radiální rychlosti a světelné křivky TESS v `data/observed/<system>/` a výstupy
fitů v `output/squadra_test/<system>/`. Obrázky a čísla vytvářejí `python tools/squadra_report.py`
a `python tools/squadra_fits.py`; oscilaci V994 Her odstraňuje `tools/remove_oscillation.py`, hodnoty korekce
SPOC na přeplněné pole čte `tools/tess_crowding.py`. Části I, III a IV série testů jsou samostatné dokumenty
(*Testování SWINGERsimu*, části I, III, IV).
