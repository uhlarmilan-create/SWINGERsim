"""Inclinations against time in the joint fit (fit.inclination, incfit.py)."""

from pathlib import Path

import numpy as np
import yaml

from swingersim.config import Config

ROOT = Path(__file__).resolve().parents[1]


def test_inclination_history_of_v907_sco():
    """The published solution of V907 Sco (Zasche et al. 2023) gives the measured inclinations of
    1900-2022 (N-body backward and forward from 2021.5) without any fit."""
    import pandas as pd

    from swingersim.incfit import load, model_inclination

    d = yaml.safe_load((ROOT / "configs/squadra_test3/v907_sco_published.yaml").read_text())
    cfg = Config.model_validate(d)
    obs = pd.read_csv(ROOT / "data/observed/v907_sco/inclination_zasche2023.csv", comment="#")
    m = model_inclination(cfg, "inner", obs.mid_hjd.to_numpy(), 2459383.9)
    assert np.all((m >= 0) & (m <= 90))
    chi2 = float(np.sum(((obs.i_deg - m) / obs.sigma_deg) ** 2))
    assert chi2 < 150  # 14 epochs; the TESS / FRAM ones to ~0.3 deg
    assert abs(m[9] - 76.61) < 0.05  # t = 0 itself

    class Spec:
        file = str(ROOT / "data/observed/v907_sco/inclination_zasche2023.csv")
        pair = "inner"
        label = None
        min_sigma_deg = 0.5

    q = load(Spec, None)
    assert len(q) == 14 and q.sigma.min() == 0.5


def test_fit_inclination_config():
    d = yaml.safe_load((ROOT / "configs/squadra_test3/v907_sco_fit.yaml").read_text())
    cfg = Config.model_validate(d)
    assert cfg.fit.inclination[0].pair == "inner"
    for name in ("qx_cas_fit", "ah_cep_fit", "qx_cas_ensemble", "v907_sco_fit_lc"):
        Config.model_validate(
            yaml.safe_load((ROOT / f"configs/squadra_test3/{name}.yaml").read_text())
        )
