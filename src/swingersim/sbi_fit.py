"""Simulation-based inference (SBI) on swingersim ensembles.

A Monte Carlo ensemble is a set of simulations (parameters theta drawn from the ensemble's
distributions, observables x computed by swingersim). Neural posterior estimation (NPE, Papamakarios
& Murray 2016; Greenberg et al. 2019) trains a conditional density estimator q(theta | x) on these
pairs with the `sbi` package (Tejero-Cantero et al. 2020). Afterwards the
posterior for a measured x_obs is sampled in a fraction of a second, without new simulations -
for the same system with new data, or for every system of a class covered by the ensemble.

Relation to the importance weighting of `swingersim ensemble-report --observed` (ensemble.py): both
give p(theta | x_obs) under the ensemble's distributions as the prior. Weighting needs no
training but its effective sample size collapses when the measurement is precise; NPE
interpolates between members and keeps working with few close members, at the cost of a
training step and of the usual neural-network caveats (check it against the weighting, and with
the coverage test of `coverage()`).

Observations use the format of `ensemble.observed`: {value, sigma} entries are observables
(the training set gets Gaussian noise of that sigma, so the posterior includes the measurement
error); {min, max} entries are hard limits that select the training members (a truncated
prior). Columns may be written `backward:<column>`; the present rates of change are central
(mean of forward and backward) unless prefixed, as in `ensemble.observation_weights`.

Needs the optional dependencies `sbi` and `torch` (`pip install -e ".[sbi]"`).

Data augmentation and validation: every member enters the training set `copies` times with
independent noise. `sbi` chooses its validation set (early stopping) at random among the
appended simulations and has no option for a split by member, so copies of one member can be on
both sides and the validation loss is somewhat optimistic. The coverage and calibration tests
(`swingersim sbi --coverage`) use members held out before the augmentation and are not affected.
"""

from __future__ import annotations

import hashlib
import json
import math
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .ensemble import member_value


def training_set(
    rows: list[dict], params: list[str], observed: dict, log_params: list[str] = ()
) -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray, np.ndarray]:
    """theta (members x params), x (members x observables) of the members that pass every hard
    limit and have finite values; also the observable names, their sigmas and observed values."""
    from .ensemble import index_runs, passes_limits

    runs = index_runs(rows)
    obs_cols = [c for c, o in observed.items() if _ob(o).get("value") is not None]
    if not obs_cols:
        raise ValueError("no observable with a value: give at least one {value, sigma}")
    limits = {c: _ob(o) for c, o in observed.items() if _ob(o).get("value") is None}
    theta, x = [], []
    not_positive = 0
    for m in sorted(set(runs["forward"]) | set(runs["backward"])):

        def get(col, m=m):
            return member_value(runs, m, col)

        if not passes_limits(runs, m, limits):
            continue
        th = [get(p) for p in params]
        xx = [get(c) for c in obs_cols]
        if not all(math.isfinite(v) for v in th + xx):
            continue
        if any(p in log_params and v <= 0 for p, v in zip(params, th)):
            not_positive += 1
            continue
        th = [math.log10(v) if p in log_params else v for p, v in zip(params, th)]
        theta.append(th)
        x.append(xx)
    if not_positive:
        warnings.warn(
            f"{not_positive} members left out: a parameter inferred in log10 is <= 0",
            stacklevel=2,
        )
    sig = np.array([_ob(observed[c])["sigma"] for c in obs_cols], float)
    val = np.array([_ob(observed[c])["value"] for c in obs_cols], float)
    return np.array(theta), np.array(x), obs_cols, sig, val


def _ob(o) -> dict:
    return o if isinstance(o, dict) else o.model_dump()


def training_fingerprint(directory, n_rows: int, observed: dict) -> str:
    """Identity of a training set: the ensemble (fingerprint of its ensemble.yaml), the number
    of its rows and the hard limits of the observations. A saved model is reused only when it
    matches (more members or other limits give another training set)."""
    import yaml

    meta = Path(directory) / "ensemble.yaml"
    try:
        ens_fp = (yaml.safe_load(meta.read_text()) or {}).get("fingerprint", "")
    except (OSError, yaml.YAMLError):
        ens_fp = ""
    limits = {c: _ob(o) for c, o in observed.items() if _ob(o).get("value") is None}
    data = {"ensemble": ens_fp, "rows": int(n_rows), "limits": limits}
    return hashlib.sha256(json.dumps(data, sort_keys=True, default=str).encode()).hexdigest()[:16]


@dataclass
class SBIModel:
    params: list[str]
    log_params: list[str]
    observables: list[str]
    sigmas: np.ndarray
    n_train: int
    low: np.ndarray
    high: np.ndarray
    posterior: object = field(repr=False)
    history: dict = field(default_factory=dict)  # training and validation loss per epoch
    fingerprint: str = ""  # training_fingerprint() of the training set

    def sample(self, x_obs, n: int = 10000, seed: int = 0) -> np.ndarray:
        import torch

        torch.manual_seed(seed)
        xo = torch.as_tensor(np.asarray(x_obs, float), dtype=torch.float32)
        s = self.posterior.sample((n,), x=xo, show_progress_bars=False)
        return s.numpy()

    def save(self, path) -> None:
        import torch

        torch.save(self, path)


def load_model(path) -> SBIModel:
    import torch

    return torch.load(path, weights_only=False)


def train(
    theta: np.ndarray,
    x: np.ndarray,
    sigmas: np.ndarray,
    params: list[str],
    observables: list[str],
    log_params: list[str] = (),
    seed: int = 0,
    copies: int = 5,
    max_epochs: int = 2000,
    fingerprint: str = "",
) -> SBIModel:
    """Train NPE on (theta, x + noise). Each member is used `copies` times with independent
    measurement noise of the given sigmas (data augmentation; see the module notes on the
    validation split)."""
    import torch
    from sbi.inference import NPE
    from sbi.utils import BoxUniform

    if len(theta) < 20:
        raise ValueError(f"only {len(theta)} usable members: too few to train")
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    th = np.repeat(theta, copies, axis=0)
    xx = np.repeat(x, copies, axis=0) + rng.normal(size=(len(theta) * copies, x.shape[1])) * sigmas
    span = theta.max(axis=0) - theta.min(axis=0)
    span = np.where(span > 0, span, 1.0)
    low, high = theta.min(axis=0) - 1e-3 * span, theta.max(axis=0) + 1e-3 * span
    prior = BoxUniform(
        low=torch.as_tensor(low, dtype=torch.float32),
        high=torch.as_tensor(high, dtype=torch.float32),
    )
    inference = NPE(prior=prior, show_progress_bars=False)
    inference.append_simulations(
        torch.as_tensor(th, dtype=torch.float32), torch.as_tensor(xx, dtype=torch.float32)
    )
    estimator = inference.train(max_num_epochs=max_epochs, show_train_summary=False)
    posterior = inference.build_posterior(estimator)
    summ = getattr(inference, "summary", None) or getattr(inference, "_summary", {}) or {}
    history = {
        k: [float(v) for v in summ.get(k, [])]
        for k in ("training_loss", "validation_loss")
        if summ.get(k)
    }
    return SBIModel(
        list(params),
        list(log_params),
        list(observables),
        np.asarray(sigmas),
        len(theta),
        low,
        high,
        posterior,
        history,
        fingerprint,
    )


def coverage(
    model: SBIModel,
    theta: np.ndarray,
    x: np.ndarray,
    n_test: int = 200,
    level: float = 0.68,
    n_samples: int = 1000,
    seed: int = 1,
) -> float:
    """Fraction of test members whose true parameters lie inside the central `level` interval
    of their posterior (all parameters, marginally, averaged). Ideally equal to `level`. The
    test members should not have been used for training."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(theta), size=min(n_test, len(theta)), replace=False)
    lo_q, hi_q = 50 * (1 - level), 50 * (1 + level)
    inside = []
    for k in idx:
        xo = x[k] + rng.normal(size=x.shape[1]) * model.sigmas
        s = model.sample(xo, n_samples, seed=int(k))
        lo, hi = np.percentile(s, [lo_q, hi_q], axis=0)
        inside.append(((theta[k] >= lo) & (theta[k] <= hi)).mean())
    return float(np.mean(inside))


def calibration(
    model: SBIModel,
    theta: np.ndarray,
    x: np.ndarray,
    n_test: int = 200,
    n_samples: int = 500,
    seed: int = 1,
) -> np.ndarray:
    """Simulation-based calibration (Talts et al. 2018): for held-out members, the fraction of
    posterior samples below the true value of every parameter (n_test x parameters). For a
    calibrated posterior these fractions are uniform on [0, 1]."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(theta), size=min(n_test, len(theta)), replace=False)
    u = []
    for k in idx:
        xo = x[k] + rng.normal(size=x.shape[1]) * model.sigmas
        smp = model.sample(xo, n_samples, seed=int(k))
        u.append((smp < theta[k]).mean(axis=0))
    return np.array(u)


def coverage_curve(u: np.ndarray, levels=None) -> tuple[np.ndarray, np.ndarray]:
    """Empirical coverage of central credible intervals at the given levels, from the
    calibration fractions u of calibration(): (levels, parameters)."""
    levels = np.linspace(0.05, 0.95, 19) if levels is None else np.asarray(levels, float)
    cov = np.array([(np.abs(u - 0.5) <= lv / 2).mean(axis=0) for lv in levels])
    return levels, cov


def plot_diagnostics(model: SBIModel, u: np.ndarray | None, path) -> Path:
    """Training and validation loss per epoch, the coverage curve and the calibration
    histograms of every parameter."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(model.params)
    ncol = 2 + (n if u is not None else 0)
    fig, axes = plt.subplots(1, ncol, figsize=(3.4 * ncol, 3.1), squeeze=False)
    axes = axes[0]
    h = model.history
    for k, lab in (("training_loss", "training"), ("validation_loss", "validation")):
        if h.get(k):
            axes[0].plot(h[k], label=lab)
    axes[0].set_xlabel("epoch")
    axes[0].set_ylabel("loss (- log q)")
    axes[0].set_title("training of the density estimator", fontsize=9)
    axes[0].legend(fontsize=8)
    if u is not None:
        lv, cov = coverage_curve(u)
        for j, p in enumerate(model.params):
            axes[1].plot(lv, cov[:, j], label=p)
        axes[1].plot([0, 1], [0, 1], "k:", lw=1)
        axes[1].set_xlabel("credible level")
        axes[1].set_ylabel("held-out members inside")
        axes[1].set_title("coverage (dotted: calibrated)", fontsize=9)
        axes[1].legend(fontsize=6)
        nb = 10
        lo, hi = _binom_band(len(u), 1 / nb)
        for j, p in enumerate(model.params):
            ax = axes[2 + j]
            ax.hist(u[:, j], bins=nb, range=(0, 1), color="#3b6ea5", alpha=0.8)
            ax.axhspan(lo, hi, color="grey", alpha=0.25, lw=0)
            ax.set_xlabel("rank of the true value", fontsize=8)
            ax.set_title(p, fontsize=9)
    else:
        axes[1].axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return Path(path)


def _binom_band(n: int, p: float) -> tuple[float, float]:
    """Central 99 % range of a bin count of a uniform histogram (normal approximation)."""
    m, s = n * p, math.sqrt(n * p * (1 - p))
    return max(m - 2.58 * s, 0.0), m + 2.58 * s


def summarize_posterior(model: SBIModel, samples: np.ndarray, x_obs) -> str:
    lines = [
        f"SBI posterior (NPE, {model.n_train} training members) for "
        + ", ".join(f"{o} = {v:g}" for o, v in zip(model.observables, x_obs)),
        f"  {'parameter':28s} {'median':>12s} {'16 %':>12s} {'84 %':>12s} {'2.5 %':>12s} {'97.5 %':>12s}",
    ]
    for k, p in enumerate(model.params):
        v = samples[:, k]
        if p in model.log_params:
            v = 10.0**v
        q = np.percentile(v, [50, 16, 84, 2.5, 97.5])
        lines.append(f"  {p:28s} " + " ".join(f"{x:12.5g}" for x in q))
    return "\n".join(lines)


def plot_posterior(model: SBIModel, samples: np.ndarray, path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(model.params)
    fig, axes = plt.subplots(n, n, figsize=(2.4 * n + 1, 2.4 * n + 1), squeeze=False)
    labels = [f"log10 {p}" if p in model.log_params else p for p in model.params]
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if j > i:
                ax.axis("off")
            elif i == j:
                ax.hist(samples[:, i], bins=40, color="#3b6ea5", histtype="stepfilled", alpha=0.8)
                ax.set_yticks([])
            else:
                ax.hist2d(samples[:, j], samples[:, i], bins=40, cmap="Blues")
            if i == n - 1:
                ax.set_xlabel(labels[j], fontsize=8)
            if j == 0 and i > 0:
                ax.set_ylabel(labels[i], fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return Path(path)
