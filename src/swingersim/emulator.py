"""Neural emulators of swingersim ensembles (F4 e): system parameters -> dynamical observables.

An ensemble (``swingersim ensemble``) maps the parameters of the members (masses, orbits, angles)
to summary columns: present rates of change of the inclinations and tilts, eclipse fractions,
eccentricity maxima, Roche-lobe filling, the stability verdict of the N-body ghost test, the
fate. Every member costs seconds to hours. An emulator learns this map from a finished
ensemble and evaluates it in microseconds, so that a prior can be sampled with 10^6-10^7
draws and weighted by the observations (``swingersim emulator-report``), or the stability of a
system classified at once.

The networks are multi-layer perceptrons (scikit-learn), as the stability classifiers of
Vynatheya et al. (2022, 2023): regression targets are standardised (optionally after log10)
and fitted jointly; every binary class (0/1 column, or ``column==value`` for text columns such
as ``fate==intact``) has its own classifier. A deep ensemble of ``n_models`` networks with
different seeds gives the spread of the predictions (epistemic uncertainty); the residual
scatter on the held-out test set (``test_fraction``) is stored per target in the space the
networks fit (log10 for log targets, divided by the scale of a scaled target) and added in
quadrature, converted to the units of the target at every prediction, when the emulator replaces
the simulations in a likelihood.
"""

from __future__ import annotations

import json
import math
import pickle
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

FORMAT = 1


def _class_spec(spec: str) -> tuple[str, str | None]:
    """'ghost_stable' -> ('ghost_stable', None); 'fate==intact' -> ('fate', 'intact')."""
    if "==" in spec:
        c, v = spec.split("==", 1)
        return c.strip(), v.strip()
    return spec, None


def dataset(
    rows: list[dict],
    inputs: list[str],
    targets: list[str] = (),
    classes: list[str] = (),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int]]:
    """X (members x inputs), Y (members x targets, NaN where missing), C (members x classes,
    NaN where missing) from ensemble rows; members with a non-finite input or an error are
    left out. Columns may carry forward:/backward:/central: prefixes (ensemble.member_value)."""
    from .ensemble import index_runs, member_value, split_column

    runs = index_runs(rows)
    members = sorted(set(runs["forward"]) | set(runs["backward"]))
    X, Y, C, keep = [], [], [], []
    for m in members:
        if any((runs[d].get(m) or {}).get("fate") == "error" for d in runs):
            continue
        x = [member_value(runs, m, c) for c in inputs]
        if not all(math.isfinite(v) for v in x):
            continue
        y = [member_value(runs, m, c) for c in targets]
        cl = []
        for spec in classes:
            col, val = _class_spec(spec)
            if val is None:
                cl.append(member_value(runs, m, col))
            else:
                d, c = split_column(col)
                r = runs[d].get(m) or runs["backward" if d == "forward" else "forward"].get(m)
                cl.append(math.nan if r is None else float(str(r.get(c)) == val))
        X.append(x)
        Y.append(y)
        C.append(cl)
        keep.append(m)
    return (
        np.array(X, float).reshape(len(X), len(inputs)),
        np.array(Y, float).reshape(len(X), len(targets)),
        np.array(C, float).reshape(len(X), len(classes)),
        keep,
    )


@dataclass
class Emulator:
    inputs: list[str]
    targets: list[str]
    classes: list[str]
    log_inputs: list[str]
    log_targets: list[str]
    x_mean: np.ndarray
    x_std: np.ndarray
    y_mean: np.ndarray
    y_std: np.ndarray
    regressors: list = field(default_factory=list, repr=False)
    classifiers: dict = field(default_factory=dict, repr=False)
    resid_std: np.ndarray | None = None  # held-out RMS residual of every target (its units)
    resid_net: np.ndarray | None = None  # the same in the space the networks fit (log10, scaled)
    metrics: dict = field(default_factory=dict)
    info: dict = field(default_factory=dict)
    test: dict = field(default_factory=dict, repr=False)  # held-out X, Y, C for diagnostics
    angle_inputs: list = field(default_factory=list)  # inputs in degrees -> (cos, sin)
    target_scales: dict = field(default_factory=dict)  # target -> [input, power, reference]

    def _scale(self, X: np.ndarray) -> np.ndarray:
        """Factor (x / reference)^power of every target (1 without a scale), n x targets."""
        f = np.ones((len(X), len(self.targets)))
        for c, (inp, power, ref) in getattr(self, "target_scales", {}).items():
            f[:, self.targets.index(c)] = (X[:, self.inputs.index(inp)] / ref) ** power
        return f

    def net_space(self, X, Y) -> np.ndarray:
        """Targets in the space the networks fit: divided by their scale, log10 for log
        targets (before standardisation)."""
        Y = np.asarray(Y, float) / self._scale(np.atleast_2d(np.asarray(X, float)))
        for k, c in enumerate(self.targets):
            if c in self.log_targets:
                Y[:, k] = np.where(Y[:, k] > 0, np.log10(np.where(Y[:, k] > 0, Y[:, k], 1)), np.nan)
        return Y

    # ------------------------------------------------------------------ transforms
    def _x(self, X: np.ndarray) -> np.ndarray:
        X = np.array(X, float, copy=True)
        for k, c in enumerate(self.inputs):
            if c in self.log_inputs:
                X[:, k] = np.log10(X[:, k])
        Z = (X - self.x_mean) / self.x_std
        angles = [k for k, c in enumerate(self.inputs) if c in getattr(self, "angle_inputs", [])]
        if angles:  # a periodic input becomes a point on the unit circle
            rad = np.radians(X[:, angles])
            Z[:, angles] = np.cos(rad)
            Z = np.column_stack([Z, np.sin(rad)])
        return Z

    def _y_from(self, Z: np.ndarray) -> np.ndarray:
        Y = Z * self.y_std + self.y_mean
        for k, c in enumerate(self.targets):
            if c in self.log_targets:
                Y[:, k] = 10.0 ** Y[:, k]
        return Y

    # ------------------------------------------------------------------ prediction
    def predict(self, X) -> dict:
        """{'mean': (n x targets), 'std': (n x targets) ensemble spread, 'prob': {class:
        probability of 1 (mean over the ensemble)}} for inputs X (n x inputs, original units)."""
        X = np.atleast_2d(np.asarray(X, float))
        Z = self._x(X)
        out = {}
        if self.targets:
            f = self._scale(X)
            preds = np.array(
                [
                    self._y_from(np.atleast_2d(r.predict(Z)).reshape(len(Z), -1)) * f
                    for r in self.regressors
                ]
            )
            out["mean"] = preds.mean(axis=0)
            out["std"] = preds.std(axis=0)
        out["prob"] = {
            c: np.mean([m.predict_proba(Z)[:, 1] for m in models], axis=0)
            for c, models in self.classifiers.items()
        }
        return out

    def save(self, path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump({"format": FORMAT, "emulator": self}, f)
        path.with_suffix(".json").write_text(
            json.dumps(
                {
                    "inputs": self.inputs,
                    "targets": self.targets,
                    "classes": self.classes,
                    "log_inputs": self.log_inputs,
                    "log_targets": self.log_targets,
                    "angle_inputs": getattr(self, "angle_inputs", []),
                    "target_scales": getattr(self, "target_scales", {}),
                    "metrics": self.metrics,
                    "info": self.info,
                },
                indent=2,
                default=float,
            )
        )
        return path


def load(path) -> Emulator:
    with open(path, "rb") as f:
        d = pickle.load(f)
    if d.get("format") != FORMAT:
        raise ValueError(f"{path}: emulator format {d.get('format')} (expected {FORMAT})")
    return d["emulator"]


def _mlp_kwargs(hidden, seed, max_iter, activation, alpha, lr, batch):
    return {
        "hidden_layer_sizes": tuple(hidden),
        "activation": activation,
        "solver": "adam",
        "alpha": alpha,
        "batch_size": batch,
        "learning_rate_init": lr,
        "max_iter": max_iter,
        "early_stopping": True,
        "validation_fraction": 0.1,
        "n_iter_no_change": 20,
        "random_state": seed,
    }


def _fit_ensemble(make, Z, y, n_models: int, seed: int, tries: int = 3):
    """Fit n_models networks with seeds seed, seed + 1, ...; a network whose best validation
    score stays far below the others (stuck on a plateau and stopped early) is fitted again
    with a new seed, at most `tries` times. Returns the networks and the number of refits."""
    nets = []
    for s in range(n_models):
        m = make(seed + s)
        m.fit(Z, y)
        nets.append(m)
    refits = 0
    for _ in range(tries):
        scores = np.array([m.best_validation_score_ for m in nets])
        bad = np.nonzero(scores < np.median(scores) - 0.25)[0]
        if not len(bad):
            break
        for k in bad:
            refits += 1
            m = make(seed + 1000 * refits + int(k))
            m.fit(Z, y)
            nets[k] = m
    return nets, refits


def train(
    X: np.ndarray,
    Y: np.ndarray,
    C: np.ndarray,
    inputs: list[str],
    targets: list[str] = (),
    classes: list[str] = (),
    log_inputs: list[str] = (),
    log_targets: list[str] = (),
    angle_inputs: list[str] = (),
    target_scales: dict | None = None,
    hidden: tuple[int, ...] = (50, 50, 50, 50),
    n_models: int = 5,
    test_fraction: float = 0.2,
    seed: int = 0,
    max_iter: int = 2000,
    activation: str = "relu",
    alpha: float = 1e-4,
    learning_rate: float = 1e-3,
    batch_size: int | str = "auto",
) -> Emulator:
    """Train regressors (all targets jointly) and one classifier per class on a random split;
    the metrics are computed on the held-out test set. Inputs listed in ``angle_inputs`` (in
    degrees) enter the networks as (cos, sin), so that 359 deg and 1 deg are neighbours.
    ``target_scales`` {target: (input, power)} makes the networks fit target / input^power
    (relative to the median of the input), e.g. a precession rate times P_out^2: a known
    scaling taken out of the target makes it far easier to learn."""
    import warnings

    from sklearn.exceptions import ConvergenceWarning
    from sklearn.neural_network import MLPClassifier, MLPRegressor

    X = np.asarray(X, float)
    n = len(X)
    if n < 20:
        raise ValueError(f"only {n} usable members: too few to train an emulator")
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    n_test = max(1, round(test_fraction * n))
    test, tr = perm[:n_test], perm[n_test:]
    Xl = np.array(X, copy=True)
    for k, c in enumerate(inputs):
        if c in log_inputs:
            if np.any(Xl[:, k] <= 0):
                raise ValueError(f"log input {c} has values <= 0")
            Xl[:, k] = np.log10(Xl[:, k])
    x_mean, x_std = Xl[tr].mean(axis=0), Xl[tr].std(axis=0)
    x_std = np.where(x_std > 0, x_std, 1.0)
    em = Emulator(
        list(inputs),
        list(targets),
        list(classes),
        list(log_inputs),
        list(log_targets),
        x_mean,
        x_std,
        np.zeros(len(targets)),
        np.ones(len(targets)),
    )
    bad = [c for c in angle_inputs if c not in inputs or c in log_inputs]
    if bad:
        raise ValueError(f"angle inputs {bad} must be inputs and not log inputs")
    em.angle_inputs = list(angle_inputs)
    for c, (inp, power) in (target_scales or {}).items():
        if c not in targets or inp not in inputs:
            raise ValueError(f"target scale {c}: {inp}^{power} needs a target and an input")
        ref = float(np.median(X[:, list(inputs).index(inp)]))
        em.target_scales[c] = [inp, float(power), ref if ref != 0 else 1.0]
    Z = em._x(X)
    metrics: dict = {"n_train": len(tr), "n_test": len(test)}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        if targets:
            Y = np.asarray(Y, float)
            Yl = em.net_space(X, Y)
            ok = np.all(np.isfinite(Yl), axis=1)
            tr_y, te_y = tr[ok[tr]], test[ok[test]]
            if len(tr_y) < 20:
                raise ValueError("too few members with finite values of all targets")
            em.y_mean, em.y_std = Yl[tr_y].mean(axis=0), Yl[tr_y].std(axis=0)
            em.y_std = np.where(em.y_std > 0, em.y_std, 1.0)
            Zy = (Yl - em.y_mean) / em.y_std
            yy = Zy[tr_y] if len(targets) > 1 else Zy[tr_y, 0]
            em.regressors, retried = _fit_ensemble(
                lambda sd: MLPRegressor(
                    **_mlp_kwargs(
                        hidden, sd, max_iter, activation, alpha, learning_rate, batch_size
                    )
                ),
                Z[tr_y],
                yy,
                n_models,
                seed,
            )
            metrics["regression:retrained_networks"] = retried
            pred = em.predict(X[te_y])["mean"]
            truth = Y[te_y]
            res = pred - truth
            em.resid_std = np.sqrt(np.mean(res**2, axis=0))
            # the residual scatter where the networks work (log10 / scaled): a log target is
            # predicted with a relative error, not with one absolute error for all its values
            z_te = em._x(X[te_y])
            net_mean = np.mean(
                [
                    np.atleast_2d(r.predict(z_te)).reshape(len(z_te), -1) * em.y_std + em.y_mean
                    for r in em.regressors
                ],
                axis=0,
            )
            em.resid_net = np.sqrt(np.nanmean((net_mean - Yl[te_y]) ** 2, axis=0))
            is_log = np.array([c in log_targets for c in targets])
            t_net = em.net_space(X[te_y], truth)
            p_net = em.net_space(X[te_y], np.where(is_log & (pred <= 0), 1e-300, pred))
            for k, c in enumerate(targets):
                # R² in the space the network fits (scaled, log10 for log targets)
                t_k, p_k = t_net[:, k], p_net[:, k]
                var = np.var(t_k)
                metrics[f"{c}:rmse"] = float(em.resid_std[k])
                metrics[f"{c}:r2"] = (
                    float(1.0 - np.mean((p_k - t_k) ** 2) / var) if var > 0 else math.nan
                )
        for k, spec in enumerate(classes):
            lab = np.asarray(C, float)[:, k]
            ok = np.isfinite(lab)
            tr_c, te_c = tr[ok[tr]], test[ok[test]]
            if len(np.unique(lab[tr_c])) < 2:
                raise ValueError(f"class {spec}: the training set has only one class")
            models, retried = _fit_ensemble(
                lambda sd: MLPClassifier(
                    **_mlp_kwargs(
                        hidden, sd, max_iter, activation, alpha, learning_rate, batch_size
                    )
                ),
                Z[tr_c],
                lab[tr_c].astype(int),
                n_models,
                seed,
            )
            metrics[f"{spec}:retrained_networks"] = retried
            em.classifiers[spec] = models
            p = np.mean([m.predict_proba(Z[te_c])[:, 1] for m in models], axis=0)
            metrics.update(
                {f"{spec}:{k2}": v for k2, v in classification_scores(lab[te_c], p >= 0.5).items()}
            )
    em.metrics = metrics
    em.test = {
        "X": X[test],
        "Y": np.asarray(Y, float)[test] if len(targets) else np.zeros((len(test), 0)),
        "C": np.asarray(C, float)[test] if len(classes) else np.zeros((len(test), 0)),
    }
    em.info = {
        "hidden": list(hidden),
        "n_models": n_models,
        "activation": activation,
        "alpha": alpha,
        "learning_rate": learning_rate,
        "test_fraction": test_fraction,
        "seed": seed,
    }
    return em


def classification_scores(truth, predicted) -> dict:
    """Overall score, precisions and recalls of a binary classification (1 = 'positive'),
    as defined by Vynatheya et al. (2023), Sect. 5.1 (there: 1 = stable)."""
    t = np.asarray(truth).astype(bool)
    p = np.asarray(predicted).astype(bool)
    tp, tn = int(np.sum(t & p)), int(np.sum(~t & ~p))
    fp, fn = int(np.sum(~t & p)), int(np.sum(t & ~p))

    def div(a, b):
        return a / b if b else math.nan

    return {
        "score": div(tp + tn, len(t)),
        "precision_1": div(tp, tp + fp),
        "precision_0": div(tn, tn + fn),
        "recall_1": div(tp, tp + fn),
        "recall_0": div(tn, tn + fp),
        "n": len(t),
    }


def diagnostics(em: Emulator) -> dict:
    """Held-out diagnostics of a trained emulator: for every target the truth, the prediction,
    the deep-ensemble spread and the z-scores (prediction - truth) / spread; for every class the
    truth, the predicted probability, the ROC curve and its area, and the reliability curve."""
    t = em.test
    if not t:
        raise ValueError(
            "the emulator has no stored test set (trained with an early development version)"
        )
    pr = em.predict(t["X"])
    out: dict = {"targets": {}, "classes": {}}
    if em.targets:
        # everything in the space the networks fit (scaled, log10 for log targets)
        Z = em._x(t["X"])
        nets = np.array(
            [
                np.atleast_2d(r.predict(Z)).reshape(len(Z), -1) * em.y_std + em.y_mean
                for r in em.regressors
            ]
        )
        truth = em.net_space(t["X"], t["Y"])
    for k, c in enumerate(em.targets):
        y = truth[:, k]
        ok = np.isfinite(y)
        mu, sd = nets[:, ok, k].mean(axis=0), nets[:, ok, k].std(axis=0)
        resid = np.sqrt(np.mean((mu - y[ok]) ** 2))
        out["targets"][c] = {
            "truth": y[ok],
            "pred": mu,
            "spread": sd,
            "z": (mu - y[ok]) / np.where(sd > 0, sd, np.nan),
            "z_total": (mu - y[ok]) / np.sqrt(sd**2 + resid**2),
            "label": _net_label(em, c),
        }
    for k, c in enumerate(em.classes):
        lab = t["C"][:, k]
        ok = np.isfinite(lab)
        y, p = lab[ok].astype(bool), pr["prob"][c][ok]
        thr = np.r_[np.inf, np.sort(np.unique(p))[::-1], -np.inf]
        tpr = np.array([np.mean(p[y] >= h) if y.any() else np.nan for h in thr])
        fpr = np.array([np.mean(p[~y] >= h) if (~y).any() else np.nan for h in thr])
        auc = (
            float((getattr(np, "trapezoid", None) or np.trapz)(tpr, fpr))
            if y.any() and (~y).any()
            else math.nan
        )
        bins = np.linspace(0, 1, 11)
        idx = np.clip(np.digitize(p, bins) - 1, 0, 9)
        rel = [
            (float(p[idx == b].mean()), float(y[idx == b].mean()), int((idx == b).sum()))
            for b in range(10)
            if (idx == b).any()
        ]
        out["classes"][c] = {
            "truth": y,
            "prob": p,
            "fpr": fpr,
            "tpr": tpr,
            "auc": auc,
            "reliability": rel,
        }
    return out


def _net_label(em: Emulator, c: str) -> str:
    lab = c
    sc = getattr(em, "target_scales", {}).get(c)
    if sc:
        lab = f"{c} / ({sc[0]}/{sc[2]:.3g})^{sc[1]:g}"
    return f"log10 {lab}" if c in em.log_targets else lab


def diagnostic_plot(em: Emulator, path, max_targets: int = 4) -> Path:
    """Figure of the training (loss and validation curves of every network), the held-out
    predictions against the truth, the calibration of the deep-ensemble spread and, for
    classifiers, the ROC and reliability curves."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = diagnostics(em)
    targets = list(d["targets"])[:max_targets]
    classes = list(d["classes"])
    nets = [("regression", r) for r in em.regressors] + [
        (c, m) for c, ms in em.classifiers.items() for m in ms
    ]
    rows = 1 + (1 if targets else 0) + len(classes)
    ncol = max(len(targets), 3 if (classes or targets) else 2)
    fig, axes = plt.subplots(rows, ncol, figsize=(3.6 * ncol, 3.2 * rows), squeeze=False)
    for ax in axes.flat:
        ax.set_visible(False)
    # row 0: training curves
    ax = axes[0, 0]
    ax.set_visible(True)
    for name, m in nets:
        ax.plot(m.loss_curve_, lw=0.8, color="C0" if name == "regression" else "C1")
    ax.set_yscale("log")
    ax.set_xlabel("epoch")
    ax.set_ylabel("training loss")
    ax.set_title("training loss of every network", fontsize=9)
    ax = axes[0, 1]
    ax.set_visible(True)
    for name, m in nets:
        if getattr(m, "validation_scores_", None) is not None:
            ax.plot(m.validation_scores_, lw=0.8, color="C0" if name == "regression" else "C1")
    ax.set_xlabel("epoch")
    ax.set_ylabel("validation score (R² or accuracy)")
    ax.set_title("score on the 10 % validation split", fontsize=9)
    if targets:
        ax = axes[0, 2] if ncol > 2 else None
        if ax is not None:
            ax.set_visible(True)
            for k, c in enumerate(targets):
                z = d["targets"][c]["z"]
                z = z[np.isfinite(z)]
                ax.hist(np.clip(z, -6, 6), bins=40, histtype="step", density=True, label=c)
            zz = np.linspace(-6, 6, 200)
            ax.plot(zz, np.exp(-0.5 * zz**2) / math.sqrt(2 * math.pi), "k--", lw=1)
            ax.set_xlabel("(prediction - truth) / ensemble spread")
            ax.set_title("calibration of the spread (dashed: ideal)", fontsize=9)
            ax.legend(fontsize=6)
        for j, c in enumerate(targets):
            t = d["targets"][c]
            ax = axes[1, j]
            ax.set_visible(True)
            ax.errorbar(
                t["truth"],
                t["pred"],
                yerr=t["spread"],
                fmt=".",
                ms=2,
                lw=0.4,
                alpha=0.6,
                color="C0",
            )
            lo = float(np.nanmin([t["truth"].min(), t["pred"].min()]))
            hi = float(np.nanmax([t["truth"].max(), t["pred"].max()]))
            ax.plot([lo, hi], [lo, hi], "k-", lw=0.8)
            ax.set_xlabel(f"true {t['label']}", fontsize=7)
            ax.set_ylabel("emulated", fontsize=8)
            ax.set_title(f"R² = {em.metrics.get(c + ':r2', math.nan):.3f}", fontsize=9)
    r0 = 2 if targets else 1
    for i, c in enumerate(classes):
        t = d["classes"][c]
        ax = axes[r0 + i, 0]
        ax.set_visible(True)
        ax.plot(t["fpr"], t["tpr"], color="C1")
        ax.plot([0, 1], [0, 1], "k:", lw=0.8)
        ax.set_xlabel("false positive rate")
        ax.set_ylabel("true positive rate")
        ax.set_title(f"ROC of {c}: area {t['auc']:.3f}", fontsize=9)
        ax = axes[r0 + i, 1]
        ax.set_visible(True)
        rel = np.array(t["reliability"])
        if len(rel):
            ax.plot(rel[:, 0], rel[:, 1], "o-", color="C1")
        ax.plot([0, 1], [0, 1], "k:", lw=0.8)
        ax.set_xlabel("predicted probability")
        ax.set_ylabel("observed frequency")
        ax.set_title("reliability (dotted: ideal)", fontsize=9)
        if ncol > 2:
            ax = axes[r0 + i, 2]
            ax.set_visible(True)
            ax.hist(t["prob"][t["truth"]], bins=20, range=(0, 1), alpha=0.6, label=f"{c} = 1")
            ax.hist(t["prob"][~t["truth"]], bins=20, range=(0, 1), alpha=0.6, label=f"{c} = 0")
            ax.set_xlabel("predicted probability")
            ax.set_title("held-out members by true class", fontsize=9)
            ax.legend(fontsize=7)
    fig.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def train_from_ensemble(out_dir, inputs, targets=(), classes=(), **kw) -> Emulator:
    from .ensemble import load_rows

    if not (Path(out_dir) / "rows.jsonl").exists():
        raise ValueError(
            f"the ensemble {out_dir} has not been computed yet (no rows.jsonl): "
            "run the ensemble job first, then this training job"
        )
    rows = load_rows(out_dir)
    X, Y, C, _ = dataset(rows, list(inputs), list(targets), list(classes))
    em = train(X, Y, C, list(inputs), list(targets), list(classes), **kw)
    em.info["ensemble"] = str(out_dir)
    return em


# ------------------------------------------------------------------ inference with an emulator
def prior_sample(cfg, inputs: list[str], n: int, seed: int = 0) -> np.ndarray:
    """n draws of the emulator inputs from the ``ensemble.vary`` distributions of cfg (inputs
    must be varied parameters; fixed values are allowed). The draws are independent: an
    ensemble that redraws members (require, ma01_window, ma01_unstable_keep) has a different
    prior, which the emulator cannot reproduce - refused."""
    from .ensemble import draw_values

    ens = cfg.ensemble
    rules = [
        name
        for name, on in (
            ("require", bool(ens.require)),
            ("ma01_window", ens.ma01_window is not None),
            ("ma01_unstable_keep", ens.ma01_unstable_keep < 1.0),
        )
        if on
    ]
    if rules:
        raise ValueError(
            f"ensemble.{', ensemble.'.join(rules)} select the members after drawing them: the "
            "emulator's prior (independent draws of ensemble.vary) would differ from the "
            "ensemble's; weight the ensemble itself (swingersim ensemble-report) instead"
        )
    vary = ens.vary
    missing = [c for c in inputs if c not in vary]
    if missing:
        raise ValueError(
            f"emulator inputs {missing} are not ensemble.vary parameters of this configuration; "
            "emulator-report can only sample varied parameters"
        )
    rng = np.random.default_rng(seed)
    return np.column_stack([draw_values(vary[c], rng, n) for c in inputs]).reshape(n, len(inputs))


def _wanted_class(col: str, ob: dict) -> int:
    """Observed class (1 or 0) of a classifier column: value rounded, or a limit that only one
    class satisfies ({min: 0.5} -> 1, {max: 0.5} -> 0)."""
    if ob.get("value") is not None:
        return round(ob["value"])
    lo, hi = ob.get("min"), ob.get("max")
    if lo is not None and lo > 0.0 and (hi is None or hi >= 1.0):
        return 1
    if hi is not None and hi < 1.0 and (lo is None or lo <= 0.0):
        return 0
    raise ValueError(
        f"observed class '{col}': give value 0 or 1, or a limit that selects one class "
        "(min: 0.5 for 1, max: 0.5 for 0)"
    )


def _resid_var(em: Emulator, X: np.ndarray, k: int, mu: np.ndarray) -> np.ndarray | float:
    """Held-out residual variance of target k in its own units at every draw: the residual of
    the network space (log10 or scaled) times the derivative of the back-transformation (a
    relative error for a log target); emulators trained by older versions (quadsim) have only the
    residual in the target's units."""
    r_net = getattr(em, "resid_net", None)
    if r_net is None:
        return em.resid_std[k] ** 2 if em.resid_std is not None else 0.0
    c = em.targets[k]
    if c in em.log_targets:
        d = math.log(10.0) * np.abs(mu)
    else:
        d = np.abs(em._scale(X)[:, k])
    return (r_net[k] * d) ** 2


def emulated_weights(em: Emulator, X: np.ndarray, observed: dict) -> np.ndarray:
    """Likelihood of every draw: Gaussian in each observed target with sigma^2 + emulator
    variance (ensemble spread^2 + held-out residual^2); limits as probabilities; a class
    'col: {equals: ...}' or 'col==value' multiplies by the predicted probability."""
    from scipy.stats import norm

    pr = em.predict(X)
    w = np.ones(len(X))
    for col, ob in observed.items():
        ob = ob if isinstance(ob, dict) else ob.model_dump()
        if col in em.classifiers:  # a 0/1 class (e.g. ghost_stable): which class is observed?
            p = pr["prob"][col]
            w *= p if _wanted_class(col, ob) else 1.0 - p
            continue
        if ob.get("equals") is not None:
            vals = ob["equals"] if isinstance(ob["equals"], list) else [ob["equals"]]
            ps = [pr["prob"].get(f"{col}=={v}") for v in vals]
            if any(p is None for p in ps):
                raise ValueError(f"no classifier for {col}=={vals} in the emulator")
            w *= np.clip(np.sum(ps, axis=0), 0.0, 1.0)
            continue
        if col not in em.targets:
            raise ValueError(f"observed column {col} is not an emulator target or class")
        k = em.targets.index(col)
        mu = pr["mean"][:, k]
        s2 = pr["std"][:, k] ** 2 + _resid_var(em, X, k, mu)
        if ob.get("value") is not None:
            s = np.sqrt(ob["sigma"] ** 2 + s2)
            w *= np.exp(-0.5 * ((mu - ob["value"]) / s) ** 2)
        sd = np.sqrt(np.maximum(s2, 1e-300))
        if ob.get("min") is not None:
            w *= norm.sf((ob["min"] - mu) / sd)
        if ob.get("max") is not None:
            w *= norm.cdf((ob["max"] - mu) / sd)
    return w


def emulated_report(em: Emulator, cfg, n: int, observed: dict | None = None, seed: int = 0) -> str:
    """Posterior summary of the varied parameters from n prior draws weighted by the emulated
    likelihood of the observations (cfg.ensemble.observed unless given)."""
    from .ensemble import _wquantile, effective_size

    observed = observed if observed is not None else cfg.ensemble.observed
    if not observed:
        raise ValueError("no observations: give ensemble.observed in the configuration")
    X = prior_sample(cfg, em.inputs, n, seed)
    w = emulated_weights(em, X, observed)
    lines = [
        (
            f"emulated posterior: {n} prior draws, effective sample size {effective_size(w):.1f}, "
            f"mean likelihood {w.mean():.3g}"
        )
    ]
    if w.sum() <= 0:
        lines.append("  no draw is compatible with the observations")
        return "\n".join(lines)
    for k, c in enumerate(em.inputs):
        q = _wquantile(X[:, k], w, [0.16, 0.5, 0.84])
        q2 = _wquantile(X[:, k], w, [0.025, 0.975])
        lines.append(
            f"  {c:28s} median {q[1]:.5g}  68 %: {q[0]:.5g} .. {q[2]:.5g}  "
            f"95 %: {q2[0]:.5g} .. {q2[1]:.5g}"
        )
    return "\n".join(lines)
