"""Explanatory figures of the Physics and Methods volume (docs/figures/methods/)."""

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = Path(__file__).resolve().parents[1] / "docs" / "figures" / "methods"
OUT.mkdir(parents=True, exist_ok=True)
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#dcdbd6"
plt.rcParams.update(
    {
        "font.size": 10,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "lines.linewidth": 2,
    }
)
G = 2.959122082855911e-4  # au^3 / (Msun d^2)
RSUN_AU = 0.00465047


# ------------------------------------------------------------------ 1. eclipses
def eclipses():
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.0))
    P = np.logspace(np.log10(0.5), np.log10(1000), 400)
    for m, r, lab, c in (
        (2.0, 2 * 0.89, "1 + 1 M☉, R = 0.89 R☉ each (ZAMS)", S1),
        (4.0, 2 * 1.70, "2 + 2 M☉, R = 1.70 R☉ each", S2),
        (2.0, 2 * 5.0, "1 + 1 M☉, R = 5 R☉ each (subgiants)", S3),
    ):
        acube = G * m * P**2 / (4 * math.pi**2)
        aau = acube ** (1 / 3)
        x = r * RSUN_AU / aau
        x = np.where(x < 0.7, x, np.nan)  # contact or overlapping stars: not shown
        a.plot(P, np.degrees(np.arcsin(x)), color=c, label=lab)
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("orbital period P [d]")
    a.set_ylabel("eclipse window |90° − i|max [deg]")
    a.set_title("(a) geometric eclipse limit, |cos i| < (R1 + R2)/a", fontsize=10)
    a.legend(fontsize=8, frameon=False)
    d = np.logspace(-4, np.log10(0.9), 400)
    for l3, c in ((0.0, S1), (0.5, S2), (0.9, S3), (0.99, S4)):
        dm = -2.5 * np.log10(1 - (1 - l3) * d)
        b.plot(d, dm, color=c, label=f"l3 = {l3:g}")
    b.axhline(1e-3, color=INK2, lw=1, ls="--")
    b.text(1.2e-4, 1.25e-3, "example detection limit 1 mmag", color=INK2, fontsize=8)
    b.set_xscale("log")
    b.set_yscale("log")
    b.set_xlabel("intrinsic fractional depth d of the binary")
    b.set_ylabel("observed depth Δm [mag]")
    b.set_title("(b) dilution by third light", fontsize=10)
    b.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "eclipse_detectability.png", dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------ 2. fit chain
def box(ax, x, y, w, h, title, body, fc):
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.03", fc=fc, ec=INK2, lw=1
        )
    )
    ax.text(x + w / 2, y + h - 0.03, title, ha="center", va="top", fontsize=9.5, weight="bold")
    ax.text(
        x + w / 2,
        y + h - 0.085,
        body,
        ha="center",
        va="top",
        fontsize=8,
        color=INK2,
        linespacing=1.35,
    )


def arrow(ax, p, q, text=None, rad=0.0, dy=0.015):
    ax.add_patch(
        FancyArrowPatch(
            p,
            q,
            arrowstyle="-|>",
            mutation_scale=12,
            color=INK,
            lw=1.2,
            connectionstyle=f"arc3,rad={rad}",
        )
    )
    if text:
        ax.text(
            (p[0] + q[0]) / 2,
            (p[1] + q[1]) / 2 + dy,
            text,
            fontsize=7.5,
            color=INK2,
            ha="center",
            va="bottom" if dy > 0 else "top",
        )


def fit_chain():
    fig, ax = plt.subplots(figsize=(10, 6.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.04, 1)
    ax.axis("off")
    light, mid = "#eef4fc", "#fdf0ea"
    w, h = 0.25, 0.25
    X = (0.02, 0.375, 0.73)
    Y = (0.71, 0.38, 0.05)
    box(
        ax,
        X[0],
        Y[0],
        w,
        h,
        "Data",
        "times of minima (per source)\nradial velocities\n"
        "light-curve segments\ninclinations / eclipse limits",
        "#f3f3f1",
    )
    box(
        ax,
        X[1],
        Y[0],
        w,
        h,
        "Stage 1 (analytic, per pair)",
        "T0 + P(n + j/2) + Q n²\n+ LTTE (Irwin 1952)\n+ apsidal term (Giménez &\nBastero 1995)",
        light,
    )
    box(
        ax,
        X[2],
        Y[0],
        w,
        h,
        "Starting configuration",
        "published elements, or\nstage-1 outer "
        "orbit; masses,\nradii, k2 from the stellar\nmodel or the literature",
        "#f3f3f1",
    )
    box(
        ax,
        X[2],
        Y[1],
        w,
        h,
        "Anchoring (optional)",
        "P and M0 of every eclipsing\npair set so "
        "that the model\nephemeris matches the data\n(removes osculating drift)",
        light,
    )
    box(
        ax,
        X[1],
        Y[1],
        w,
        h,
        "Stage 2: N-body fit",
        "one integration per direction\n"
        "residuals of all data + priors\ntrust-region least squares,\ncentral-difference Jacobian",
        mid,
    )
    box(
        ax,
        X[0],
        Y[1],
        w,
        h,
        "Fits in a row",
        "per-source factor or jitter\nσ-clipping "
        "(robust σ)\nred-noise factor β\nrepeat until nothing changes",
        mid,
    )
    box(
        ax,
        X[1],
        Y[2],
        w,
        h,
        "Uncertainties",
        "covariance from the SVD of the\nscaled Jacobian, "
        "× √max(χ²red, 1)\noptional emcee: walkers, τ;\nchain ≥ 50 τ recommended",
        light,
    )
    box(
        ax,
        X[2],
        Y[2],
        w,
        h,
        "Products",
        "fitted configuration (runs,\nstability, light "
        "curves)\npredicted minima with\n16–84 % bands",
        "#f3f3f1",
    )
    ym0, ym1, ym2 = Y[0] + h / 2, Y[1] + h / 2, Y[2] + h / 2
    arrow(ax, (X[0] + w, ym0), (X[1], ym0))
    arrow(ax, (X[1] + w, ym0), (X[2], ym0))
    arrow(ax, (X[2] + w / 2, Y[0]), (X[2] + w / 2, Y[1] + h))
    arrow(ax, (X[2], ym1), (X[1] + w, ym1))
    arrow(ax, (X[1], ym1 + 0.04), (X[0] + w, ym1 + 0.04), "residuals")
    arrow(ax, (X[0] + w, ym1 - 0.04), (X[1], ym1 - 0.04), "new σ, kept points", dy=-0.015)
    arrow(ax, (X[1] + w / 2, Y[1]), (X[1] + w / 2, Y[2] + h))
    arrow(ax, (X[1] + w, ym2), (X[2], ym2))
    ax.text(
        0.02,
        -0.02,
        "Blue boxes: analytic or linear steps; orange: steps that call the N-body "
        "model many times.",
        fontsize=8,
        color=INK2,
    )
    fig.savefig(OUT / "fit_chain.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ 3. weights
def weights():
    rng = np.random.default_rng(3)
    n = 4000
    theta = 10 ** rng.uniform(np.log10(40), np.log10(3000), n)  # outer period [d]
    # toy observable: a precession rate ~ P_out^-2 with 30 % intrinsic spread (angles)
    x = -1.4 * (120.0 / theta) ** 2 * rng.lognormal(0.0, 0.3, n)
    xo, so = -1.4, 0.2
    w = np.exp(-0.5 * ((x - xo) / so) ** 2)
    neff = w.sum() ** 2 / (w**2).sum()
    Z = w.mean()
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.0))
    o = np.argsort(w)
    sc = a.scatter(theta[o], x[o], c=w[o], s=9, cmap="Blues", vmin=0, vmax=1, edgecolors="none")
    a.axhspan(xo - so, xo + so, color=S2, alpha=0.18, lw=0)
    a.axhline(xo, color=S2, lw=1.2)
    a.set_xscale("log")
    a.set_ylim(-6, 0.2)
    a.set_xlabel("drawn parameter θ (toy: outer period [d])")
    a.set_ylabel("simulated observable x (toy rate [deg/yr])")
    a.set_title("(a) members coloured by their weight w", fontsize=10)
    cb = fig.colorbar(sc, ax=a)
    cb.set_label("weight w = exp(−χ²/2)")
    bins = np.logspace(np.log10(40), np.log10(3000), 36)
    b.hist(theta, bins=bins, density=True, color="#b9b8b2", label="prior (all members)")
    b.hist(
        theta,
        bins=bins,
        weights=w,
        density=True,
        histtype="step",
        color=S1,
        lw=2,
        label="conditioned (weighted)",
    )
    b.set_xscale("log")
    b.set_xlabel("drawn parameter θ (toy: outer period [d])")
    b.set_ylabel("probability density [1/d]")
    b.set_title("(b) prior and conditioned distribution", fontsize=10)
    b.legend(fontsize=8, frameon=False)
    b.text(
        0.97,
        0.62,
        f"n = {n}\nn_eff = {neff:.0f}\nmean weight Z = {Z:.3f}",
        transform=b.transAxes,
        ha="right",
        fontsize=8.5,
        color=INK,
    )
    fig.tight_layout()
    fig.savefig(OUT / "weights_toy.png", dpi=150)
    plt.close(fig)
    print("toy n_eff", neff, "Z", Z)


eclipses()
fit_chain()
weights()
