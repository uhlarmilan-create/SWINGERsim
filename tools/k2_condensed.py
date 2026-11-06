"""Apsidal-motion constant k2 of a core-dominated (stripped) giant: a compact core of mass Mc and
radius rc under a polytropic envelope of mass Me reaching the radius R (Part IV, CzeV343).

The envelope is a polytrope of index n in the gravity of the whole star (valid for Me << Mc):
hydrostatic equilibrium gives rho(r) ~ (1/r - 1/R)^n. The Clairaut-Radau equation
    r d(eta)/dr = 6 - eta^2 + eta - 6 (rho / rho_mean)(eta + 1),   eta = d ln(epsilon) / d ln r,
is integrated outwards from eta = 0 inside the (uniform) core, and
    k2 = (3 - eta(R)) / (2 (2 + eta(R)))
(Kopal 1959; Claret & Gimenez 1992). Checks: a uniform sphere gives 0.75.

    python tools/k2_condensed.py [--mc 0.258] [--me 0.004] [--radius 14] [--rc 0.03] [--n 1.5]
"""

import argparse

import numpy as np
from scipy.integrate import quad, solve_ivp


def k2_condensed(mc: float, me: float, radius: float, rc: float, n: float = 1.5) -> float:
    def shape(r):
        return (1.0 / r - 1.0 / radius) ** n

    norm = me / quad(lambda r: 4 * np.pi * r * r * shape(r), rc, radius, limit=400)[0]

    def rho(r):
        return norm * shape(r)

    def mass(r):
        return mc + quad(lambda x: 4 * np.pi * x * x * rho(x), rc, r, limit=400)[0]

    def rhs(r, y):
        eta = y[0]
        mean = mass(r) / (4.0 / 3.0 * np.pi * r**3)
        return [(6.0 - eta * eta + eta - 6.0 * rho(r) / mean * (eta + 1.0)) / r]

    sol = solve_ivp(rhs, [rc, radius], [0.0], rtol=1e-8, atol=1e-10, max_step=radius / 2000)
    eta = float(sol.y[0, -1])
    return (3.0 - eta) / (2.0 * (2.0 + eta))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mc", type=float, default=0.258, help="core mass [Msun]")
    ap.add_argument("--me", type=float, nargs="+", default=[0.004, 0.01, 0.02, 0.05, 0.07, 0.1])
    ap.add_argument("--radius", type=float, default=14.0, help="stellar radius [Rsun]")
    ap.add_argument("--rc", type=float, default=0.03, help="core radius [Rsun]")
    ap.add_argument(
        "--n", type=float, nargs="+", default=[1.5, 3.0], help="envelope polytrope index"
    )
    a = ap.parse_args(argv)
    print("n     M_env [Msun]  M_env/M   k2")
    for n in a.n:
        for me in a.me:
            k2 = k2_condensed(a.mc, me, a.radius, a.rc, n)
            print(f"{n:<5g} {me:<13g} {me / (a.mc + me):<9.3f} {k2:.2e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
