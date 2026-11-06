"""Fast kernel of the secular 2+2 equations (the conservative part of SecularTwoTwo.rhs).

The same scalar code runs as plain Python or, when numba is installed, compiled with
numba.njit (about two orders of magnitude faster); results are identical to rounding. The
equations are those of swingersim.secular.twotwo (Vokrouhlicky 2016, eqs. 12-25, plus the e-only
GR and stellar-distortion terms); the test suite compares both implementations.

State y (18): K_A, pi_A, K_B, pi_B, K, pi (vectorial elements of A, B and the outer orbit).
c = (c2A, c2B, c3A, c3B), prec[x] = (L, n, alpha, beta, gamma) for x = A, B, out.
Returns (H, dy/dt).
"""

from __future__ import annotations

import math
import os

import numpy as np

try:  # optional acceleration; SWINGERSIM_NO_NUMBA=1 forces the plain Python version
    if os.environ.get("SWINGERSIM_NO_NUMBA"):
        raise ImportError
    from numba import njit

    HAVE_NUMBA = True
except ImportError:  # pragma: no cover - depends on the installation
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f


@njit(cache=True)
def _dot(y, i, j):
    return y[i] * y[j] + y[i + 1] * y[j + 1] + y[i + 2] * y[j + 2]


_PK_22 = np.array([0, 1])  # 2+2: A (orbit 0) and B (orbit 1) inside the outer orbit (2)
_PL_22 = np.array([2, 2])


@njit(cache=True)
def secular_rhs(y, c, prec, octupole):
    """2+2 system (orbits A, B, out; c = (c2_A, c2_B, c3_A, c3_B)): secular_rhs_pairs with the
    pairs A-out and B-out (one implementation of the pair terms; bit-identical results)."""
    return secular_rhs_pairs(y, _PK_22, _PL_22, c[:2], c[2:], prec, octupole)


@njit(cache=True)
def secular_rhs_pairs(y, pk, pl, c2s, c3s, prec, octupole):
    """General nested hierarchy: the same pair terms for every pair of orbits (k inside l).

    y holds (K, pi) of every orbit, 6 numbers each; pair p couples orbit pk[p] (inside) with
    orbit pl[p] (outside) with the quadrupole coefficient c2s[p] and octupole c3s[p];
    prec[x] = (L, n, alpha, beta, gamma) for orbit x. Returns (H, dy/dt)."""
    n_orb = prec.shape[0]
    g = np.zeros(6 * n_orb)
    out = np.zeros(6 * n_orb)
    H = 0.0
    for q in range(pk.shape[0]):
        ik = 6 * pk[q]
        ip = ik + 3
        jk = 6 * pl[q]
        jp = jk + 3
        c2 = c2s[q]
        c3 = c3s[q]
        e2 = _dot(y, jp, jp)
        eta = math.sqrt(1.0 - e2)
        em3 = eta**-3
        em5 = eta**-5
        em7 = eta**-7
        em9 = eta**-9
        ei2 = _dot(y, ip, ip)
        a_ = _dot(y, ik, jk)
        b_ = _dot(y, ip, jk)
        q1 = 6.0 * ei2 - 1.0
        q2 = 3.0 * a_ * a_ - 15.0 * b_ * b_
        H += -c2 * (em3 * q1 + em5 * q2)
        fp = -c2 * (3.0 * em5 * q1 + 5.0 * em7 * q2)
        for k in range(3):
            Kk = y[jk + k]
            Kik = y[ik + k]
            pik = y[ip + k]
            g[ip + k] += -c2 * (12.0 * em3 * pik - 30.0 * em5 * b_ * Kk)
            g[ik + k] += -c2 * (6.0 * em5 * a_ * Kk)
            g[jk + k] += -c2 * em5 * (6.0 * a_ * Kik - 30.0 * b_ * pik)
            g[jp + k] += fp * y[jp + k]
        if octupole and c3 != 0.0:
            d = _dot(y, ip, jp)
            cc = _dot(y, jp, ik)
            w = 35.0 * b_ * b_ - 5.0 * a_ * a_
            r = 1.0 - 8.0 * ei2
            H += -c3 * (em5 * r * d + em7 * d * w - 10.0 * em7 * a_ * b_ * cc)
            for k in range(3):
                Kk = y[jk + k]
                Kik = y[ik + k]
                pik = y[ip + k]
                pk_ = y[jp + k]
                g[ip + k] += -c3 * (
                    em5 * (-16.0 * d * pik + r * pk_)
                    + em7 * (w * pk_ + 70.0 * d * b_ * Kk)
                    - 10.0 * em7 * a_ * cc * Kk
                )
                g[ik + k] += -c3 * (
                    -10.0 * em7 * d * a_ * Kk - 10.0 * em7 * (b_ * cc * Kk + a_ * b_ * pk_)
                )
                g[jk + k] += -c3 * (
                    em7 * d * (-10.0 * a_ * Kik + 70.0 * b_ * pik)
                    - 10.0 * em7 * (b_ * cc * Kik + a_ * cc * pik)
                )
                g[jp + k] += -c3 * (
                    r * (5.0 * em7 * d * pk_ + em5 * pik)
                    + w * (7.0 * em9 * d * pk_ + em7 * pik)
                    - 10.0 * a_ * b_ * (7.0 * em9 * cc * pk_ + em7 * Kik)
                )
    for x in range(n_orb):
        ik = 6 * x
        L = prec[x, 0]
        n = prec[x, 1]
        alpha = prec[x, 2]
        beta = prec[x, 3]
        gamma = prec[x, 4]
        ex = math.sqrt(_dot(y, ik, ik))
        I1 = -2.625 / 9.0 * ex**-9 + 1.75 / 7.0 * ex**-7 - 0.125 / 5.0 * ex**-5
        dI1 = (2.625 - 1.75 * ex**2 + 0.125 * ex**4) * ex**-10
        I2 = -(ex**-3) / 3.0
        dI2 = ex**-4
        H += -alpha / ex + L * n * (beta * I1 + gamma * I2)
        dH = alpha / ex**2 + L * n * (beta * dI1 + gamma * dI2)
        for k in range(3):
            g[ik + k] += dH * y[ik + k] / ex
    for x in range(n_orb):
        ik = 6 * x
        ip = ik + 3
        L = prec[x, 0]
        for k in range(3):
            k1 = (k + 1) % 3
            k2 = (k + 2) % 3
            p_gp = y[ip + k1] * g[ip + k2] - y[ip + k2] * g[ip + k1]
            K_gK = y[ik + k1] * g[ik + k2] - y[ik + k2] * g[ik + k1]
            K_gp = y[ik + k1] * g[ip + k2] - y[ik + k2] * g[ip + k1]
            p_gK = y[ip + k1] * g[ik + k2] - y[ip + k2] * g[ik + k1]
            out[ik + k] = -(p_gp + K_gK) / L
            out[ip + k] = -(K_gp + p_gK) / L
    return H, out


@njit(cache=True)
def leconte_pair(G, m0, m1, r0, r1, k0, k1, g0, g1, lag0, lag1, a, e, w0, w1):
    """Constant-time-lag tidal rates of a binary of two stars with aligned spins (Hut 1981;
    Leconte et al. 2010; the equations of swingersim.tides.leconte_rates):
    (da/dt, de/dt, dspin0/dt, dspin1/dt).

    One rule for every caller: a star with lag <= 0 does not dissipate (no contribution, its
    spin constant); a star with dissipation but no moment of inertia (rg2 <= 0) is taken in
    the limit rg2 -> 0, i.e. rotating pseudo-synchronously (spin/n = N/Omega, Hut 1981, the
    equilibrium of its spin equation), and its spin rate in the state stays constant."""
    n = math.sqrt(G * (m0 + m1) / a**3)
    e2 = e * e
    bb = 1.0 - e2
    N = (1 + 7.5 * e2 + 5.625 * e2**2 + 0.3125 * e2**3) / bb**6
    Na1 = (1 + 15.5 * e2 + 31.875 * e2**2 + 11.5625 * e2**3 + 0.390625 * e2**4) / bb**7.5
    Om_e = (1 + 1.5 * e2 + 0.125 * e2**2) / bb**5
    Ne1 = (1 + 3.75 * e2 + 1.875 * e2**2 + 0.078125 * e2**3) / bb**6.5
    Om = (1 + 3 * e2 + 0.375 * e2**2) / bb**4.5
    da = 0.0
    de = 0.0
    dw0 = 0.0
    dw1 = 0.0
    for i in range(2):
        if i == 0:
            mi, mj, r, k2, rg2, lag, w = m0, m1, r0, k0, g0, lag0, w0
        else:
            mi, mj, r, k2, rg2, lag, w = m1, m0, r1, k1, g1, lag1, w1
        if lag <= 0.0:
            continue
        K = 1.5 * (2.0 * k2) * lag * (G * mj**2 / r) * (r / a) ** 6 * n**2
        C = rg2 * mi * r**2
        xs = w / n if C > 0.0 else N / Om
        da += 4.0 * a**2 / (G * m0 * m1) * K * (N * xs - Na1)
        de += 11.0 * a * e / (G * m0 * m1) * K * (Om_e * xs - 18.0 / 11.0 * Ne1)
        if C > 0.0:
            dw = -2.0 * K / (C * n) * (Om * xs - N)
            if i == 0:
                dw0 = dw
            else:
                dw1 = dw
    return da, de, dw0, dw1


@njit(cache=True)
def _add_e_rate(y, out, sk, de):
    """Add the change de/dt of the eccentricity of the orbit whose vectors start at sk (K at
    sk, pi at sk + 3), at constant directions of k and p."""
    e = math.sqrt(y[sk + 3] ** 2 + y[sk + 4] ** 2 + y[sk + 5] ** 2)
    if e > 1e-12:
        eta = math.sqrt(y[sk] ** 2 + y[sk + 1] ** 2 + y[sk + 2] ** 2)
        for k in range(3):
            out[sk + 3 + k] += de * y[sk + 3 + k] / e
            out[sk + k] += -(e / eta) * de * y[sk + k] / eta


@njit(cache=True)
def _lock_spins(y, out, sk, g_orb, s_spin):
    """Aligned spins locked to the orbit whose vectors start at sk (K at sk, e at sk + 3):
    the torque of the other orbits turns the orbit together with the spins of its stars, so
    the turning of K (the part of dK/dt perpendicular to K) is slowed by G / (G + S), G the
    orbital angular momentum L sqrt(1 - e^2) and S the spin angular momentum along it, and
    de/dt gets the component along k that keeps e perpendicular to K. The limit of fast
    spin-orbit precession (adiabatic: Anderson, Storch & Lai 2016, Sect. 4) of the
    conservative equations of Correia et al. (2011, eqs. 8-15): d(G + S)/dt = torque, total
    angular momentum conserved."""
    if s_spin <= 0.0 or g_orb <= 0.0:
        return
    eta = math.sqrt(y[sk] ** 2 + y[sk + 1] ** 2 + y[sk + 2] ** 2)
    if eta <= 0.0:
        return
    f = g_orb / (g_orb + s_spin)
    k0, k1, k2 = y[sk] / eta, y[sk + 1] / eta, y[sk + 2] / eta
    par = out[sk] * k0 + out[sk + 1] * k1 + out[sk + 2] * k2
    dk = np.empty(3)
    for i, kk in enumerate((k0, k1, k2)):
        perp = out[sk + i] - par * kk
        out[sk + i] = par * kk + f * perp
        dk[i] = f * perp / eta  # d(k hat)/dt
    e_dk = y[sk + 3] * dk[0] + y[sk + 4] * dk[1] + y[sk + 5] * dk[2]
    de_k = out[sk + 3] * k0 + out[sk + 4] * k1 + out[sk + 5] * k2
    corr = -e_dk - de_k
    out[sk + 3] += corr * k0
    out[sk + 4] += corr * k1
    out[sk + 5] += corr * k2


@njit(cache=True)
def dissipative_rhs(y, cst, bin_r, bin_k2, bin_rg2, bin_lag, bin_tides, a_out, gr, octupole,
                    bin_slaved):  # fmt: skip
    """Full right-hand side of the 2+2 system with tidal dissipation (state of 24: the 18
    vectorial elements, a_A, a_B, the four spins). Same equations as SecularTwoTwo.rhs
    (coefficients of _coeffs and _precession_terms, Leconte et al. 2010 tidal rates).

    cst = (G, c_light, m1, m2, m3, m4); bin_* have shape (2, 2): binary (A, B) x star;
    bin_slaved = 1: the spin is held pseudo-synchronous (backward runs, swingersim.tides.slaved_spin):
    the tidal rates as for rg2 -> 0, the spin in the state constant."""
    G, cl = cst[0], cst[1]
    m = cst[2:6]
    MA, MB = m[0] + m[1], m[2] + m[3]
    MAB = MA + MB
    muA, muB, mu = m[0] * m[1] / MA, m[2] * m[3] / MB, MA * MB / MAB
    XAa, XAb, XBa, XBb = m[0] / MA, m[1] / MA, m[2] / MB, m[3] / MB
    aA, aB = y[18], y[19]
    c = np.empty(4)
    c[0] = G * muA * MB / (8.0 * a_out) * (aA / a_out) ** 2
    c[1] = G * MA * muB / (8.0 * a_out) * (aB / a_out) ** 2
    c[2] = 15.0 / 64.0 * G * muA * MB / a_out * (aA / a_out) ** 3 * (XAa - XAb)
    c[3] = 15.0 / 64.0 * G * MA * muB / a_out * (aB / a_out) ** 3 * (XBb - XBa)
    prec = np.zeros((3, 5))
    for x in range(3):
        if x == 0:
            M, mu_x, a, b = MA, muA, aA, 0
        elif x == 1:
            M, mu_x, a, b = MB, muB, aB, 1
        else:
            M, mu_x, a, b = MAB, mu, a_out, -1
        L = mu_x * math.sqrt(G * M * a)
        n = math.sqrt(G * M / a**3)
        alpha = L * 3.0 * (G * M) ** 1.5 / (cl**2 * a**2.5) if gr else 0.0
        beta = 0.0
        gamma = 0.0
        if b >= 0 and bin_tides[b]:
            mb0 = m[2 * b]
            mb1 = m[2 * b + 1]
            for i in range(2):
                mi = mb0 if i == 0 else mb1
                mj = mb1 if i == 0 else mb0
                q = mj / mi
                x5 = bin_k2[b, i] * (bin_r[b, i] / a) ** 5
                beta += 15.0 * q * x5
                gamma += (1.0 + q) * (y[20 + 2 * b + i] / n) ** 2 * x5
        prec[x, 0] = L
        prec[x, 1] = n
        prec[x, 2] = alpha
        prec[x, 3] = beta
        prec[x, 4] = gamma
    _H, d18 = secular_rhs(y[:18].copy(), c, prec, octupole)
    out = np.zeros(24)
    out[:18] = d18
    for b in range(2):  # the spins locked to their binary (orbit b of the state)
        s_b = 0.0
        for i in range(2):
            s_b += bin_rg2[b, i] * m[2 * b + i] * bin_r[b, i] ** 2 * y[20 + 2 * b + i]
        sk = 6 * b
        eta = math.sqrt(y[sk] ** 2 + y[sk + 1] ** 2 + y[sk + 2] ** 2)
        _lock_spins(y, out, sk, prec[b, 0] * eta, s_b)
    for b in range(2):
        if not bin_tides[b] or max(bin_lag[b, 0], bin_lag[b, 1]) <= 0.0:
            continue
        sk = 6 * b
        e = math.sqrt(y[sk + 3] ** 2 + y[sk + 4] ** 2 + y[sk + 5] ** 2)
        da, de, dw0, dw1 = leconte_pair(
            G,
            m[2 * b],
            m[2 * b + 1],
            bin_r[b, 0],
            bin_r[b, 1],
            bin_k2[b, 0],
            bin_k2[b, 1],
            bin_rg2[b, 0] * (1.0 - bin_slaved[b, 0]),
            bin_rg2[b, 1] * (1.0 - bin_slaved[b, 1]),
            bin_lag[b, 0],
            bin_lag[b, 1],
            y[18 + b],
            e,
            y[20 + 2 * b],
            y[21 + 2 * b],
        )
        out[18 + b] = da
        out[20 + 2 * b] = dw0
        out[21 + 2 * b] = dw1
        _add_e_rate(y, out, sk, de)
    return out


@njit(cache=True)
def peters_pair(G, cl, m0, m1, a, e):
    """(da/dt, de/dt) of gravitational-wave emission (Peters 1964; swingersim.gw.peters_rates)."""
    M = m0 + m1
    k = G**3 * m0 * m1 * M / cl**5
    e2 = e * e
    b = 1.0 - e2
    da = -64.0 / 5.0 * k / (a**3 * b**3.5) * (1 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e2 * e2)
    de = -304.0 / 15.0 * e * k / (a**4 * b**2.5) * (1 + 121.0 / 304.0 * e2)
    return da, de


@njit(cache=True)
def nested_dissipative_rhs(y, cst, orb, a_idx, pairs, pair_f, tb_i, tb_f, d_i, d_f):
    """Right-hand side of a nested hierarchy with tidal dissipation (swingersim.secular.nested,
    SecularNested._rhs_dissipative): the coefficients of SecularNested.coefficients at the
    semi-major axes and spins of the state, the pair terms of secular_rhs_pairs and, in every
    dissipative two-star orbit, the Leconte et al. (2010) tidal rates and (optionally) the
    Peters (1964) gravitational-wave rates. The arguments are packed into a few arrays (the
    call overhead grows with their number):

    cst = (G, c_light, gr, octupole, gw); orb (n_orb, 3): masses of the two components and
    the semi-major axis in force of every orbit; a_idx: index of a in the state (-1: orb[:, 2]);
    pairs (n_pair, 2): inner and outer orbit of every pair, pair_f (n_pair, 2): mass of the
    other component of the outer orbit and the signed mass-fraction difference of the inner;
    tb_i (n_tb, 2) and tb_f (n_tb, 5): the rows of SecularNested._two_body_arrays (orbit, index
    of the spin in the state or -1; mass of the star and of its companion, R, k2, spin);
    d_i (n_d, 6): dissipative orbits (orbit, index of a and of the first spin in the state,
    tides on, the two stars' spins held pseudo-synchronous: backward runs,
    swingersim.tides.slaved_spin), d_f (n_d, 10): masses, R, k2, rg2 and lag of their two stars."""
    G, cl = cst[0], cst[1]
    gr, octupole, gw = cst[2] != 0.0, cst[3] != 0.0, cst[4] != 0.0
    n_orb = orb.shape[0]
    n6 = 6 * n_orb
    a = orb[:, 2].copy()
    for q in range(n_orb):
        if a_idx[q] >= 0:
            a[q] = y[a_idx[q]]
    prec = np.zeros((n_orb, 5))
    mu = np.empty(n_orb)
    for q in range(n_orb):
        M = orb[q, 0] + orb[q, 1]
        mu[q] = orb[q, 0] * orb[q, 1] / M
        L = mu[q] * math.sqrt(G * M * a[q])
        n = math.sqrt(G * M / a[q] ** 3)
        prec[q, 0] = L
        prec[q, 1] = n
        prec[q, 2] = L * 3.0 * (G * M) ** 1.5 / (cl**2 * a[q] ** 2.5) if gr else 0.0
    for r in range(tb_i.shape[0]):
        q = tb_i[r, 0]
        w = y[tb_i[r, 1]] if tb_i[r, 1] >= 0 else tb_f[r, 4]
        qq = tb_f[r, 1] / tb_f[r, 0]
        x5 = tb_f[r, 3] * (tb_f[r, 2] / a[q]) ** 5
        prec[q, 3] += 15.0 * qq * x5
        prec[q, 4] += (1.0 + qq) * (w / prec[q, 1]) ** 2 * x5
    npair = pairs.shape[0]
    pk = np.empty(npair, dtype=np.int64)
    pl = np.empty(npair, dtype=np.int64)
    c2 = np.empty(npair)
    c3 = np.empty(npair)
    for p in range(npair):
        pk[p] = pairs[p, 0]
        pl[p] = pairs[p, 1]
        ak = a[pk[p]]
        al = a[pl[p]]
        c2[p] = G * mu[pk[p]] * pair_f[p, 0] * ak**2 / (8.0 * al**3)
        c3[p] = 15.0 / 64.0 * G * mu[pk[p]] * pair_f[p, 0] * ak**3 / al**4 * pair_f[p, 1]
    out = np.zeros(y.shape[0])
    _H, d = secular_rhs_pairs(y[:n6].copy(), pk, pl, c2, c3, prec, octupole)
    out[:n6] = d
    for i in range(d_i.shape[0]):  # the spins locked to their orbit (with tides on)
        if d_i[i, 3]:
            q, isp, f = d_i[i, 0], d_i[i, 2], d_f[i]
            sk = 6 * q
            s_q = f[6] * f[0] * f[2] ** 2 * y[isp] + f[7] * f[1] * f[3] ** 2 * y[isp + 1]
            eta = math.sqrt(y[sk] ** 2 + y[sk + 1] ** 2 + y[sk + 2] ** 2)
            _lock_spins(y, out, sk, prec[q, 0] * eta, s_q)
    for i in range(d_i.shape[0]):
        q = d_i[i, 0]
        ia = d_i[i, 1]
        isp = d_i[i, 2]
        sk = 6 * q
        e = math.sqrt(y[sk + 3] ** 2 + y[sk + 4] ** 2 + y[sk + 5] ** 2)
        f = d_f[i]
        da = 0.0
        de = 0.0
        dw0 = 0.0
        dw1 = 0.0
        if d_i[i, 3]:
            da, de, dw0, dw1 = leconte_pair(
                G, f[0], f[1], f[2], f[3], f[4], f[5], f[6] * (1 - d_i[i, 4]),
                f[7] * (1 - d_i[i, 5]), f[8], f[9], a[q], e,
                y[isp], y[isp + 1],
            )  # fmt: skip
        if gw:
            ga, ge = peters_pair(G, cl, f[0], f[1], a[q], e)
            da, de = da + ga, de + ge
        out[ia] = da
        out[isp] = dw0
        out[isp + 1] = dw1
        _add_e_rate(y, out, sk, de)
    return out
