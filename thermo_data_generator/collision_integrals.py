"""Generate include/collision_integral_data.h: Omega*(2,2) and d ln Omega*/d ln T*.

Run from the repository root with the main environment (needs numpy/scipy):

    uv run python thermo_data_generator/collision_integrals.py

The delta* = 0 column (Lennard-Jones 12-6) is computed here by quadrature of
the classical collision integrals, NOT copied from the Monchick-Mason table:

  chi(E, b)  = pi - 2 b int_{r0}^inf dr / (r^2 sqrt(1 - b^2/r^2 - V*(r)/E))
  Q2*(E)     = 3 int_0^inf (1 - cos^2 chi) b db          (rigid sphere -> 1)
  Om*(T)     = 1/(6 T^4) int E^3 exp(-E/T) Q2*(E) dE
  dOm*/dT    = (-4 I3 + I4) / (6 T),  In = int_0^inf x^n exp(-x) Q2*(T x) dx

(reduced units: r in sigma, E and V* = 4(r^-12 - r^-6) in epsilon). The
outermost turning point r0 is taken, and the impact-parameter integral is split
at the orbiting b below the orbiting energy. Measured: converged to <= 1e-5
in Q2* (inner rule x5, tolerance x100), and Om* agrees with Kim & Monroe
(2014, J. Comput. Phys. 273:358) to <= 5e-5 at every node in their range
(0.3 <= T* <= 400). The Monchick-Mason delta* = 0 column, by contrast, is
off by up to 1.3e-3 at T* = 0.3 and 6.1e-3 at T* = 100 (#485).

The delta* > 0 (Stockmayer) columns are the Monchick-Mason (1961) values as
transcribed by Cantera (MMCollisionInt.cpp); their slopes are ESTIMATED with
a natural cubic spline in (ln T*, ln Omega), the best estimated-slope scheme
on the delta* = 0 column (0.093% max, against 0.032% for exact slopes). No
exact Stockmayer reference is computed here.
"""

from __future__ import annotations

import math
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import integrate, interpolate, optimize

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "include" / "collision_integral_data.h"

# Monchick-Mason axes and table (Cantera MMCollisionInt.cpp transcription).
T_STAR = [
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    1.2,
    1.4,
    1.6,
    1.8,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    5.0,
    6.0,
    7.0,
    8.0,
    9.0,
    10.0,
    12.0,
    14.0,
    16.0,
    18.0,
    20.0,
    25.0,
    30.0,
    35.0,
    40.0,
    50.0,
    75.0,
    100.0,
]
DELTA = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]
MM_OMEGA22 = [
    4.1005,
    4.266,
    4.833,
    5.742,
    6.729,
    8.624,
    10.34,
    11.89,
    3.2626,
    3.305,
    3.516,
    3.914,
    4.433,
    5.57,
    6.637,
    7.618,
    2.8399,
    2.836,
    2.936,
    3.168,
    3.511,
    4.329,
    5.126,
    5.874,
    2.531,
    2.522,
    2.586,
    2.749,
    3.004,
    3.64,
    4.282,
    4.895,
    2.2837,
    2.277,
    2.329,
    2.46,
    2.665,
    3.187,
    3.727,
    4.249,
    2.0838,
    2.081,
    2.13,
    2.243,
    2.417,
    2.862,
    3.329,
    3.786,
    1.922,
    1.924,
    1.97,
    2.072,
    2.225,
    2.614,
    3.028,
    3.435,
    1.7902,
    1.795,
    1.84,
    1.934,
    2.07,
    2.417,
    2.788,
    3.156,
    1.6823,
    1.689,
    1.733,
    1.82,
    1.944,
    2.258,
    2.596,
    2.933,
    1.5929,
    1.601,
    1.644,
    1.725,
    1.838,
    2.124,
    2.435,
    2.746,
    1.4551,
    1.465,
    1.504,
    1.574,
    1.67,
    1.913,
    2.181,
    2.451,
    1.3551,
    1.365,
    1.4,
    1.461,
    1.544,
    1.754,
    1.989,
    2.228,
    1.28,
    1.289,
    1.321,
    1.374,
    1.447,
    1.63,
    1.838,
    2.053,
    1.2219,
    1.231,
    1.259,
    1.306,
    1.37,
    1.532,
    1.718,
    1.912,
    1.1757,
    1.184,
    1.209,
    1.251,
    1.307,
    1.451,
    1.618,
    1.795,
    1.0933,
    1.1,
    1.119,
    1.15,
    1.193,
    1.304,
    1.435,
    1.578,
    1.0388,
    1.044,
    1.059,
    1.083,
    1.117,
    1.204,
    1.31,
    1.428,
    0.99963,
    1.004,
    1.016,
    1.035,
    1.062,
    1.133,
    1.22,
    1.319,
    0.96988,
    0.9732,
    0.983,
    0.9991,
    1.021,
    1.079,
    1.153,
    1.236,
    0.92676,
    0.9291,
    0.936,
    0.9473,
    0.9628,
    1.005,
    1.058,
    1.121,
    0.89616,
    0.8979,
    0.903,
    0.9114,
    0.923,
    0.9545,
    0.9955,
    1.044,
    0.87272,
    0.8741,
    0.878,
    0.8845,
    0.8935,
    0.9181,
    0.9505,
    0.9893,
    0.85379,
    0.8549,
    0.858,
    0.8632,
    0.8703,
    0.8901,
    0.9164,
    0.9482,
    0.83795,
    0.8388,
    0.8414,
    0.8456,
    0.8515,
    0.8678,
    0.8895,
    0.916,
    0.82435,
    0.8251,
    0.8273,
    0.8308,
    0.8356,
    0.8493,
    0.8676,
    0.8901,
    0.80184,
    0.8024,
    0.8039,
    0.8065,
    0.8101,
    0.8201,
    0.8337,
    0.8504,
    0.78363,
    0.784,
    0.7852,
    0.7872,
    0.7899,
    0.7976,
    0.8081,
    0.8212,
    0.76834,
    0.7687,
    0.7696,
    0.7712,
    0.7733,
    0.7794,
    0.7878,
    0.7983,
    0.75518,
    0.7554,
    0.7562,
    0.7575,
    0.7592,
    0.7642,
    0.7711,
    0.7797,
    0.74364,
    0.7438,
    0.7445,
    0.7455,
    0.747,
    0.7512,
    0.7569,
    0.7642,
    0.71982,
    0.72,
    0.7204,
    0.7211,
    0.7221,
    0.725,
    0.7289,
    0.7339,
    0.70097,
    0.7011,
    0.7014,
    0.7019,
    0.7026,
    0.7047,
    0.7076,
    0.7112,
    0.68545,
    0.6855,
    0.6858,
    0.6861,
    0.6867,
    0.6883,
    0.6905,
    0.6932,
    0.67232,
    0.6724,
    0.6726,
    0.6728,
    0.6733,
    0.6743,
    0.6762,
    0.6784,
    0.65099,
    0.651,
    0.6512,
    0.6513,
    0.6516,
    0.6524,
    0.6534,
    0.6546,
    0.61397,
    0.6141,
    0.6143,
    0.6145,
    0.6147,
    0.6148,
    0.6148,
    0.6147,
    0.5887,
    0.5889,
    0.5894,
    0.59,
    0.5903,
    0.5901,
    0.5895,
    0.5885,
]

# Kim & Monroe (2014) interpolation for Omega*(2,2), 0.3 <= T* <= 400, with
# the scaling the authors' coefficient table uses (as in the chapensk package).
KM_A = -0.92032979
KM_B = [
    2.3508044,
    0.50110649,
    -0.1 * 4.7193769,
    0.1 * 1.5806367,
    -0.01 * 2.6367184,
    0.001 * 1.8120118,
]
KM_C = [
    1.6330213,
    -0.1 * 6.9795156,
    0.1 * 1.6096572,
    -0.01 * 2.2109440,
    0.001 * 1.7031434,
    -0.0001 * 0.56699986,
]


def kim_monroe_omega22(T: float) -> float:
    return KM_A + sum(KM_B[k] / T ** (k + 1) + KM_C[k] * math.log(T) ** (k + 1) for k in range(6))


# --- classical Lennard-Jones collision integrals -----------------------------

_XG, _WG = np.polynomial.legendre.leggauss(400)


def _V(r):
    return 4.0 * (r**-12 - r**-6)


def _F(r, E, b):
    return 1.0 - (b / r) ** 2 - _V(r) / E


def _r0_outer(E: float, b: float) -> float:
    hi = max(2.0 * b, 3.0) + 2.0
    while _F(hi, E, b) <= 0:
        hi *= 2
    rs = np.geomspace(hi, 0.3, 4000)
    f = _F(rs, E, b)
    for i in range(len(rs) - 1):
        if f[i] > 0 and f[i + 1] <= 0:
            return optimize.brentq(_F, rs[i + 1], rs[i], args=(E, b), xtol=1e-15, rtol=1e-15)
    raise RuntimeError("no turning point")


def _chi(E: float, b: float) -> float:
    if b == 0.0:
        return math.pi
    r0 = _r0_outer(E, b)
    s = 0.5 * (_XG + 1.0)
    u = 1.0 - s**2  # u = r0/r; removes the inverse-sqrt endpoint at u = 1
    g = np.maximum(_F(r0 / u, E, b), 1e-300)
    return math.pi - 2.0 * (b / r0) * float(np.dot(0.5 * _WG, 2.0 * s / np.sqrt(g)))


def _b_orbit(E: float) -> float | None:
    def peak(b):
        res = optimize.minimize_scalar(
            lambda r: -(E * b * b / r**2 + _V(r)), bounds=(1.0, 10.0), method="bounded"
        )
        return -res.fun - E

    try:
        if peak(0.5) * peak(10.0) > 0:
            return None
        return optimize.brentq(peak, 0.5, 10.0, xtol=1e-14)
    except (ValueError, RuntimeError):
        return None


def q2_star(E: float) -> float:
    import warnings

    warnings.simplefilter("ignore")
    f = lambda b: (1.0 - math.cos(_chi(E, b)) ** 2) * b  # noqa: E731
    bmax = max(6.0, 3.0 * E ** (-1.0 / 6.0) + 3.0)
    bo = _b_orbit(E) if E < 0.85 else None
    val, _ = integrate.quad(
        f, 0.0, bmax, points=[bo] if bo else None, limit=2000, epsabs=1e-13, epsrel=1e-10
    )
    tail, _ = integrate.quad(f, bmax, 4 * bmax, limit=400, epsabs=1e-15, epsrel=1e-9)
    return 3.0 * (val + tail)


def omega22_exact(T_nodes: list[float], workers: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """(Omega*(2,2), dOmega*/dT*) at each T*, from Q2* on an energy grid."""
    E = np.unique(np.concatenate([np.geomspace(1e-4, 1e4, 480), np.linspace(0.6, 1.0, 81)]))
    with Pool(workers) as p:
        Q = np.array(p.map(q2_star, E, chunksize=4))
    s = interpolate.PchipInterpolator(np.log(E), np.log(Q))

    def qf(e):
        return np.exp(s(np.log(e)))

    om, dom = [], []
    for T in T_nodes:
        lo, hi = math.log(1e-4 / T), math.log(min(1e4 / T, 60.0))
        f3 = lambda lx, T=T: math.exp(4 * lx - math.exp(lx)) * qf(T * math.exp(lx))  # noqa: E731
        f4 = lambda lx, T=T: math.exp(5 * lx - math.exp(lx)) * qf(T * math.exp(lx))  # noqa: E731
        I3 = integrate.quad(f3, lo, hi, limit=400, epsabs=0, epsrel=1e-11)[0]
        I4 = integrate.quad(f4, lo, hi, limit=400, epsabs=0, epsrel=1e-11)[0]
        om.append(I3 / 6.0)
        dom.append((-4.0 * I3 + I4) / (6.0 * T))
    return np.array(om), np.array(dom)


# Grid extension beyond Monchick-Mason's T* = 100, so hydrogen (eps 38 K)
# stays on exact nodes up to ~15000 K instead of a power-law extrapolation
# that overshoots by 0.5% at T* = 400.
T_EXTRA = [150.0, 200.0, 300.0, 400.0]


def build_tables() -> tuple[np.ndarray, np.ndarray, list[float]]:
    """(Omega table, d ln Omega / d ln T* table, T* axis), row-major [T*][delta*].

    delta* = 0: exact quadrature, values and slopes.
    delta* > 0: Monchick-Mason scaled by the delta* = 0 column's exact/MM
    ratio at the same T*. Their table carries the SAME error as its LJ
    column wherever the dipole contribution is small -- at T* = 100 every
    column sits 0.6% above the exact LJ value they all converge to -- so
    the correction keeps the columns mutually consistent (no 0.6% step in
    delta* between the exact and the tabulated columns). At low T*, where the
    dipole matters, the ratio is an assumption, but it changes values by at
    most 0.13% there. Beyond T* = 100 the polar columns are the exact LJ value
    times their (~1) ratio to it at T* = 100.
    """
    axis = T_STAR + T_EXTRA
    T = np.array(axis)
    mm = np.array(MM_OMEGA22).reshape(37, 8)
    om0, dom0 = omega22_exact(axis)
    km_dev = max(
        abs(o / kim_monroe_omega22(t) - 1) for t, o in zip(T, om0, strict=True) if 0.3 <= t <= 400
    )
    if km_dev > 1e-4:
        raise RuntimeError(f"quadrature disagrees with Kim & Monroe by {km_dev:.1e}")
    table = np.zeros((len(axis), 8))
    table[:, 0] = om0
    ratio = om0[:37] / mm[:, 0]
    for c in range(1, 8):
        table[:37, c] = mm[:, c] * ratio
        table[37:, c] = om0[37:] * (mm[36, c] / mm[36, 0])
    slopes = np.zeros_like(table)
    slopes[:, 0] = T * dom0 / om0
    lx = np.log(T)
    for c in range(1, 8):
        spl = interpolate.CubicSpline(lx, np.log(table[:, c]), bc_type="natural")
        slopes[:, c] = spl(lx, 1)
    return table, slopes, axis


def write_header(table: np.ndarray, slopes: np.ndarray, axis: list[float]) -> None:
    def rows(a: np.ndarray) -> str:
        return "\n".join("    " + ", ".join(f"{v:.10g}" for v in r) + "," for r in a)

    OUT.write_text(f"""#pragma once
// !! DO NOT EDIT THIS FILE MANUALLY !!
// Generated by thermo_data_generator/collision_integrals.py
// Regenerate with (repository root):
//   uv run python thermo_data_generator/collision_integrals.py
//
// Reduced collision integral Omega*(2,2)(T*, delta*) and its logarithmic slope
// d ln Omega* / d ln T* on the Monchick-Mason grid, for cubic Hermite
// interpolation in (ln T*, ln Omega*) (src/transport.cpp).
//
// delta* = 0: computed by quadrature of the Lennard-Jones (12-6) collision
//   integrals, values AND slopes; agrees with Kim & Monroe (2014) to <= 1e-4.
// delta* > 0: Monchick & Mason (1961), as transcribed by Cantera, scaled by
//   the delta* = 0 column's exact/tabulated ratio (consistent columns; exact
//   where the dipole contribution vanishes); slopes ESTIMATED with a natural
//   cubic spline in (ln T*, ln Omega*).
// T* grid: Monchick-Mason's 0.1-100 plus 150, 200, 300, 400 (exact LJ).

namespace combaero::collision {{

constexpr int kNumTStar = {len(axis)};
constexpr int kNumDelta = 8;

constexpr double kTStar[kNumTStar] = {{
    {", ".join(f"{v:g}" for v in axis)}}};

constexpr double kDelta[kNumDelta] = {{{", ".join(f"{v:g}" for v in DELTA)}}};

// Omega*(2,2), row-major [T*][delta*].
constexpr double kOmega22[kNumTStar * kNumDelta] = {{
{rows(table)}
}};

// d ln Omega*(2,2) / d ln T*, row-major [T*][delta*].
constexpr double kDlnOmega22DlnT[kNumTStar * kNumDelta] = {{
{rows(slopes)}
}};

}}  // namespace combaero::collision
""")


if __name__ == "__main__":
    t, s, axis = build_tables()
    write_header(t, s, axis)
    print(f"Wrote {OUT}")
