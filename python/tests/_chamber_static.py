"""Independent reference for a momentum chamber's static state (#576).

A chamber node carries its STAGNATION temperature T0. Its static T solves the
energy equation h(T0) = h(T) + v^2 / 2 with v = m / (rho(T, P) A); the
stagnation pressure is then Pt_from_static(P, T, M). Solved here by bracketing
on the energy equation and the library's iterative P0 -- a different route
from the closed form the C++ closure uses -- so a test comparing the two
checks the closure rather than restating it.
"""

from __future__ import annotations

from scipy.optimize import brentq

import combaero as cb


def chamber_static_T(P: float, T0: float, m: float, area: float, X: list[float]) -> float:
    if m == 0.0:
        return T0

    def energy(T: float) -> float:
        v = abs(m) / (float(cb.density(T, P, X)) * area)
        return float(cb.h_mass(T0, X)) - float(cb.h_mass(T, X)) - 0.5 * v * v

    return brentq(energy, 0.2 * T0, T0, xtol=1e-12, rtol=1e-14)


def chamber_isentropic_rise(P: float, T0: float, m: float, area: float, X: list[float]) -> float:
    """P0 - P of the chamber's static state: what its closure calls dynamic head."""
    T = chamber_static_T(P, T0, m, area, X)
    v = abs(m) / (float(cb.density(T, P, X)) * area)
    M = v / float(cb.speed_of_sound(T, X))
    return float(cb.Pt_from_static(P, T, M, X)) - P


def chamber_static_density(P: float, T0: float, m: float, area: float, X: list[float]) -> float:
    return float(cb.density(chamber_static_T(P, T0, m, area, X), P, X))
