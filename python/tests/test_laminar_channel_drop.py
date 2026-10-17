"""A channel's laminar drop is Hagen-Poiseuille, linear in its flow.

The network friction kernels floored the laminar branch at Re 64 with the
turbulent fits' guard, so below that f saturated at 1 and the drop went as
m^2, with zero slope at zero flow; Colebrook had no laminar branch at all.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero import _solver_tools

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
L, D, P, TT = 0.5, 0.05, 2e5, 300.0


def _drop(m: float, model: str) -> tuple[float, float]:
    r = _solver_tools.channel_residuals_and_jacobian(
        m, P, P, TT, Y_AIR, P - 1e3, L, D, 0.0, model, 1.0, False
    )
    return float(r.dP_calc), float(r.d_dP_d_mdot)


@pytest.mark.parametrize("model", ["haaland", "serghides", "colebrook"])
def test_the_low_flow_drop_is_hagen_poiseuille(model: str) -> None:
    X = cb.mass_to_mole(Y_AIR)
    rho = cb.density(TT, P, X)
    mu = cb.transport_state(TT, P, X).mu
    A = math.pi * D * D / 4.0
    slope_hp = 32.0 * mu * L / (rho * A * D * D)  # dP/dm, f Re = 64
    for m in (1e-8, 1e-6, 1e-4):
        dP, dP_dm = _drop(m, model)
        Re = m * D / (mu * A)
        assert Re < 300.0
        assert dP == pytest.approx(slope_hp * m, rel=2e-3), (model, m)
        assert dP_dm == pytest.approx(slope_hp, rel=2e-3), (model, m)
