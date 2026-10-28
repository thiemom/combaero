"""Exceptions the solver's safety net caught in the stress sweep, fixed at source.

- An orifice whose crossflow Cd underflows to exactly 0 (McGreehan-Schotsch
  at U1/Vi -> inf, where the hole's pressure drop vanishes) made the nozzle
  kernel refuse A_eff <= 0: 35 caught probes in the effusion family.
- A pressure loss's composition column took a central step of 0.5 y, which
  is 0 at a denormal mass fraction (5e-324), and divided by it: random 37.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

import combaero as cb
from combaero.network import NetworkSolver, OrificeElement
from combaero.network.components import _DM_DCD, NetworkMixtureState

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _states() -> tuple[NetworkMixtureState, NetworkMixtureState]:
    up = NetworkMixtureState(P=1.9e5, Pt=2e5, Tt=600.0, m_dot=0.1, Y=list(Y_AIR))
    dn = NetworkMixtureState(P=1.8e5, Pt=1.8e5, Tt=598.0, m_dot=0.1, Y=list(Y_AIR))
    return up, dn


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_a_vanishing_cd_is_a_vanishing_flow_not_an_exception(regime: str) -> None:
    """The kernels are linear in Cd, so a trace Cd is the Cd = 1 flow scaled,
    to the last bit; Cd = 0 is zero flow with zero partials. d(m)/d(Cd)
    (_DM_DCD) is the Cd = 1 flow throughout: the slope, not a partial that
    scales."""
    o = OrificeElement("o", "a", "b", Cd=0.6, area=1e-3, regime=regime)
    up, dn = _states()
    m1, d1 = o._kernel_from(up, dn, "a", "b", 1.0)
    assert d1.pop(_DM_DCD) == m1
    m0, d0 = o._kernel_from(up, dn, "a", "b", 0.0)
    assert d0.pop(_DM_DCD) == m1
    assert m0 == 0.0 and all(v == 0.0 for v in d0.values())
    cd = 1e-13
    m, d = o._kernel_from(up, dn, "a", "b", cd)
    assert d.pop(_DM_DCD) == m1
    assert m == m1 * cd
    assert d == {k: v * cd for k, v in d1.items()}
    # Above the threshold the kernel itself runs, and agrees with the scaling.
    cd = 1e-11
    m, _ = o._kernel_from(up, dn, "a", "b", cd)
    assert m == pytest.approx(m1 * cd, rel=1e-12)


def _caught(family: str, seed: int) -> tuple[dict, list[dict]]:
    net, _ = FAMILIES[family](seed)
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=60.0)
    return r, list(r.get("__caught_exceptions__") or [])


def test_the_effusion_liner_that_raised_a_zero_cd_now_converges_cleanly() -> None:
    """Stress effusion seed 5: one caught A_eff probe before, now none."""
    r, caught = _caught("effusion", 5)
    assert not caught, caught
    assert r["__success__"], r["__message__"]


def test_a_denormal_mass_fraction_no_longer_divides_by_zero() -> None:
    """Stress random seed 37: three ZeroDivisionErrors in a pressure loss's
    composition column before. Its non-convergence is a separate matter."""
    _, caught = _caught("random", 37)
    assert not [c for c in caught if c.get("type") == "ZeroDivisionError"], caught
