"""The Cd(Re) consistency solve of a discharge-hole orifice returns a root (#861).

OrificeElement solves m* = Cd(Re(m*)) |K| itself for a discharge-hole
correlation (#582), and the row is m - m*. Where |K| dCd/da sits just under 1
(Lichtarowicz, a long hole at low Re, Cd rising steeply), the loop takes the
contracting fixed-point step at a rate of 1 - g' per iteration. Its 50
iterations then end short of the root, and the stop was returned silently.
On the stress family `random` with Lichtarowicz forced on its orifices, 73
such stops were more than 1e-6 off (up to 1.7%). The state below is one of
them, from seed 4: m* was 2.3356e-5 kg/s and the root is 1.9346e-5.
"""

from __future__ import annotations

import copy
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
    components,
)

Y_AIR = list(cb.species.dry_air_mass())
CD_ITERATE = 0.35529840493010123


def _states() -> tuple[OrificeElement, object, object]:
    """The element and (state_in, state_out) objects, set to the stress iterate."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=210304.6, Tt=828.5, Y=Y_AIR))
    g.add_node(PressureBoundary("b", Pt=210304.6, Tt=537.9, Y=Y_AIR))
    g.add_element(
        OrificeElement(
            "o",
            "a",
            "b",
            diameter=0.3,
            correlation="Lichtarowicz",
            plate_thickness=0.6,
            regime="incompressible",
            in_pipe=False,
        )
    )
    s = NetworkSolver(g)
    s.network.resolve_all_topology()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s._build_x0()
    seen = []
    orig = components.OrificeElement._self_consistent_residuals

    def spy(self, si, so, flows, cd):
        seen.append((copy.copy(si), copy.copy(so)))
        return orig(self, si, so, flows, cd)

    components.OrificeElement._self_consistent_residuals = spy
    try:
        s._residuals_and_jacobian(np.array([1e-4]))
    finally:
        components.OrificeElement._self_consistent_residuals = orig
    si, so = seen[-1]
    si.Pt, si.P, si.Tt, si.m_dot = (
        210304.62477269676,
        210311.42965211833,
        828.5067017707306,
        6.0913e-4,
    )
    so.Pt, so.P, so.Tt = 210304.63561925475, 210317.3632033477, 537.8886721668993
    return g.elements["o"], si, so


def test_the_row_is_m_minus_a_root() -> None:
    el, si, so = _states()
    assert el._cd_follows_own_flow() and not el._approach_reads_flow()  # the path under test
    res, _ = el._self_consistent_residuals(si, so, {}, CD_ITERATE)
    m_star = abs(float(si.m_dot) - res[0])
    aK = abs(el._flow_and_partials(si, so, CD_ITERATE)[0] / CD_ITERATE)
    probe = copy.copy(si)
    probe.m_dot = m_star
    assert el._effective_Cd(probe, so, {}) * aK == pytest.approx(m_star, rel=1e-9)
    assert m_star == pytest.approx(1.9346e-5, rel=1e-4)
