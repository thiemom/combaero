"""The tuned constants in the MPCEv2 junction closure, and their switches.

Under this repo's policy a tuned constant is labelled at its definition and
must be switchable so it can be scored on and off against the validation
data. These tests guard three things: the constants exist with their
documented values (a silent edit is a retune without a table), the switches
are live (a measurement with a dead switch proves nothing), and the defaults
are behaviour-preserving (the faithful Mynard port must not move when the
switches are added).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb
from combaero.network import _mynard2010 as mynard
from combaero.network.components import NetworkMixtureState
from combaero.network.mpce_v2_element import MPCEv2Element

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_A = np.array([0.01, 0.01, 0.01])
_THETA = np.array([math.pi, 0.0, math.pi / 2.0])  # common, straight, lateral
_U_DIVIDING = np.array([10.0, -6.0, -4.0])  # supplier +, collectors -


def test_tuned_constants_carry_their_documented_values():
    """Mynard Eq 36 coefficients and the Matlab regulariser, as labelled.

    Changing any of these is a retune and needs an on/off table on #271
    before it lands -- not a silent edit.
    """
    assert mynard.MYNARD_ETA_A0 == 0.8
    assert mynard.MYNARD_ETA_A1 == -0.2
    assert mynard.FLOW_RATIO_DAMPING == 0.02
    assert MPCEv2Element.DEFAULT_JOINING_ETRANSFER_ALPHA == 0.2


def test_eta_scale_default_reproduces_the_faithful_port():
    """Adding the switch must not move the closure at its default."""
    explicit = mynard.junction_loss_coefficient(_U_DIVIDING, _A, _THETA, eta_scale=1.0)
    implicit = mynard.junction_loss_coefficient(_U_DIVIDING, _A, _THETA)

    assert implicit.K is not None and explicit.K is not None
    np.testing.assert_array_equal(explicit.K, implicit.K)
    np.testing.assert_array_equal(explicit.C, implicit.C)


def test_eta_switch_is_live():
    """eta_scale = 0 removes the energy-transfer term and must change K.

    A dividing tee at 90 deg has a non-zero eta on both collectors, so a dead
    switch would show up here as identical coefficients.
    """
    with_eta = mynard.junction_loss_coefficient(_U_DIVIDING, _A, _THETA, eta_scale=1.0)
    without = mynard.junction_loss_coefficient(_U_DIVIDING, _A, _THETA, eta_scale=0.0)

    assert with_eta.K is not None and without.K is not None
    assert not np.allclose(with_eta.K, without.K), "eta_scale=0 changed nothing"
    # With eta off the pseudosupplier area is the plain flow-weighted one.
    assert np.allclose(without.pseudo_area, without.pseudo_area[0])


def test_eta_off_matches_the_formula_with_a0_a1_zeroed():
    """eta_scale=0 must be exactly what setting a0 = a1 = 0 would give.

    Guards the switch against being implemented as anything other than a
    multiplier on the eta term itself.
    """
    off = mynard.junction_loss_coefficient(_U_DIVIDING, _A, _THETA, eta_scale=0.0)

    a0, a1 = mynard.MYNARD_ETA_A0, mynard.MYNARD_ETA_A1
    try:
        mynard.MYNARD_ETA_A0 = 0.0
        mynard.MYNARD_ETA_A1 = 0.0
        zeroed = mynard.junction_loss_coefficient(_U_DIVIDING, _A, _THETA, eta_scale=1.0)
    finally:
        mynard.MYNARD_ETA_A0, mynard.MYNARD_ETA_A1 = a0, a1

    assert zeroed.K is not None and off.K is not None
    np.testing.assert_allclose(off.K, zeroed.K, rtol=1e-12)


def _element(**kwargs) -> MPCEv2Element:
    element = MPCEv2Element.__new__(MPCEv2Element)
    element.id = "jct"
    element.N = 3
    element.port_nodes = ["p0", "p1", "p2"]
    element.port_areas = list(_A)
    element.port_angles_deg = [180.0, 0.0, 90.0]
    element._port_signs = [-1.0, 1.0, 1.0]
    element._port_element_ids = ["e0", "e1", "e2"]
    element.flow_direction = "branch"
    element.strict = False
    element.joining_etransfer_alpha = 0.2
    element.jacobian_method = "sympy"
    element.penalty_alpha = 0.0
    element.eta_scale = kwargs.get("eta_scale", 1.0)
    return element


def _states():
    return [
        NetworkMixtureState(P=1.0e5, Pt=100_300.0, T=300.0, Tt=300.5, m_dot=0.0, Y=_Y)
        for _ in range(3)
    ]


def test_element_plumbs_eta_scale_into_its_residual():
    """The switch has to reach the residual, not just the closure."""
    mdots = [-0.10, 0.06, 0.04]
    on, _ = _element(eta_scale=1.0).residuals(_states(), 100_300.0, mdots)
    off, _ = _element(eta_scale=0.0).residuals(_states(), 100_300.0, mdots)

    assert on[0] == pytest.approx(off[0])  # common row carries no K
    assert on[1] != pytest.approx(off[1]) or on[2] != pytest.approx(off[2])


def test_element_constructor_accepts_and_defaults_eta_scale():
    element = MPCEv2Element(
        id="jct",
        inlet_nodes=["a"],
        outlet_nodes=["b", "c"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0],
        port_areas=[0.01, 0.01, 0.01],
    )
    assert element.eta_scale == 1.0

    scaled = MPCEv2Element(
        id="jct",
        inlet_nodes=["a"],
        outlet_nodes=["b", "c"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0],
        port_areas=[0.01, 0.01, 0.01],
        eta_scale=0.0,
    )
    assert scaled.eta_scale == 0.0
