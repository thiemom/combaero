"""EjectorElement: 3-port supersonic ejector on the momentum-CV junction
topology.

Physics is Huang (1999) Eqs. 1-8 (entrainment ratio) + Kracik & Dvorak
(2016) Eqs. 7-13 (critical back pressure) -- see
combaero.network._ejector_huang1999 and validation/ejector/README.md for
the derivation, paper references, and accuracy comparison against
alternatives. This file tests the ELEMENT (residuals/diagnostics/topology/
solver integration), not the underlying physics functions, which are
already covered by python/tests/test_ejector_huang.py against Huang's
digitized Table 3/4 data.

References: docs/ejector/Huang_1d_analysis_ejector.pdf,
docs/ejector/Development_of_an_Analytical_Method_for.pdf
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
import combaero.species as sp
from combaero.network import ChannelElement, FlowNetwork, MomentumChamberNode, PressureBoundary
from combaero.network.components import NetworkMixtureState
from combaero.network.ejector_element import EjectorElement, choke_plane_gamma
from combaero.network.solver import NetworkSolver

_DRY_AIR_Y = list(cb.mole_to_mass(sp.dry_air()))

_THROAT_AREA = 3.14e-5  # ~6.3 mm diameter
_NOZZLE_EXIT_AREA = 1.0e-4
_MIXING_AREA = 8.0e-4
_CHANNEL_D = 0.05


def _states() -> list[NetworkMixtureState]:
    primary = NetworkMixtureState(
        P=700000.0, Pt=700000.0, T=440.0, Tt=440.0, m_dot=0.0, Y=_DRY_AIR_Y
    )
    secondary = NetworkMixtureState(
        P=22850.0, Pt=22850.0, T=280.0, Tt=280.0, m_dot=0.0, Y=_DRY_AIR_Y
    )
    outlet = NetworkMixtureState(P=40000.0, Pt=40000.0, T=300.0, Tt=300.0, m_dot=0.0, Y=_DRY_AIR_Y)
    return [primary, secondary, outlet]


def _element(**overrides) -> EjectorElement:
    kwargs = {
        "id": "ej1",
        "primary_node": "p",
        "secondary_node": "s",
        "outlet_node": "o",
        "throat_area": _THROAT_AREA,
        "nozzle_exit_area": _NOZZLE_EXIT_AREA,
        "mixing_area": _MIXING_AREA,
    }
    kwargs.update(overrides)
    return EjectorElement(**kwargs)


# ---------------------------------------------------------------------------
# Construction / validation
# ---------------------------------------------------------------------------


def test_construction_sets_area_ratios():
    e = _element()
    assert e.area_ratio_nozzle == pytest.approx(_NOZZLE_EXIT_AREA / _THROAT_AREA)
    assert e.area_ratio_mix == pytest.approx(_MIXING_AREA / _THROAT_AREA)
    assert e.geometry.area_ratio_nozzle == e.area_ratio_nozzle
    assert e.geometry.area_ratio_mix == e.area_ratio_mix


def test_construction_rejects_nonpositive_throat_area():
    with pytest.raises(ValueError, match="throat_area"):
        _element(throat_area=0.0)


def test_construction_rejects_nozzle_exit_not_exceeding_throat():
    with pytest.raises(ValueError, match="nozzle_exit_area"):
        _element(nozzle_exit_area=_THROAT_AREA * 0.5)


def test_construction_rejects_mixing_area_not_exceeding_nozzle_exit():
    with pytest.raises(ValueError, match="mixing_area"):
        _element(mixing_area=_NOZZLE_EXIT_AREA * 0.5)


def test_construction_rejects_nonpositive_recovery_efficiency():
    with pytest.raises(ValueError, match="recovery_efficiency"):
        _element(recovery_efficiency=0.0)


def test_ports_are_primary_then_secondary_then_outlet():
    e = _element()
    assert e.inlet_nodes == ["p", "s"]
    assert e.outlet_nodes == ["o"]
    assert e.port_nodes == ["p", "s", "o"]
    assert e._port_signs == [-1.0, -1.0, 1.0]


def test_row_scale_kinds_is_mdot_mdot_mdot_p():
    # Opposite of the base class's own N-pressure-rows + 1-mass-row pattern
    # -- see solver.py's scaling block and the module docstring.
    assert _element().row_scale_kinds() == ["mdot", "mdot", "mdot", "p"]


# ---------------------------------------------------------------------------
# Residuals: unit-level, against a hand-guessed (not necessarily converged) state
# ---------------------------------------------------------------------------


def test_residuals_cross_check_against_reference_physics():
    """residuals() must reproduce the SAME numbers as calling the reference
    physics functions directly with the live-evaluated gamma/R -- this is
    the cross-check against combaero.network._ejector_huang1999 (which is
    itself validated against Huang's Table 3/4 in test_ejector_huang.py)."""
    from combaero.network._ejector_huang1999 import (
        ETA_P,
        choked_mass_flow,
        critical_back_pressure,
        entrainment_ratio,
    )

    e = _element()
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    primary, secondary, outlet = _states()

    mp_guess, ms_guess = 0.5, 0.3
    mdot_out_guess = mp_guess + ms_guess
    port_mdots = [-mp_guess, -ms_guess, mdot_out_guess]
    P_jct_guess = 100000.0

    res, jac = e.residuals([primary, secondary, outlet], P_jct_guess, port_mdots)
    assert len(res) == 4
    assert set(jac.keys()) == {2}  # only the mass row is analytic

    gamma = choke_plane_gamma(secondary.Tt, secondary.X)
    r_gas = cb.specific_gas_constant(secondary.X)
    geom = e.geometry
    mp_choked = choked_mass_flow(primary.Pt, primary.Tt, e.throat_area, gamma, r_gas, eta=ETA_P)
    omega = entrainment_ratio(primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma).omega
    pc = critical_back_pressure(
        primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma, r_gas
    ).p_c_pa

    assert res[0] == pytest.approx(mp_guess - mp_choked)
    assert res[1] == pytest.approx(ms_guess - omega * mp_guess)
    assert res[2] == pytest.approx(0.0, abs=1e-9)  # mdot_out_guess == mp+ms by construction
    assert res[3] == pytest.approx(P_jct_guess - pc)


def test_residuals_all_zero_at_the_converged_point():
    """Solve the 4 equations directly (no full network) for the point where
    all residuals vanish, then verify residuals() actually returns zero
    there -- a tighter check than the cross-check above, which only proves
    internal consistency, not that R=0 has the physically expected root."""
    from combaero.network._ejector_huang1999 import (
        ETA_P,
        choked_mass_flow,
        critical_back_pressure,
        entrainment_ratio,
    )

    e = _element()
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    primary, secondary, outlet = _states()

    gamma = choke_plane_gamma(secondary.Tt, secondary.X)
    r_gas = cb.specific_gas_constant(secondary.X)
    geom = e.geometry
    mp = choked_mass_flow(primary.Pt, primary.Tt, e.throat_area, gamma, r_gas, eta=ETA_P)
    omega = entrainment_ratio(primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma).omega
    ms = omega * mp
    mdot_out = mp + ms
    pc = critical_back_pressure(
        primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma, r_gas
    ).p_c_pa

    port_mdots = [-mp, -ms, mdot_out]
    res, _jac = e.residuals([primary, secondary, outlet], pc, port_mdots)
    for r in res:
        assert r == pytest.approx(0.0, abs=1e-6)


def test_diagnostics_reports_omega_and_critical_mode():
    e = _element()
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    primary, secondary, outlet = _states()
    port_mdots = [-0.5, -0.3, 0.8]
    diag = e.diagnostics([primary, secondary, outlet], 100000.0, port_mdots)
    assert 0.0 < diag["omega"] < 5.0
    assert diag["gamma"] == pytest.approx(1.4, abs=0.05)  # air
    assert diag["r_gas"] == pytest.approx(287.0, abs=1.0)  # air
    assert diag["m_dot_primary"] == pytest.approx(0.5)
    assert diag["m_dot_secondary"] == pytest.approx(0.3)
    assert diag["m_dot_outlet"] == pytest.approx(0.8)
    # outlet Pt (40 kPa) vs P_jct guess (100 kPa): still "critical" by this metric.
    assert diag["critical_mode"] == 1.0


def test_verify_solution_consistent_flags_subcritical():
    e = _element()
    sol_ok = {"ej1.P_jct": 100000.0, "o.Pt": 40000.0}
    sol_bad = {"ej1.P_jct": 40000.0, "o.Pt": 100000.0}
    assert e.verify_solution_consistent(sol_ok) is True
    assert e.verify_solution_consistent(sol_bad) is False


def test_verify_solution_consistent_true_on_missing_keys():
    e = _element()
    assert e.verify_solution_consistent({}) is True


# ---------------------------------------------------------------------------
# choke_plane_gamma: sanity against known air properties
# ---------------------------------------------------------------------------


def test_choke_plane_gamma_matches_known_air_value():
    X = sp.dry_air()
    gamma = choke_plane_gamma(281.15, X)
    assert gamma == pytest.approx(1.4, abs=0.02)


# ---------------------------------------------------------------------------
# Full network integration: build a real graph, solve it, check the result.
# ---------------------------------------------------------------------------


def _expected_operating_point() -> dict[str, float]:
    """mp / ms / mdot_out / P_c* for the network's boundary conditions, from
    the same reference physics the element itself uses -- a good Newton
    starting point for a system with no analytic Jacobian yet (see the
    FD-fallback note in ejector_element.py's module docstring). Port static
    P/Pt are seeded at the boundary Pt (negligible drop expected over the
    short, large-diameter test channels).
    """
    from combaero.network._ejector_huang1999 import (
        ETA_P,
        EjectorGeometry,
        choked_mass_flow,
        critical_back_pressure,
        entrainment_ratio,
    )

    geom = EjectorGeometry(
        area_ratio_nozzle=_NOZZLE_EXIT_AREA / _THROAT_AREA,
        area_ratio_mix=_MIXING_AREA / _THROAT_AREA,
    )
    gamma = choke_plane_gamma(280.0, sp.dry_air())
    r_gas = cb.specific_gas_constant(sp.dry_air())
    mp = choked_mass_flow(700000.0, 440.0, _THROAT_AREA, gamma, r_gas, eta=ETA_P)
    omega = entrainment_ratio(700000.0, 440.0, 22850.0, 280.0, geom, gamma).omega
    ms = omega * mp
    pc = critical_back_pressure(700000.0, 440.0, 22850.0, 280.0, geom, gamma, r_gas).p_c_pa
    return {"mp": mp, "ms": ms, "mout": mp + ms, "pc": pc}


def _build_network() -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_primary", Pt=700000.0, Tt=440.0, Y=_DRY_AIR_Y))
    net.add_node(PressureBoundary("pb_secondary", Pt=22850.0, Tt=280.0, Y=_DRY_AIR_Y))
    net.add_node(PressureBoundary("pb_outlet", Pt=40000.0, Tt=300.0, Y=_DRY_AIR_Y))

    op = _expected_operating_point()

    mc_primary = MomentumChamberNode("mc_primary", area=math.pi * (_CHANNEL_D / 2.0) ** 2)
    mc_primary.initial_guess = {"mc_primary.P": 700000.0, "mc_primary.Pt": 700000.0}
    net.add_node(mc_primary)
    mc_secondary = MomentumChamberNode("mc_secondary", area=math.pi * (_CHANNEL_D / 2.0) ** 2)
    mc_secondary.initial_guess = {"mc_secondary.P": 22850.0, "mc_secondary.Pt": 22850.0}
    net.add_node(mc_secondary)
    mc_outlet = MomentumChamberNode("mc_outlet", area=math.pi * (_CHANNEL_D / 2.0) ** 2)
    mc_outlet.initial_guess = {"mc_outlet.P": 40000.0, "mc_outlet.Pt": 40000.0}
    net.add_node(mc_outlet)

    ch_primary = ChannelElement(
        "ch_primary",
        "pb_primary",
        "mc_primary",
        length=0.2,
        diameter=_CHANNEL_D,
        regime="incompressible",
    )
    ch_primary.initial_guess = {"ch_primary.m_dot": op["mp"]}
    net.add_element(ch_primary)

    ch_secondary = ChannelElement(
        "ch_secondary",
        "pb_secondary",
        "mc_secondary",
        length=0.2,
        diameter=_CHANNEL_D,
        regime="incompressible",
    )
    ch_secondary.initial_guess = {"ch_secondary.m_dot": op["ms"]}
    net.add_element(ch_secondary)

    ch_outlet = ChannelElement(
        "ch_outlet",
        "mc_outlet",
        "pb_outlet",
        length=0.2,
        diameter=_CHANNEL_D,
        regime="incompressible",
    )
    ch_outlet.initial_guess = {"ch_outlet.m_dot": op["mout"]}
    net.add_element(ch_outlet)

    ej1 = EjectorElement(
        id="ej1",
        primary_node="mc_primary",
        secondary_node="mc_secondary",
        outlet_node="mc_outlet",
        throat_area=_THROAT_AREA,
        nozzle_exit_area=_NOZZLE_EXIT_AREA,
        mixing_area=_MIXING_AREA,
    )
    ej1.initial_guess = {"ej1.P_jct": op["pc"]}
    net.add_element(ej1)
    return net


def test_full_network_solves_and_conserves_mass():
    net = _build_network()
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0)
    assert sol["__success__"], sol.get("__message__")

    mp = sol["ch_primary.m_dot"]
    ms = sol["ch_secondary.m_dot"]
    mout = sol["ch_outlet.m_dot"]
    assert mp > 0.0
    assert ms > 0.0
    assert mout == pytest.approx(mp + ms, rel=1e-4)

    # Entrainment ratio should be a plausible value, not a degenerate root.
    omega = ms / mp
    assert 0.0 < omega < 3.0

    # Critical mode: the actual outlet back pressure (imposed by pb_outlet,
    # 40 kPa) must not exceed the ejector's own P_c* (ej1.P_jct).
    assert sol["mc_outlet.Pt"] <= sol["ej1.P_jct"] * 1.01
