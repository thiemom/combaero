"""#892: FuelBoundary, a fuel controller -- the fuel flow that puts a
combustor's exit at a target total temperature, on the branch asked for."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    CombustorNode,
    FlowNetwork,
    FuelBoundary,
    LosslessConnectionElement,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.components import OrificeElement

AIR = list(cb.species.humid_air_mass(288.15, 101325.0, 0.6))
CH4 = list(cb.species.to_mass(cb.species.pure_species("CH4")))
T_AIR = 389.2


def _net(target: float, branch: str = "lean", **kw) -> FlowNetwork:
    n = FlowNetwork()
    n.add_node(PressureBoundary("air", Pt=2.536e5, Tt=T_AIR, Y=AIR))
    n.add_node(
        FuelBoundary("fuel", combustor="burner", Tt_target=target, branch=branch, Y=CH4, **kw)
    )
    n.add_node(CombustorNode("burner"))
    n.add_node(PressureBoundary("out", Pt=2.3e5, Tt=800.0, Y=AIR))
    n.add_element(
        OrificeElement(
            "inlet", "air", "burner", diameter=0.06, Cd=0.8, correlation="fixed", in_pipe=False
        )
    )
    n.add_element(LosslessConnectionElement("inj", "fuel", "burner"))
    n.add_element(
        OrificeElement(
            "exit", "burner", "out", diameter=0.05, Cd=0.8, correlation="fixed", in_pipe=False
        )
    )
    return n


def _solve(net: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return s, s.solve()


@pytest.mark.parametrize(
    ("target", "branch"),
    [(808.0, "lean"), (1500.0, "lean"), (2300.0, "lean"), (1500.0, "rich"), (2000.0, "rich")],
)
def test_the_target_is_met_on_the_branch_asked_for(target: float, branch: str) -> None:
    _, r = _solve(_net(target, branch))
    assert r["__success__"], r["__message__"]
    assert r["burner.Tt"] == pytest.approx(target, abs=1e-3)
    assert (r["fuel.phi"] < 1.0) if branch == "lean" else (r["fuel.phi"] > 1.0)
    assert r["fuel.m_dot"] > 0.0
    # the mass row, to solver tolerance
    assert r["inj.m_dot"] == pytest.approx(r["fuel.m_dot"], rel=1e-4)


def test_lean_is_the_default() -> None:
    assert FuelBoundary("f", combustor="c", Tt_target=1000.0).branch == "lean"


def test_its_fuel_burns_to_the_target_by_the_combustors_own_combustion() -> None:
    """Adiabatic, at a fixed air flow: the network's fuel, mixed with the air
    and burnt by complete_combustion (the combustor's model, smooth phi
    blending), gives the target. set_fuel_stream_for_Tad burns without the
    blending and lands within 1e-3 of it."""
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("air", m_dot=0.3, Tt=T_AIR, Y=AIR))
    n.add_node(FuelBoundary("fuel", combustor="burner", Tt_target=1400.0, Tt=300.0, Y=CH4))
    n.add_node(CombustorNode("burner"))
    n.add_node(PressureBoundary("out", Pt=2.0e5, Tt=1400.0, Y=AIR))
    n.add_element(LosslessConnectionElement("a_in", "air", "burner"))
    n.add_element(LosslessConnectionElement("inj", "fuel", "burner"))
    n.add_element(
        OrificeElement(
            "exit", "burner", "out", diameter=0.05, Cd=0.8, correlation="fixed", in_pipe=False
        )
    )
    _, r = _solve(n)
    assert r["__success__"], r["__message__"]

    def stream(T: float, Y: list[float], m: float) -> cb.Stream:
        st = cb.Stream()
        st.set_T(T)
        st.set_P(2.0e5)
        st.set_X(list(cb.mass_to_mole(Y)))
        st.set_mdot(m)
        return st

    mixed = cb.mix([stream(T_AIR, AIR, 0.3), stream(300.0, CH4, r["fuel.m_dot"])])
    burnt = cb.complete_combustion(mixed.T, list(mixed.X), mixed.P)
    assert pytest.approx(1400.0, abs=1e-3) == burnt.T
    ref = cb.set_fuel_stream_for_Tad(
        1400.0, stream(300.0, CH4, 0.0), stream(T_AIR, AIR, 0.3), tol=1e-6
    )
    assert r["fuel.m_dot"] == pytest.approx(ref.mdot, rel=1e-3)


def test_the_other_root_is_refused_with_the_reason() -> None:
    """808 K has no rich root within phi 10 (complete combustion): the solve
    finds the lean one and refuses it, saying why."""
    _, r = _solve(_net(808.0, "rich"))
    assert not r["__success__"]
    assert "not the rich one asked for" in r["__message__"]
    assert "No rich root at the start" in r["__message__"]


def test_a_target_below_the_oxidizer_is_refused() -> None:
    _, r = _solve(_net(350.0))
    assert not r["__success__"]
    assert "not above the oxidizer's 389.2 K" in r["__message__"]


def test_a_target_above_the_peak_is_refused() -> None:
    _, r = _solve(_net(2500.0))
    assert not r["__success__"]
    assert "above the peak" in r["__message__"]


def test_its_rows_are_exact() -> None:
    """Every Jacobian row against central differences at the root, columns
    scaled by |x| (the fuel's ln m_dot column through the combustor's Tt
    relay and the mass rows)."""
    s, r = _solve(_net(1500.0))
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    j_fuel = s.unknown_names.index("fuel.ln_m_dot")
    for i in range(J.shape[0]):
        scale = (np.abs(J[i]) * np.abs(x)).max()
        for j, name in enumerate(s.unknown_names):
            h = 1e-6 * max(abs(x[j]), 1e-3)
            e = np.zeros_like(x)
            e[j] = h
            fd = (F(x + e)[i] - F(x - e)[i]) / (2 * h)
            assert J[i, j] * abs(x[j]) == pytest.approx(
                fd * abs(x[j]), rel=1e-5, abs=1e-7 * scale
            ), (
                i,
                name,
            )
    kinds = s._row_kinds(x)
    assert ("node 'fuel' Tt_target", "T") in kinds
    assert np.any(J[:, j_fuel] != 0.0)


def test_it_must_feed_a_combustor_directly() -> None:
    n = _net(1000.0)
    n.nodes["fuel"].combustor = "out"
    with pytest.raises(ValueError, match="is not a CombustorNode"):
        NetworkSolver(n).solve()
    with pytest.raises(ValueError, match="branch must be"):
        FuelBoundary("f", combustor="c", Tt_target=1000.0, branch="stoich")
