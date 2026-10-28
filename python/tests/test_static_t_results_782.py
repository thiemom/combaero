"""``T`` is the static temperature in results and on the state (#782, 0.8).

``<node>.T`` and an element's ``T_in``/``T_out`` held the stagnation T
until 0.7; they are static now, as ``P``/``P_in`` already were, and the
stagnation value is ``<node>.Tt``/``Tt_in``. The static T is regime-aware:
an incompressible node's is its ``Tt`` (the low-Mach limit its closure
solves), a compressible one's lies on the isentrope through (Tt, Pt) at P.
"""

from __future__ import annotations

import warnings

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkMixtureState,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)
from combaero.network import solver as solver_module

Y_AIR = list(cb.species.dry_air_mass())
X_AIR = list(cb.species.dry_air())


def _chamber_network(regime: str) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("inlet", m_dot=0.5, Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("chamber", area=3e-3, regime=regime))
    g.add_node(PressureBoundary("outlet", Pt=101325.0))
    for eid, a, b in (("orf1", "inlet", "chamber"), ("orf2", "chamber", "outlet")):
        g.add_element(
            OrificeElement(
                eid, a, b, Cd=0.6, diameter=0.04, correlation="fixed", regime="incompressible"
            )
        )
    return g


@pytest.fixture
def fresh_warning(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(solver_module, "_STATIC_T_WARNED", False)


def test_a_moving_compressible_node_reports_a_static_T_below_its_Tt(fresh_warning: None) -> None:
    sol = NetworkSolver(_chamber_network("compressible")).solve()
    assert sol["__success__"]
    with pytest.warns(FutureWarning, match="STATIC temperature since combaero 0.8"):
        T = sol["chamber.T"]
    Tt, P, Pt = sol["chamber.Tt"], sol["chamber.P"], sol["chamber.Pt"]
    assert Tt - 1.0 > T
    # On the isentrope through (Tt, Pt): the solved closure's static state.
    assert cb.s(T, X_AIR, P) == pytest.approx(cb.s(Tt, X_AIR, Pt), rel=1e-8)
    # The element leaving the chamber reports the chamber's static and total T.
    diag = sol["__element_diag__"]["orf2"]
    assert diag["T_in"] == pytest.approx(T, rel=1e-12)  # isentrope vs closure: round-off
    assert diag["Tt_in"] == Tt


def test_the_warning_is_shown_once(fresh_warning: None) -> None:
    sol = NetworkSolver(_chamber_network("compressible")).solve()
    with pytest.warns(FutureWarning):
        sol["chamber.T"]
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        sol["chamber.T"]
        sol["__element_diag__"]["orf2"]["T_in"]


def test_the_old_T_static_key_is_the_new_T_with_a_deprecation(fresh_warning: None) -> None:
    sol = NetworkSolver(_chamber_network("compressible")).solve()
    with pytest.warns(DeprecationWarning, match=r"result\['chamber.T'\]"):
        T_static = sol["chamber.T_static"]
    assert T_static == dict.__getitem__(sol, "chamber.T")


def test_an_incompressible_node_and_a_plenum_report_T_equal_to_Tt(fresh_warning: None) -> None:
    sol = NetworkSolver(_chamber_network("incompressible")).solve()
    assert sol["__success__"]
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        assert sol["chamber.T"] == sol["chamber.Tt"]
        assert sol["inlet.T"] == sol["inlet.Tt"]
    assert sol["chamber.P"] < sol["chamber.Pt"]  # it moves; the regime decides


def test_the_state_computes_its_static_T() -> None:
    plenum = NetworkMixtureState(P=2.0e5, Pt=2.0e5, Tt=600.0, Y=Y_AIR)
    assert plenum.T == 600.0
    moving = NetworkMixtureState(P=1.8e5, Pt=2.0e5, Tt=600.0, Y=Y_AIR)
    assert moving.T < 600.0
    assert cb.s(moving.T, X_AIR, 1.8e5) == pytest.approx(cb.s(600.0, X_AIR, 2.0e5), rel=1e-10)
    low_mach = NetworkMixtureState(P=1.8e5, Pt=2.0e5, Tt=600.0, Y=Y_AIR, incompressible=True)
    assert low_mach.T == 600.0


def test_the_0_7_constructor_and_setter_still_work_and_say_so() -> None:
    with pytest.warns(DeprecationWarning, match="0.9.0"):
        old = NetworkMixtureState(2.0e5, 2.0e5, 600.0, 600.0, 0.1, Y_AIR)
    assert (old.P, old.Pt, old.Tt, old.m_dot) == (2.0e5, 2.0e5, 600.0, 0.1)
    with pytest.warns(DeprecationWarning, match="0.9.0"):
        only_T = NetworkMixtureState(P=2.0e5, Pt=2.0e5, T=650.0, Y=Y_AIR)
    assert only_T.Tt == 650.0
    with pytest.warns(DeprecationWarning, match="0.9.0"):
        old.T = 610.0
    assert old.Tt == 610.0
    with pytest.raises(TypeError, match="Tt"):
        NetworkMixtureState(P=2.0e5, Pt=2.0e5, Y=Y_AIR)
