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
    PlenumNode,
    PressureBoundary,
)
from combaero.network.components import ChannelElement, OrificeElement

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
    assert "reaches with no injected fuel, 389.2 K" in r["__message__"]


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


# -- every species as the fuel (user sweep, 2026-11-01) ----------------------

_FUELS = [n for n in cb.species.names if n not in ("N2", "O2", "AR", "CO2", "H2O")]


def _adiabatic(Y_fuel: list[float], target: float, branch: str, Y_ox: list[float] = AIR):
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("air", m_dot=0.3, Tt=T_AIR, Y=Y_ox))
    n.add_node(FuelBoundary("fuel", combustor="burner", Tt_target=target, branch=branch, Y=Y_fuel))
    n.add_node(CombustorNode("burner"))
    n.add_node(PressureBoundary("out", Pt=2.0e5, Tt=target, Y=AIR))
    n.add_element(LosslessConnectionElement("a_in", "air", "burner"))
    n.add_element(LosslessConnectionElement("inj", "fuel", "burner"))
    n.add_element(
        OrificeElement(
            "exit", "burner", "out", diameter=0.05, Cd=0.8, correlation="fixed", in_pipe=False
        )
    )
    return _solve(n)[1]


def _burnt(Y_ox: list[float], Y_fuel: list[float], m_fuel: float) -> float:
    def stream(T: float, Y: list[float], m: float) -> cb.Stream:
        st = cb.Stream()
        st.set_T(T)
        st.set_P(2.0e5)
        st.set_X(list(cb.mass_to_mole(Y)))
        st.set_mdot(m)
        return st

    mixed = cb.mix([stream(T_AIR, Y_ox, 0.3), stream(300.0, Y_fuel, m_fuel)])
    return float(cb.complete_combustion(mixed.T, list(mixed.X), mixed.P).T)


@pytest.mark.parametrize("name", _FUELS)
@pytest.mark.parametrize(("target", "branch"), [(1200.0, "lean"), (1500.0, "rich")])
def test_every_fuel_species_meets_the_target(name: str, target: float, branch: str) -> None:
    """Hydrocarbons, CO, H2 and NH3, both branches: the target, on its side of
    stoichiometric, and the combustor's own model burns that fuel to it."""
    Y = list(cb.species.to_mass(cb.species.pure_species(name)))
    r = _adiabatic(Y, target, branch)
    assert r["__success__"], (name, r["__message__"])
    assert r["burner.Tt"] == pytest.approx(target, abs=1e-3)
    assert (r["fuel.phi"] < 1.0) if branch == "lean" else (r["fuel.phi"] > 1.0)
    assert _burnt(AIR, Y, r["fuel.m_dot"]) == pytest.approx(target, abs=1e-3)


@pytest.mark.parametrize(
    ("name", "what"),
    [
        ("N2", "nothing combustible"),
        ("AR", "nothing combustible"),
        ("CO2", "nothing combustible"),
        ("H2O", "nothing combustible"),
        ("O2", "an oxidizer, not a fuel"),
    ],
)
def test_a_non_fuel_is_refused_by_name(name: str, what: str) -> None:
    Y = list(cb.species.to_mass(cb.species.pure_species(name)))
    with pytest.raises(ValueError, match=f"\\({name}\\) is {what}"):
        _adiabatic(Y, 1200.0, "lean")


def test_phi_does_not_count_water_vapour_as_fuel() -> None:
    """Humid air carries H in its water; the reported phi is the reactants'
    combustibles over their O2, so it starts from 0, not from 0.026."""
    Y = list(cb.species.to_mass(cb.species.pure_species("CH4")))
    r = _adiabatic(Y, 1200.0, "lean")
    names = list(cb.species.names)
    Y_u = [
        (0.3 * a + r["fuel.m_dot"] * f) / (0.3 + r["fuel.m_dot"])
        for a, f in zip(AIR, Y, strict=True)
    ]
    need = cb.oxygen_required_per_kg_mixture(list(cb.mass_to_mole(Y_u)))
    assert r["fuel.phi"] == pytest.approx(need / Y_u[names.index("O2")], rel=1e-9)


def _premixed(x_ch4: float) -> list[float]:
    X = [v * (1.0 - x_ch4) for v in cb.mass_to_mole(AIR)]
    X[list(cb.species.names).index("CH4")] += x_ch4
    return list(cb.mole_to_mass(X))


def test_air_carrying_fuel_still_meets_a_target_above_its_own_burn() -> None:
    """3% CH4 in the air burns alone to ~1130 K: a 1430 K target needs only
    the rest, and the model burns air + injected fuel to it."""
    CH4 = list(cb.species.to_mass(cb.species.pure_species("CH4")))
    Y_ox = _premixed(0.03)
    r = _adiabatic(CH4, 1430.0, "lean", Y_ox)
    assert r["__success__"], r["__message__"]
    assert r["fuel.Tt_without_fuel"] == pytest.approx(1130.0, abs=1.0)
    assert _burnt(Y_ox, CH4, r["fuel.m_dot"]) == pytest.approx(1430.0, abs=1e-3)


def test_air_carrying_fuel_refuses_a_target_below_its_own_burn() -> None:
    CH4 = list(cb.species.to_mass(cb.species.pure_species("CH4")))
    r = _adiabatic(CH4, 1030.0, "lean", _premixed(0.03))
    assert not r["__success__"]
    assert "reaches with no injected fuel" in r["__message__"]
    assert "the oxidizer carries fuel of its own" in r["__message__"]


# -- edge cases (2026-11-01) --------------------------------------------------


def _with_fuel(n: FlowNetwork, fid: str, target: float = 1200.0, **kw) -> FlowNetwork:
    n.add_node(FuelBoundary(fid, combustor="burner", Tt_target=target, Y=CH4, **kw))
    n.add_element(LosslessConnectionElement(f"{fid}_inj", fid, "burner"))
    return n


def _adiabatic_base(m_air: float = 0.3, method: str = "complete") -> FlowNetwork:
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("air", m_dot=m_air, Tt=T_AIR, Y=AIR))
    n.add_node(CombustorNode("burner", method=method))
    n.add_node(PressureBoundary("out", Pt=2.0e5, Tt=1200.0, Y=AIR))
    n.add_element(LosslessConnectionElement("a_in", "air", "burner"))
    n.add_element(
        OrificeElement(
            "exit", "burner", "out", diameter=0.05, Cd=0.8, correlation="fixed", in_pipe=False
        )
    )
    return n


def test_two_controllers_on_one_combustor_are_refused() -> None:
    """Two targets for one Tt would leave the pilot/main split free: a
    'converged' answer with an arbitrary ratio."""
    n = _with_fuel(_with_fuel(_adiabatic_base(), "fuel"), "pilot")
    with pytest.raises(ValueError, match="already has a fuel controller"):
        NetworkSolver(n).solve()


def test_a_controller_and_a_tt_target_heat_boundary_are_refused() -> None:
    from combaero.network import EnergyBoundary

    n = _with_fuel(_adiabatic_base(), "fuel")
    n.nodes["burner"].add_energy_boundary(EnergyBoundary("q", Tt_target=1150.0))
    with pytest.raises(ValueError, match="also has a Tt_target heat boundary"):
        NetworkSolver(n).solve()


@pytest.mark.parametrize("above", [0.01, 1.0])
def test_a_target_just_above_the_no_fuel_tt_is_met(above: float) -> None:
    """Micrograms per second of fuel: the start is found by the combustor's own
    model, not a 1 K search that returns no fuel."""
    _, r = _solve(_with_fuel(_adiabatic_base(), "fuel", target=T_AIR + above))
    assert r["__success__"], r["__message__"]
    assert r["burner.Tt"] == pytest.approx(T_AIR + above, abs=1e-4)
    assert 0.0 < r["fuel.m_dot"] < 1e-4


@pytest.mark.parametrize(
    ("target", "branch"), [(1200.0, "lean"), (1500.0, "rich"), (2000.0, "rich")]
)
def test_an_equilibrium_combustor_is_controlled_too(target: float, branch: str) -> None:
    _, r = _solve(
        _with_fuel(_adiabatic_base(method="equilibrium"), "fuel", target=target, branch=branch)
    )
    assert r["__success__"], r["__message__"]
    assert r["burner.Tt"] == pytest.approx(target, abs=1e-3)
    assert (r["fuel.phi"] < 1.0) if branch == "lean" else (r["fuel.phi"] > 1.0)


@pytest.mark.parametrize("guess", [None, 1e-6, 1e-10, 1e-14])
def test_no_oxidizer_flow_is_named(guess: float | None) -> None:
    """With no air every flow can vanish and the Tt row is then met by an
    undetermined state: from a 1e-10 start (and on Linux CI from the default
    one) the solve reaches tolerance. Refused at any start, never returned."""
    _, r = _solve(_with_fuel(_adiabatic_base(m_air=0.0), "fuel", m_dot_guess=guess))
    assert not r["__success__"]
    assert "gets no oxidizer flow" in r["__message__"]


# -- fuel routes: the boundary names its combustor, like a probe -------------


def _routed(route: str, target: float = 1400.0) -> FlowNetwork:
    """Adiabatic, fixed air flow, the fuel reaching the burner by `route`:
    "direct" (one lossless injector), "premix" (a plenum the air also enters,
    then the burner), "line" (a friction fuel line into the burner)."""
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("air", m_dot=0.3, Tt=T_AIR, Y=AIR))
    n.add_node(FuelBoundary("fuel", combustor="burner", Tt_target=target, Tt=300.0, Y=CH4))
    n.add_node(CombustorNode("burner"))
    n.add_node(PressureBoundary("out", Pt=2.0e5, Tt=target, Y=AIR))
    if route == "premix":
        n.add_node(PlenumNode("premix"))
        n.add_element(LosslessConnectionElement("a_in", "air", "premix"))
        n.add_element(LosslessConnectionElement("inj", "fuel", "premix"))
        n.add_element(LosslessConnectionElement("mix_out", "premix", "burner"))
    else:
        n.add_element(LosslessConnectionElement("a_in", "air", "burner"))
        if route == "line":
            n.add_node(PlenumNode("manifold"))
            n.add_element(ChannelElement("line", "fuel", "manifold", length=0.5, diameter=0.004))
            n.add_element(LosslessConnectionElement("inj", "manifold", "burner"))
        else:
            n.add_element(LosslessConnectionElement("inj", "fuel", "burner"))
    n.add_element(
        OrificeElement(
            "exit", "burner", "out", diameter=0.05, Cd=0.8, correlation="fixed", in_pipe=False
        )
    )
    return n


@pytest.mark.parametrize("route", ["premix", "line"])
@pytest.mark.parametrize(("target", "branch"), [(1400.0, "lean"), (2000.0, "rich")])
def test_the_fuel_route_does_not_change_the_answer(route: str, target: float, branch: str) -> None:
    """Adiabatic at a fixed air flow the fuel flow is set by the energy
    balance alone: through a premix plenum or a lossy fuel line it is the
    direct-injection answer."""
    ref_net = _routed("direct", target)
    ref_net.nodes["fuel"].branch = branch
    _, ref = _solve(ref_net)
    assert ref["__success__"], ref["__message__"]
    net = _routed(route, target)
    net.nodes["fuel"].branch = branch
    _, r = _solve(net)
    assert r["__success__"], r["__message__"]
    assert r["burner.Tt"] == pytest.approx(target, abs=1e-3)
    assert r["fuel.m_dot"] == pytest.approx(ref["fuel.m_dot"], rel=1e-6)
    assert r["fuel.phi"] == pytest.approx(ref["fuel.phi"], rel=1e-6)
    assert r["fuel.Tt_without_fuel"] == pytest.approx(ref["fuel.Tt_without_fuel"], abs=1e-6)


@pytest.mark.parametrize(
    ("target", "says"),
    [(350.0, "reaches with no injected fuel, 389.2 K"), (3000.0, "is above the peak")],
)
def test_a_premix_plenum_names_an_unreachable_target(target: float, says: str) -> None:
    """Through a premix plenum the failed point need not close its mass rows;
    the reason still names the network's oxidizer, not 'no oxidizer'."""
    _, r = _solve(_routed("premix", target))
    assert not r["__success__"]
    assert says in r["__message__"]
    assert "no oxidizer" not in r["__message__"]


def test_fuel_that_can_bypass_the_combustor_is_refused() -> None:
    """A split that sends part of the fuel past the burner leaves the fuel
    the target controls undetermined by the boundary alone: refused."""
    n = _routed("premix")
    n.add_element(LosslessConnectionElement("leak", "premix", "out"))
    with pytest.raises(ValueError, match="part of its fuel can reach"):
        NetworkSolver(n).solve()


def test_fuel_with_no_path_to_its_combustor_is_refused() -> None:
    n = _routed("direct")
    n.elements.pop("inj")
    n.add_node(PlenumNode("dump"))
    n.add_element(LosslessConnectionElement("inj", "fuel", "dump"))
    n.add_element(LosslessConnectionElement("dump_out", "dump", "out"))
    with pytest.raises(ValueError, match="fuel"):
        NetworkSolver(n).solve()
