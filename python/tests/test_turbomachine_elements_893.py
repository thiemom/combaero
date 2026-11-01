"""#893: CompressorElement and TurbineElement in a network -- the outlet
state they deliver, the work they report, the ledger and second law that
count it, their Jacobian, and what they refuse."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    CombustorNode,
    CompressorElement,
    FlowNetwork,
    FuelBoundary,
    LosslessConnectionElement,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    TurbineElement,
)
from combaero.network.components import OrificeElement

AIR = list(cb.mole_to_mass(cb.species.dry_air()))
X_AIR = list(cb.species.dry_air())
CH4 = list(cb.species.to_mass(cb.species.pure_species("CH4")))


def _orifice(eid: str, a: str, b: str, d: float = 0.05) -> OrificeElement:
    return OrificeElement(eid, a, b, diameter=d, Cd=0.8, correlation="fixed", in_pipe=False)


def _between(cls, PR: float, Pin: float, Tin: float, Pout: float, **eta) -> FlowNetwork:
    """in -> orifice -> a -> machine -> b -> orifice -> out."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=Pin, Tt=Tin, Y=AIR))
    n.add_node(PlenumNode("a"))
    n.add_node(PlenumNode("b"))
    n.add_node(PressureBoundary("out", Pt=Pout, Tt=Tin, Y=AIR))
    n.add_element(_orifice("o1", "in", "a"))
    n.add_element(cls("m", "a", "b", PR, **eta))
    n.add_element(_orifice("o2", "b", "out"))
    return n


def _solve(net: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return s, s.solve()


CASES = [
    (CompressorElement, 4.0, 1.0e5, 288.15, 3.0e5, cb.compression),
    (TurbineElement, 4.0, 1.0e6, 1400.0, 2.0e5, cb.expansion),
]


@pytest.mark.parametrize("eta", [{"eta_poly": 0.9}, {"eta_is": 0.85}])
@pytest.mark.parametrize(("cls", "PR", "Pin", "Tin", "Pout", "kernel"), CASES)
def test_it_delivers_the_path_outlet_and_reports_its_work(
    cls, PR, Pin, Tin, Pout, kernel, eta
) -> None:
    _, r = _solve(_between(cls, PR, Pin, Tin, Pout, **eta))
    assert r["__success__"], r["__message__"]
    p = kernel(r["a.Tt"], X_AIR, PR, **eta)
    assert r["b.Tt"] == pytest.approx(p.T2, rel=1e-10)
    ratio = r["b.Pt"] / r["a.Pt"]
    assert ratio == pytest.approx(PR if cls is CompressorElement else 1.0 / PR, rel=1e-9)
    assert r["m.dh"] == pytest.approx(p.dh, rel=1e-10)
    assert r["m.power"] == pytest.approx(abs(r["m.m_dot"] * p.dh), rel=1e-10)
    assert r["m.power"] > 0.0
    assert (r["m.dh"] > 0.0) if cls is CompressorElement else (r["m.dh"] < 0.0)
    assert r["m.eta_poly"] == pytest.approx(p.eta_poly)
    assert r["m.eta_is"] == pytest.approx(p.eta_is)


@pytest.mark.parametrize(("cls", "PR", "Pin", "Tin", "Pout", "kernel"), CASES)
def test_the_ledger_counts_the_shaft_work(cls, PR, Pin, Tin, Pout, kernel) -> None:
    _, r = _solve(_between(cls, PR, Pin, Tin, Pout, eta_poly=0.9))
    led = r["__ledger__"]
    assert abs(led["energy_closure_rel"]) < 1e-9
    assert led["W_shaft"] == pytest.approx(r["m.m_dot"] * r["m.dh"], rel=1e-9)


def test_a_machine_into_an_outlet_boundary_leaves_at_its_outlet_enthalpy() -> None:
    """The ledger prices a stream leaving through the machine itself at the
    machine's outlet Tt: priced at its inlet's, the missing work closed it by
    accident."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=1.0e5, Tt=288.15, Y=AIR))
    n.add_node(PlenumNode("a"))
    n.add_node(PressureBoundary("out", Pt=2.5e5, Tt=288.15, Y=AIR))
    n.add_element(_orifice("o1", "in", "a"))
    n.add_element(CompressorElement("m", "a", "out", 3.0, eta_poly=0.88))
    _, r = _solve(n)
    assert r["__success__"], r["__message__"]
    assert abs(r["__ledger__"]["energy_closure_rel"]) < 1e-9
    assert r["__ledger__"]["W_shaft"] > 0.0


def _brayton() -> FlowNetwork:
    """Ambient air -> compressor PR 8 -> burner (CH4 to 1400 K, a FuelBoundary)
    -> turbine PR 4 -> plenum -> nozzle -> ambient."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("amb", Pt=1.01325e5, Tt=288.15, Y=AIR))
    n.add_node(PlenumNode("c_in"))
    n.add_node(CombustorNode("burner"))
    n.add_node(FuelBoundary("fuel", combustor="burner", Tt_target=1400.0, Tt=300.0, Y=CH4))
    n.add_node(PlenumNode("t_out"))
    n.add_node(PressureBoundary("exit", Pt=1.01325e5, Tt=800.0, Y=AIR))
    n.add_element(_orifice("intake", "amb", "c_in", d=0.12))
    n.add_element(CompressorElement("comp", "c_in", "burner", 8.0, eta_poly=0.88))
    n.add_element(LosslessConnectionElement("inj", "fuel", "burner"))
    n.add_element(TurbineElement("turb", "burner", "t_out", 4.0, eta_poly=0.9))
    n.add_element(_orifice("nozzle", "t_out", "exit", d=0.07))
    return n


def test_a_brayton_cycle_closes() -> None:
    """Compressor, burner on a fuel controller, turbine: every Tt is the
    path's, the energy balance closes with both machines' work, and the
    turbine delivers more than the compressor absorbs."""
    _, r = _solve(_brayton())
    assert r["__success__"], r["__message__"]
    assert r["burner.Tt"] == pytest.approx(1400.0, abs=1e-3)
    X_b = list(cb.mass_to_mole(_state_Y(r, "burner")))
    t = cb.expansion(r["burner.Tt"], X_b, 4.0, eta_poly=0.9)
    assert r["t_out.Tt"] == pytest.approx(t.T2, rel=1e-9)
    c = cb.compression(r["c_in.Tt"], X_AIR, 8.0, eta_poly=0.88)
    assert r["comp.Tt_out"] == pytest.approx(c.T2, rel=1e-10)
    led = r["__ledger__"]
    assert abs(led["energy_closure_rel"]) < 1e-9
    assert r["turb.power"] > r["comp.power"]


def _state_Y(r: dict, nid: str) -> list[float]:
    n = len(cb.species.names)
    return [float(r[f"{nid}.Y[{k}]"]) for k in range(n)]


@pytest.mark.parametrize("build", ["compressor", "turbine", "brayton"])
def test_its_rows_are_exact(build: str) -> None:
    """Every Jacobian row against central differences at the root, columns
    scaled by |x|: the outlet Tt chained through the inlet's Tt and
    composition (the burner's products into the turbine)."""
    if build == "brayton":
        net = _brayton()
    else:
        c = CASES[0] if build == "compressor" else CASES[1]
        net = _between(*c[:5], eta_poly=0.9)
    s, r = _solve(net)
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    for i in range(J.shape[0]):
        scale = (np.abs(J[i]) * np.abs(x)).max()
        for j, name in enumerate(s.unknown_names):
            h = 1e-6 * max(abs(x[j]), 1e-3)
            e = np.zeros_like(x)
            e[j] = h
            fd = (F(x + e)[i] - F(x - e)[i]) / (2 * h)
            assert J[i, j] * abs(x[j]) == pytest.approx(
                fd * abs(x[j]), rel=1e-5, abs=1e-7 * scale
            ), (i, name)


def test_a_reversed_flow_is_refused() -> None:
    """The outlet held above what PR 2 can reach from the inlet: the network
    drives the flow backwards through the compressor. Refused, never
    returned."""
    _, r = _solve(_between(CompressorElement, 2.0, 1.0e5, 288.15, 3.0e5, eta_poly=0.9))
    assert not r["__success__"]
    assert r["__converged__"]
    assert r["m.m_dot"] < 0.0
    assert "reversed" in r["__message__"]


def test_between_two_pressure_boundaries_is_refused() -> None:
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=1.0e5, Tt=300.0, Y=AIR))
    n.add_node(PressureBoundary("out", Pt=4.0e5, Tt=300.0, Y=AIR))
    n.add_element(CompressorElement("m", "in", "out", 4.0, eta_poly=0.9))
    with pytest.raises(ValueError, match="nothing sets its flow"):
        NetworkSolver(n).solve()


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        ({"pressure_ratio": 4.0}, "exactly one"),
        ({"pressure_ratio": 4.0, "eta_poly": 0.9, "eta_is": 0.9}, "exactly one"),
        ({"pressure_ratio": 4.0, "eta_poly": 1.1}, "outside"),
        ({"pressure_ratio": 0.8, "eta_poly": 0.9}, "< 1"),
    ],
)
def test_invalid_machines_are_refused(kw: dict, match: str) -> None:
    for cls in (CompressorElement, TurbineElement):
        with pytest.raises(ValueError, match=match):
            cls("m", "a", "b", **kw)


def test_the_second_law_sees_the_work_as_entropy_free() -> None:
    """Shaft work carries no entropy: an ideal machine (eta 1) generates none,
    a real one some. Counted as heat or not at all, a compressor read as
    entropy destroyed (Pt up at constant Tt) and was demoted."""
    sig = {}
    for eta in (1.0, 0.8):
        _, r = _solve(_between(CompressorElement, 4.0, 1.0e5, 288.15, 3.0e5, eta_poly=eta))
        assert r["__success__"], r["__message__"]
        sig[eta] = r["__second_law__"]["m"]["sigma_norm"]
    assert abs(sig[1.0]) < 1e-8
    assert sig[0.8] > 1e-3
