"""GUI wiring tests for the ejector element type.

Covers (1) the EjectorData schema, (2) that graph_builder converts an
`ejector` node into an EjectorElement with correct primary/secondary/outlet
wiring, and (3) an end-to-end runner solve whose element diagnostics
(omega, p_c_star_pa) round-trip through the ElementResult serialization.
Mirrors test_gui_mpce_tee.py.
"""

from __future__ import annotations

import pytest

from combaero.network.components import MultiPortChamberElement
from combaero.network.ejector_element import EjectorElement
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import EjectorData, NetworkGraphSchema

_THROAT = 3.14e-5
_NOZZLE_EXIT = 1.0e-4
_MIXING = 8.0e-4


def _node(nid, ntype, data, x=0, y=0):
    return {"id": nid, "type": ntype, "position": {"x": x, "y": y}, "data": data}


def _edge(eid, src, tgt, src_handle=None, tgt_handle=None):
    e = {"id": eid, "source": src, "target": tgt, "data": {}}
    if src_handle:
        e["sourceHandle"] = src_handle
    if tgt_handle:
        e["targetHandle"] = tgt_handle
    return e


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


def test_ejector_data_defaults():
    d = EjectorData()
    assert d.recovery_efficiency == 1.0
    assert d.throat_area > 0.0
    assert d.nozzle_exit_area > d.throat_area
    assert d.mixing_area > d.nozzle_exit_area
    assert d.initial_guess == {}


def test_ejector_data_ignores_unknown_keys():
    # extra="ignore" so stray GUI fields (e.g. a stale param) do not crash.
    d = EjectorData(throat_area=2e-5, bogus_field=123, label="ej")
    assert d.throat_area == 2e-5
    assert d.label == "ej"
    assert not hasattr(d, "bogus_field")


# ---------------------------------------------------------------------------
# Dispatch / wiring
# ---------------------------------------------------------------------------


def _three_plenum_ejector_schema() -> NetworkGraphSchema:
    """Three plena on the ejector's ports (primary/secondary inlets, outlet).
    Plena connect directly to the port handles, so no auto-MCN is inserted and
    the port node ids are exactly the plenum ids -- easy to assert."""
    return NetworkGraphSchema(
        nodes=[
            _node("p_pri", "plenum", {}),
            _node("p_sec", "plenum", {}, y=200),
            _node("p_out", "plenum", {}, x=400),
            _node(
                "ej1",
                "ejector",
                {
                    "throat_area": _THROAT,
                    "nozzle_exit_area": _NOZZLE_EXIT,
                    "mixing_area": _MIXING,
                },
                x=200,
            ),
        ],
        edges=[
            _edge("e_pri", "p_pri", "ej1", tgt_handle="port-primary-target"),
            _edge("e_sec", "p_sec", "ej1", tgt_handle="port-secondary-target"),
            _edge("e_out", "ej1", "p_out", src_handle="port-outlet-source"),
        ],
    )


def test_ejector_dispatch_builds_element_with_correct_ports():
    net = build_network_from_schema(_three_plenum_ejector_schema())
    e = net.elements["ej1"]
    assert isinstance(e, EjectorElement)
    assert isinstance(e, MultiPortChamberElement)
    assert e.primary_node == "p_pri"
    assert e.secondary_node == "p_sec"
    assert e.outlet_node == "p_out"
    assert e.inlet_nodes == ["p_pri", "p_sec"]
    assert e.outlet_nodes == ["p_out"]
    assert e.throat_area == _THROAT
    assert e.nozzle_exit_area == _NOZZLE_EXIT
    assert e.mixing_area == _MIXING
    assert e.recovery_efficiency == 1.0


def test_ejector_resolves_with_bare_boundary_fed_port_mcns():
    """The ejector supplies its own port_areas from geometry, so a network
    whose ports are explicit MomentumChamberNodes fed straight from
    boundaries (via the GUI's auto lossless connection -- no ChannelElement
    with a diameter to infer an area from) still resolves. Regression for the
    "could not infer area at port" failure."""
    schema = NetworkGraphSchema(
        nodes=[
            _node("mb_p", "mass_boundary", {"m_dot": 0.04, "Tt": 440.0}),
            _node("pb_s", "pressure_boundary", {"Pt": 22850.0, "Tt": 280.0}, y=200),
            _node("pb_o", "pressure_boundary", {"Pt": 40000.0, "Tt": 300.0}, x=400),
            _node("mc_p", "momentum_chamber", {"area": 0.1}, x=100),
            _node("mc_s", "momentum_chamber", {"area": 0.1}, x=100, y=200),
            _node("mc_o", "momentum_chamber", {"area": 0.15}, x=300),
            _node(
                "ej1",
                "ejector",
                {
                    "throat_area": _THROAT,
                    "nozzle_exit_area": _NOZZLE_EXIT,
                    "mixing_area": _MIXING,
                },
                x=200,
            ),
        ],
        edges=[
            _edge("e1", "mb_p", "mc_p"),
            _edge("e2", "mc_p", "ej1", tgt_handle="port-primary-target"),
            _edge("e3", "pb_s", "mc_s"),
            _edge("e4", "mc_s", "ej1", tgt_handle="port-secondary-target"),
            _edge("e5", "ej1", "mc_o", src_handle="port-outlet-source"),
            _edge("e6", "mc_o", "pb_o"),
        ],
    )
    net = build_network_from_schema(schema)
    net.resolve_all_topology()  # would raise "could not infer area" before the fix
    e = net.elements["ej1"]
    assert e.port_areas == [_NOZZLE_EXIT, _MIXING - _NOZZLE_EXIT, _MIXING]


def test_ejector_dispatch_auto_inserts_port_momentum_chambers():
    """When a channel (element) connects to an ejector port, an implicit
    MomentumChamberNode is auto-inserted at that port -- the same behaviour
    the tee relies on, generalized to the ejector's junction_ids."""
    schema = NetworkGraphSchema(
        nodes=[
            _node("pb_p", "pressure_boundary", {"Pt": 7e5, "Tt": 440.0}),
            _node("p_sec", "plenum", {}, y=200),
            _node("p_out", "plenum", {}, x=400),
            _node("ch_p", "channel", {"L": 0.2, "D": 0.05}, x=100),
            _node(
                "ej1",
                "ejector",
                {
                    "throat_area": _THROAT,
                    "nozzle_exit_area": _NOZZLE_EXIT,
                    "mixing_area": _MIXING,
                },
                x=200,
            ),
        ],
        edges=[
            _edge("e1", "pb_p", "ch_p"),
            _edge("e2", "ch_p", "ej1", tgt_handle="port-primary-target"),
            _edge("e_sec", "p_sec", "ej1", tgt_handle="port-secondary-target"),
            _edge("e_out", "ej1", "p_out", src_handle="port-outlet-source"),
        ],
    )
    net = build_network_from_schema(schema)
    e = net.elements["ej1"]
    # primary port resolves to the auto-inserted MomentumChamberNode, not the channel.
    assert e.primary_node == "__tee_jct__ej1_primary"
    assert e.primary_node in net.nodes


# ---------------------------------------------------------------------------
# End-to-end runner solve: diagnostics serialize through ElementResult.
# ---------------------------------------------------------------------------


def _solvable_ejector_schema() -> dict:
    """Cold reference network (boundary -> port MCN -> ejector), the Huang
    T3-EH operating point, with NO manual initial guesses. graph_builder's
    ejector warm-start seeds the port pressures + P_jct so the stiff
    double-choked solve converges from cold -- which it does NOT without the
    seed (confirmed: no init_strategy cracks it cold)."""
    return {
        "nodes": [
            _node("pb_p", "pressure_boundary", {"Pt": 7e5, "Tt": 440.0}),
            _node("pb_s", "pressure_boundary", {"Pt": 22850.0, "Tt": 280.0}, y=200),
            _node("pb_o", "pressure_boundary", {"Pt": 40000.0, "Tt": 300.0}, x=600),
            _node("mc_p", "momentum_chamber", {"area": 0.05}, x=150),
            _node("mc_s", "momentum_chamber", {"area": 0.05}, x=150, y=200),
            _node("mc_o", "momentum_chamber", {"area": 0.05}, x=450),
            _node(
                "ej1",
                "ejector",
                {
                    "throat_area": _THROAT,
                    "nozzle_exit_area": _NOZZLE_EXIT,
                    "mixing_area": _MIXING,
                },
                x=300,
            ),
        ],
        "edges": [
            _edge("e1", "pb_p", "mc_p"),
            _edge("e2", "mc_p", "ej1", tgt_handle="port-primary-target"),
            _edge("e3", "pb_s", "mc_s"),
            _edge("e4", "mc_s", "ej1", tgt_handle="port-secondary-target"),
            _edge("e5", "ej1", "mc_o", src_handle="port-outlet-source"),
            _edge("e6", "mc_o", "pb_o"),
        ],
    }


def test_ejector_warmstart_seeds_port_pressures():
    """graph_builder should seed the port MCN pressures from the feeding
    boundaries (primary=700k, secondary=22.85k, outlet=40k) when the user
    provided none."""
    net = build_network_from_schema(NetworkGraphSchema(**_solvable_ejector_schema()))
    assert net.nodes["mc_p"].initial_guess.get("mc_p.Pt") == 7e5
    assert net.nodes["mc_s"].initial_guess.get("mc_s.Pt") == 22850.0
    assert net.nodes["mc_o"].initial_guess.get("mc_o.Pt") == 40000.0
    # P_jct (== P_c*) seeded on the element from the ejector's own physics.
    assert net.elements["ej1"].initial_guess.get("ej1.P_jct", 0.0) > 0.0


def test_ejector_runner_solve_converges_cold_and_exposes_diagnostics():
    runner = NetworkRunner.from_dict(_solvable_ejector_schema())
    result = runner.solve(timeout=120.0)
    df = result.to_dataframe()
    row = df[df["id"] == "ej1"]
    assert not row.empty, "ejector element row missing from results"
    assert "omega" in df.columns and "p_c_star_pa" in df.columns
    assert bool(row["success"].iloc[0]), "cold ejector solve did not converge"
    omega = float(row["omega"].iloc[0])
    p_c = float(row["p_c_star_pa"].iloc[0])
    assert 0.0 < omega < 3.0
    assert p_c > 0.0


def _jetpump_ejector_schema() -> dict:
    """Low primary mass flow (0.0051 kg/s) into the same geometry against a
    101.325 kPa suction/outlet. That mass flow is well below the primary's
    choke threshold at this back pressure, so the physical operating point is
    the UNCHOKED subsonic jet-pump regime (primary Pt floats to ~the back
    pressure, dp >= 0, omega ~ 7) -- the reported case that previously failed
    to converge (converged instead to a spurious double-choked root with
    Pt_primary below the back pressure). See OPERATING_REGIMES_DESIGN.md."""
    return {
        "nodes": [
            _node("mb_p", "mass_boundary", {"m_dot": 0.0051, "Tt": 300.0}),
            _node("pb_s", "pressure_boundary", {"Pt": 101325.0, "Tt": 300.0}, y=200),
            _node("pb_o", "pressure_boundary", {"Pt": 101325.0, "Tt": 300.0}, x=600),
            _node("mc_p", "momentum_chamber", {"area": 0.1}, x=150),
            _node("mc_s", "momentum_chamber", {"area": 0.1}, x=150, y=200),
            _node("mc_o", "momentum_chamber", {"area": 0.15}, x=450),
            _node(
                "ej1",
                "ejector",
                {
                    "throat_area": _THROAT,
                    "nozzle_exit_area": _NOZZLE_EXIT,
                    "mixing_area": _MIXING,
                },
                x=300,
            ),
        ],
        "edges": [
            _edge("e1", "mb_p", "mc_p"),
            _edge("e2", "mc_p", "ej1", tgt_handle="port-primary-target"),
            _edge("e3", "pb_s", "mc_s"),
            _edge("e4", "mc_s", "ej1", tgt_handle="port-secondary-target"),
            _edge("e5", "ej1", "mc_o", src_handle="port-outlet-source"),
            _edge("e6", "mc_o", "pb_o"),
        ],
    }


def test_ejector_runner_solves_unchoked_jet_pump_regime():
    """The reported low-flow case must converge to the physical jet-pump root:
    unchoked primary (s_choke -> 0), forward flow (dp >= 0), and a large
    entrainment ratio -- NOT the spurious double-choked artifact."""
    runner = NetworkRunner.from_dict(_jetpump_ejector_schema())
    result = runner.solve(timeout=120.0)
    df = result.to_dataframe()
    row = df[df["id"] == "ej1"]
    assert not row.empty, "ejector element row missing from results"
    assert bool(row["success"].iloc[0]), "jet-pump ejector solve did not converge"

    mp = float(row["m_dot_primary"].iloc[0])
    ms = float(row["m_dot_secondary"].iloc[0])
    assert mp == pytest.approx(0.0051, rel=1e-3)
    assert ms > 0.0  # entrained flow is drawn in (positive suction)
    omega = float(row["omega"].iloc[0])
    assert 5.0 < omega < 9.0, f"expected jet-pump omega ~ 7, got {omega}"
    # Unchoked regime, not the double-choked artifact.
    assert float(row["s_choke"].iloc[0]) == pytest.approx(0.0, abs=1e-6)
    assert float(row["critical_mode"].iloc[0]) == 0.0
    # Primary floats up to ~the back pressure: forward flow, dp >= 0.
    mc_p = df[df["id"] == "mc_p"]
    if not mc_p.empty and "Pt" in df.columns:
        assert float(mc_p["Pt"].iloc[0]) >= 101325.0 * 0.999


def test_ejector_element_m_dot_reports_total_throughflow():
    """The ejector has no own `.m_dot` unknown (a junction). Its ElementResult
    m_dot -- shown on the node card and telemetry -- must report the total
    outlet flow (mp + ms), not 0.0 (the old fallback only knew the tee's
    `m_dot_com`)."""
    result = NetworkRunner.from_dict(_solvable_ejector_schema()).solve(timeout=120.0)
    er = result._element_results["ej1"]
    d = er.model_dump()
    assert er.m_dot > 0.0
    assert er.m_dot == pytest.approx(d["m_dot_outlet"], rel=1e-9)
    assert er.m_dot == pytest.approx(d["m_dot_primary"] + d["m_dot_secondary"], rel=1e-6)
