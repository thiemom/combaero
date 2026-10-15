"""GUI wiring tests for the ejector element type.

Covers (1) the EjectorData schema, (2) that graph_builder converts an
`ejector` node into an EjectorElement with correct primary/secondary/outlet
wiring, and (3) an end-to-end runner solve whose element diagnostics
(omega, p_c_star_pa) round-trip through the ElementResult serialization.
Mirrors test_gui_mpce_tee.py.
"""

from __future__ import annotations

import pytest

from combaero.network.components import MomentumChamberNode, MultiPortChamberBase
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
    # None = inherit A_3 from the outlet pipe, else the element default (#559).
    assert d.mixing_area is None
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
    """Three plena wired straight to the ejector's ports (primary/secondary
    inlets, outlet). Each port gets an auto-inserted MomentumChamberNode joined
    to its plenum by a lossless link (a plenum cannot be the port node)."""
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
    assert isinstance(e, MultiPortChamberBase)
    assert e.primary_node == "__tee_jct__ej1_primary"
    assert e.secondary_node == "__tee_jct__ej1_secondary"
    assert e.outlet_node == "__tee_jct__ej1_outlet"
    for port in ("primary", "secondary", "outlet"):
        assert isinstance(net.nodes[f"__tee_jct__ej1_{port}"], MomentumChamberNode)
    # Flow follows the canvas edge: plenum -> port on inlets, port -> plenum out.
    assert "__auto_link__p_pri____tee_jct__ej1_primary" in net.elements
    assert "__auto_link__p_sec____tee_jct__ej1_secondary" in net.elements
    assert "__auto_link____tee_jct__ej1_outlet__p_out" in net.elements
    net.resolve_all_topology()  # a plenum on a port used to be refused here
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


# ---------------------------------------------------------------------------
# Mixing-area (A_3) inheritance with the outlet pipe, both directions (#559).
# ---------------------------------------------------------------------------


def _ejector_with_outlet_pipe(base: dict, mixing: float | None, pipe_d: float | None) -> dict:
    """Replace the outlet MCN of a reference schema by a channel ch_o."""
    g = base
    next(n for n in g["nodes"] if n["id"] == "ej1")["data"]["mixing_area"] = mixing
    g["nodes"] = [n for n in g["nodes"] if n["id"] != "mc_o"]
    g["nodes"].append(_node("ch_o", "channel", {"L": 0.1, "D": pipe_d}, x=450))
    g["edges"] = [e for e in g["edges"] if e["id"] not in ("e5", "e6")]
    g["edges"] += [
        _edge("e5", "ej1", "ch_o", src_handle="port-outlet-source"),
        _edge("e6", "ch_o", "pb_o"),
    ]
    return g


def _preview(g: dict):
    from gui.backend.resolve_preview import resolve_preview

    return resolve_preview(build_network_from_schema(NetworkGraphSchema(**g)))


def test_outlet_pipe_inherits_the_ejector_mixing_area():
    rep = _preview(_ejector_with_outlet_pipe(_solvable_ejector_schema(), _MIXING, None))
    d = rep["ch_o"].fields["diameter"]
    assert rep["ch_o"].status == "resolved"
    assert d.value == pytest.approx((4 * _MIXING / 3.141592653589793) ** 0.5, rel=1e-12)
    assert d.source == "EjectorElement 'ej1'"
    assert rep["ej1"].fields["mixing_area"].source == "user"


def test_ejector_inherits_its_mixing_area_from_the_outlet_pipe():
    pipe_d = 0.04
    g = _ejector_with_outlet_pipe(_solvable_ejector_schema(), None, pipe_d)
    net = build_network_from_schema(NetworkGraphSchema(**g))
    from gui.backend.resolve_preview import resolve_preview

    m = resolve_preview(net)["ej1"].fields["mixing_area"]
    a3 = 3.141592653589793 / 4 * pipe_d**2
    assert m.value == pytest.approx(a3, rel=1e-12)
    assert m.source == "channel 'ch_o'"
    e = net.elements["ej1"]
    # All three port faces follow: nozzle exit, suction annulus, mixing duct.
    assert e.port_areas == pytest.approx([_NOZZLE_EXIT, a3 - _NOZZLE_EXIT, a3])


def test_ejector_without_outlet_pipe_uses_the_default_mixing_area():
    from combaero.network.ejector_element import DEFAULT_EJECTOR_MIXING_AREA

    g = _solvable_ejector_schema()
    next(n for n in g["nodes"] if n["id"] == "ej1")["data"]["mixing_area"] = None
    m = _preview(g)["ej1"].fields["mixing_area"]
    assert m.value == DEFAULT_EJECTOR_MIXING_AREA
    assert m.source == "default"


def test_outlet_pipe_narrower_than_the_nozzle_exit_is_refused():
    rep = _preview(_ejector_with_outlet_pipe(_solvable_ejector_schema(), None, 0.008))
    assert rep["ej1"].status == "unresolved"
    assert "nozzle_exit_area" in rep["ej1"].reason


def test_warmstart_seeds_an_outlet_port_behind_a_pipe():
    """A pipe on the outlet makes the port an auto-inserted MCN with no canvas
    edges; the warm-start must still trace to the outlet boundary."""
    g = _ejector_with_outlet_pipe(_solvable_ejector_schema(), _MIXING, None)
    net = build_network_from_schema(NetworkGraphSchema(**g))
    port = "__tee_jct__ej1_outlet"
    assert net.nodes[port].initial_guess.get(f"{port}.Pt") == 40000.0


def test_jet_pump_with_an_outlet_pipe_converges():
    """Regression: with the outlet seed missing, the jet pump behind an outlet
    pipe stalled (it converged only with a bare MCN outlet)."""
    g = _ejector_with_outlet_pipe(_jetpump_ejector_schema(), _MIXING, None)
    df = NetworkRunner.from_dict(g).solve(timeout=120.0).to_dataframe()
    row = df[df["id"] == "ej1"]
    assert bool(row["success"].iloc[0])
    assert float(row["omega"].iloc[0]) > 3.0


def test_cleared_throat_and_nozzle_exit_mean_the_default():
    """The Inspector shows an unset throat / nozzle exit as its default; a
    cleared field arrives as null and must build with that default, not fail."""
    d = EjectorData(throat_area=None, nozzle_exit_area=None)
    assert d.throat_area == EjectorData.model_fields["throat_area"].default
    assert d.nozzle_exit_area == EjectorData.model_fields["nozzle_exit_area"].default


def test_boundary_wired_straight_to_a_port_resolves():
    """A pressure boundary on the suction port (no channel, no chamber) gets a
    port node + lossless link instead of failing 'must be a
    MomentumChamberNode'; the warm-start still traces to the boundary."""
    g = _solvable_ejector_schema()
    g["nodes"] = [n for n in g["nodes"] if n["id"] != "mc_s"]
    g["edges"] = [e for e in g["edges"] if e["id"] not in ("e3", "e4")]
    g["edges"].append(_edge("e4", "pb_s", "ej1", tgt_handle="port-secondary-target"))
    net = build_network_from_schema(NetworkGraphSchema(**g))
    net.resolve_all_topology()
    port = "__tee_jct__ej1_secondary"
    assert net.elements["ej1"].secondary_node == port
    assert net.nodes[port].initial_guess.get(f"{port}.Pt") == 22850.0


# ---------------------------------------------------------------------------
# Double-choked ejector behind an outlet pipe (#567): port-chain flow seeding.
# ---------------------------------------------------------------------------

_A3_PIPE_D = (4 * _MIXING / 3.141592653589793) ** 0.5  # pipe with area A_3


def _critical_omega_bare_outlet() -> float:
    df = NetworkRunner.from_dict(_solvable_ejector_schema()).solve(timeout=120.0).to_dataframe()
    return float(df[df["id"] == "ej1"]["omega"].iloc[0])


def test_warmstart_seeds_port_chain_flows_in_the_critical_regime():
    """The cold start used the generic flow defaults (28 kg/s on the outlet
    pipe for a ~0.07 kg/s ejector). The warm-start now seeds each port chain
    from the closed-form operating point: mp, omega*mp, (1+omega)*mp."""
    g = _ejector_with_outlet_pipe(_solvable_ejector_schema(), _MIXING, 2**0.5 * _A3_PIPE_D)
    net = build_network_from_schema(NetworkGraphSchema(**g))
    seeded = net.elements["ch_o"].initial_guess["ch_o.m_dot"]
    assert seeded == pytest.approx(0.0739, rel=0.15)
    # An A_3 inherited from the outlet pipe no longer skips the P_jct seed.
    g = _ejector_with_outlet_pipe(_solvable_ejector_schema(), None, 2**0.5 * _A3_PIPE_D)
    net = build_network_from_schema(NetworkGraphSchema(**g))
    assert net.elements["ej1"].initial_guess.get("ej1.P_jct", 0.0) > 0.0


@pytest.mark.parametrize("area_factor", [2.0, 5.0])
def test_double_choked_ejector_behind_an_outlet_pipe_converges_to_the_critical_root(area_factor):
    """With A_3 fixed, an outlet pipe wider than the mixing duct must leave the
    critical entrainment ratio unchanged. Before the flow seed the 2x case did
    not converge and the 5x case "converged" to a wrong root (omega 0.29)."""
    g = _ejector_with_outlet_pipe(
        _solvable_ejector_schema(), _MIXING, area_factor**0.5 * _A3_PIPE_D
    )
    df = NetworkRunner.from_dict(g).solve(timeout=120.0).to_dataframe()
    row = df[df["id"] == "ej1"]
    assert bool(row["success"].iloc[0])
    assert float(row["omega"].iloc[0]) == pytest.approx(_critical_omega_bare_outlet(), rel=1e-6)


def test_outlet_port_chamber_takes_a3_when_the_pipe_inherits_it():
    """With both A_3 and the outlet pipe's D left to inherit, the pipe is
    still unsized when the (auto) outlet port chamber resolves. Its area
    search used to hop across the ejector to the primary inlet chamber and
    size the outlet port at 0.05 m^2, running it as a near-plenum."""
    g = _ejector_with_outlet_pipe(_solvable_ejector_schema(), None, None)
    net = build_network_from_schema(NetworkGraphSchema(**g))
    net.resolve_all_topology()
    port = net.nodes["__tee_jct__ej1_outlet"]
    a3 = net.elements["ej1"].mixing_area
    assert port.area == pytest.approx(a3, rel=1e-12)
    assert "ej1" in port._area_source
    # An explicitly sized pipe keeps sizing the port chamber.
    g = _ejector_with_outlet_pipe(_solvable_ejector_schema(), _MIXING, 0.04)
    net = build_network_from_schema(NetworkGraphSchema(**g))
    net.resolve_all_topology()
    assert net.nodes["__tee_jct__ej1_outlet"].area == pytest.approx(3.141592653589793 * 0.02**2)
