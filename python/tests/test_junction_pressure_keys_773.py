"""Junction pressures follow one convention (#773): ``P_jct`` static,
``Pt_jct`` total, on every junction model.

Before #773 the Mynard and constant-K tees reported their owned unknown, a
TOTAL pressure, as ``P_jct``, while the ideal momentum tee reported a static
one under the same key (and again as ``p_jn``). The renaming is reporting
only: these tests check that the reported values are the right physical
quantities at a converged solution, on the palette networks the GUI drops.

* Mynard / constant K: ``Pt_jct`` is the owned unknown; the closure has no
  loss on the common port, so the common port's total equals it and the
  static ``P_jct`` is the common port's static pressure.
* Ideal momentum: ``P_jct`` is the owned unknown ``p_jn``; ``Pt_jct`` adds
  the common port's dynamic head.
* Every model: ``Pt_jct - P_jct`` is the common port's dynamic head, which
  at these Mach numbers is 1/2 rho v^2 to within the compressible M^2/4.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "scripts"))

import gui_result_keys as keys  # noqa: E402

from gui.backend.graph_builder import build_network_from_schema  # noqa: E402
from gui.backend.schemas import NetworkGraphSchema  # noqa: E402
from gui.backend.solver_setup import make_solver  # noqa: E402

_TEES = [
    (variant, graph)
    for kind, variant, graph in keys.networks()
    if kind == "mpce_tee" and variant.split(", ")[0] in ("branch", "merge")
]


def _solve(graph: dict):
    schema = NetworkGraphSchema(**graph)
    net = build_network_from_schema(schema)
    solver = make_solver(net, schema)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = solver.solve(timeout=60.0)
    assert sol["__success__"], sol.get("__message__")
    return sol, solver, net


def test_every_junction_model_is_covered() -> None:
    models = {v.split(", ")[1] for v, _ in _TEES}
    assert models == {"mynard", "constant_k", "ideal_momentum"}
    assert {v.split(", ")[0] for v, _ in _TEES} == {"branch", "merge"}


@pytest.mark.parametrize(("variant", "graph"), _TEES, ids=[v for v, _ in _TEES])
def test_P_jct_is_static_and_Pt_jct_total(variant: str, graph: dict) -> None:
    direction, model = variant.split(", ")
    sol, solver, net = _solve(graph)
    jct = net.elements["x"]
    diag = sol["__element_diag__"]["x"]
    assert "p_jn" not in diag
    x = float(solver._last_x[solver._unknown_indices["x"][0]])
    common = 0 if direction == "branch" else jct.N - 1
    states = [solver._get_node_state(net.nodes[n], solver._last_x) for n in jct.port_nodes]
    st = states[common]

    if model == "ideal_momentum":
        assert diag["P_jct"] == x  # the owned unknown is the static p_jn
    else:
        assert diag["Pt_jct"] == x  # the owned unknown is the total
        assert st.Pt == pytest.approx(x, abs=1e-6 * x)  # K = 0 on the common port
        assert diag["P_jct"] == pytest.approx(st.P, abs=1e-6 * x)
    if model == "ideal_momentum" and direction == "branch":
        # The lone supplier sits at p_jn, so its total is the junction's.
        assert diag["P_jct"] == pytest.approx(st.P, abs=1e-6 * x)
        assert diag["Pt_jct"] == pytest.approx(st.Pt, abs=1e-6 * x)

    # Pt - P = the common port's dynamic head, 1/2 rho v^2 at low Mach.
    rho = jct.port_densities(states)[common][0]
    m_c = abs(float(sol[f"{jct._port_element_ids[common]}.m_dot"]))
    v = m_c / (rho * float(jct.port_areas[common]))
    q_low_mach = 0.5 * rho * v * v
    mach = float(sol[f"{jct.port_nodes[common]}.mach"])
    assert 0.0 < mach < 0.4
    assert diag["Pt_jct"] - diag["P_jct"] > 0.0
    assert diag["Pt_jct"] - diag["P_jct"] == pytest.approx(q_low_mach, rel=0.3 * mach * mach + 1e-9)
    assert diag["Pt_jct"] - diag["P_jct"] == pytest.approx(st.Pt - st.P, abs=1e-6 * x)


def test_the_ejector_reports_its_mixing_plane_static_once() -> None:
    """P_jct is the mixing-plane static P_py; the duplicate p_py_pa is gone,
    and there is no Pt_jct (the plane has two streams at one static pressure,
    hence no single total)."""
    graph = next(g for k, v, g in keys.networks() if k == "ejector" and v == "palette")
    sol, solver, _ = _solve(graph)
    diag = sol["__element_diag__"]["x"]
    x = float(solver._last_x[solver._unknown_indices["x"][0]])
    assert diag["P_jct"] == x
    for old in ("p_py_pa", "outlet_pt_pa", "p_c_star_pa", "p_b0_pa", "omega", "Pt_jct"):
        assert old not in diag, old
    assert diag["Pt_out"] > 0.0
    assert diag["entrainment_ratio"] == pytest.approx(
        diag["m_dot_secondary"] / diag["m_dot_primary"], rel=1e-12
    )


@pytest.mark.parametrize(("variant", "graph"), _TEES, ids=[v for v, _ in _TEES])
def test_the_unknown_is_named_after_what_it_holds(variant: str, graph: dict) -> None:
    """The owned unknown is Pt_jct where it is a total (Mynard, constant K)
    and P_jct where it is the static (ideal momentum), so the solver variable
    and the result key of one name are one quantity."""
    model = variant.split(", ")[1]
    sol, solver, net = _solve(graph)
    owned = [n for n in solver.unknown_names if n.startswith("x.")]
    expected = "x.P_jct" if model == "ideal_momentum" else "x.Pt_jct"
    assert owned == [expected]
    x = float(solver._last_x[solver._unknown_indices["x"][0]])
    assert sol[expected] == x
    assert sol["x.Pt_jct"] >= sol["x.P_jct"]


def _mynard_solver(guess: dict):
    graph = next(g for v, g in _TEES if v == "branch, mynard")
    graph = {**graph, "nodes": [dict(n) for n in graph["nodes"]]}
    schema = NetworkGraphSchema(**graph)
    net = build_network_from_schema(schema)
    net.elements["x"].initial_guess = guess
    solver = make_solver(net, schema)
    net.resolve_all_topology()
    return solver


def test_an_old_P_jct_guess_no_longer_seeds_the_mynard_tee() -> None:
    """Before #773 the Mynard tee's total was the unknown 'x.P_jct'. The old
    key seeded 'x.Pt_jct', with a DeprecationWarning, until 0.8.0; now only
    the current key does."""
    new = _mynard_solver({"x.Pt_jct": 1.234e5})
    assert new._build_x0()[new.unknown_names.index("x.Pt_jct")] == 1.234e5
    old = _mynard_solver({"x.P_jct": 1.234e5})
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        x0 = old._build_x0()
    assert x0[old.unknown_names.index("x.Pt_jct")] != 1.234e5
    assert not hasattr(old.network.elements["x"], "legacy_unknown_names")


def test_a_gui_file_with_the_old_junction_guess_loads_silently() -> None:
    """A saved GUI tee carries the short key; the backend maps P_jct to
    Pt_jct for the Mynard and constant-K models without a warning, and leaves
    the ideal tee's static P_jct alone."""
    for model, key in (
        ("mynard", "Pt_jct"),
        ("constant_k", "Pt_jct"),
        ("ideal_momentum", "P_jct"),
    ):
        graph = next(g for v, g in _TEES if v == f"branch, {model}")
        nodes = []
        for n in graph["nodes"]:
            n = dict(n)
            if n["id"] == "x":
                n["data"] = {**n["data"], "initial_guess": {"P_jct": 1.1e5}}
            nodes.append(n)
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            net = build_network_from_schema(NetworkGraphSchema(**{**graph, "nodes": nodes}))
        assert net.elements["x"].initial_guess == {f"x.{key}": 1.1e5}, model
