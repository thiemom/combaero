"""Global energy and node mass closure of a converged network solve (#481).

The network is the control volume. Enthalpy enters with every stream that
leaves a boundary and leaves with every stream that enters one, at the state
of the node the stream comes FROM -- by the sign of its flow, not the
element's declared direction. Heat enters through the user's
EnergyBoundaries and leaves through walls that put heat on a boundary
(``{node}.Q_wall_out``):

    sum_out m h + sum Q_wall_out - sum_in m h - sum Q_user = 0

conftest.py runs this on every converged ``NetworkSolver.solve`` in the
suite, so every test network is also a conservation test.
"""

from __future__ import annotations

from typing import Any

import numpy as np

import combaero as cb
from combaero.network import MassFlowBoundary, PressureBoundary


def _h(state: Any) -> float:
    return float(cb.h_mass(float(state.Tt), cb.mass_to_mole(list(state.Y))))


def _cpT(state: Any) -> float:
    """Sensible magnitude of a stream's enthalpy, the scale its balance is
    judged on: absolute h is near zero for air at ambient, whatever the flow."""
    return float(cb.cp_mass(float(state.Tt), cb.mass_to_mole(list(state.Y)))) * float(state.Tt)


def closure(solver: Any, result: dict) -> dict[str, float] | None:
    """Imbalances of a converged solve, or None where the check does not apply."""
    net = solver.network
    nodes = net.nodes
    bnd = {n for n, o in nodes.items() if isinstance(o, (PressureBoundary, MassFlowBoundary))}
    for nid, node in nodes.items():
        # A MassFlowBoundary between elements is an injection (#481, A2);
        # 'fraction' energy boundaries are on absolute h (#481, A3).
        if isinstance(node, MassFlowBoundary) and (
            net.get_upstream_elements(nid) and net.get_downstream_elements(nid)
        ):
            return None
        if any(getattr(eb, "fraction", 0.0) for eb in getattr(node, "energy_boundaries", []) or []):
            return None

    x = np.asarray(result["__x_solution__"], dtype=float)
    idx = solver._unknown_indices
    two = solver._two_port_mdot_indices()
    state = {n: solver._get_node_state(o, x) for n, o in nodes.items()}
    h = {n: _h(st) for n, st in state.items()}
    cpT = {n: _cpT(st) for n, st in state.items()}

    e_in = e_out = 0.0
    thermal = 0.0
    imbalance = dict.fromkeys(nodes, 0.0)
    for eid, e in net.elements.items():
        if eid in two:
            m = float(x[two[eid]])
            up, down = (e.from_node, e.to_node) if m >= 0.0 else (e.to_node, e.from_node)
            q = abs(m)
            if up in bnd:
                e_in += q * h[up]
                thermal += q * cpT[up]
            if down in bnd:
                e_out += q * h[up]
                thermal += q * cpT[up]
            imbalance[down] += q
            imbalance[up] -= q
            continue
        ind = idx.get(eid, [])

        def flow(n: str, _e: Any = e, _ind: list[int] = ind) -> float:
            try:
                return float(_e.flow_at_node(n, x, _ind))
            except (IndexError, TypeError):
                return 0.0

        q_src = {n: flow(n) for n in e.all_source_nodes()}
        tot = sum(q_src.values())
        H = sum(q * h[n] for n, q in q_src.items())
        h_e = H / tot if abs(tot) > 1e-300 else 0.0
        for n, q in q_src.items():
            if n in bnd:
                e_in += q * h[n]
                thermal += abs(q) * cpT[n]
            imbalance[n] -= q
        for n in e.all_sink_nodes():
            q = flow(n)
            if n in bnd:
                e_out += q * h_e
            imbalance[n] += q

    q_user = sum(
        eb.Q
        for nid, node in nodes.items()
        for eb in getattr(node, "energy_boundaries", []) or []
        if eb.id != f"_wall_{nid}"
    )
    q_wall_out = sum(
        v for k, v in result.items() if isinstance(k, str) and k.endswith(".Q_wall_out")
    )
    mass = sum(abs(v) for n, v in imbalance.items() if n not in bnd)
    h_span = max((abs(v) for v in h.values()), default=0.0)
    scale = max(abs(e_in), abs(e_out), abs(q_user), thermal, 1.0)
    return {
        "energy": e_out + q_wall_out - e_in - q_user,
        "energy_tol": 1e-7 * scale + mass * (h_span + 1.0),
        "mass": mass,
        "scale": scale,
    }
