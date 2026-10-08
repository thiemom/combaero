"""Global energy and node mass closure of a converged network solve (#481).

The network is the control volume. Enthalpy enters with every stream that
leaves a boundary and leaves with every stream that enters one, at the state
of the node the stream comes FROM -- by the sign of its flow, not the
element's declared direction. Heat enters through the user's
EnergyBoundaries (Q, and ``{node}.Q_fraction`` for 'fraction'), leaves
through walls that put heat on a boundary (``{node}.Q_wall_out``), and is
not taken up where a node cannot absorb it (``{node}.Q_withheld``):

    sum_out m h + sum Q_wall_out + sum Q_withheld - sum_in m h - sum Q_user = 0

conftest.py runs this on every converged ``NetworkSolver.solve`` in the
suite, so every test network is also a conservation test.
"""

from __future__ import annotations

from typing import Any

import numpy as np

import combaero as cb
from combaero.network import CombustorNode, MassFlowBoundary, PressureBoundary


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
    # A MassFlowBoundary between elements is an injection: an interior node
    # with a source stream of its own (m_dot at its Tt and Y).
    injections = {
        n: o
        for n, o in nodes.items()
        if isinstance(o, MassFlowBoundary)
        and net.get_upstream_elements(n)
        and net.get_downstream_elements(n)
    }
    bnd = {
        n
        for n, o in nodes.items()
        if isinstance(o, (PressureBoundary, MassFlowBoundary)) and n not in injections
    }
    for node in nodes.values():
        # A combustor's 'fraction' acts on its products' sensible enthalpy,
        # which the result does not report.
        if isinstance(node, CombustorNode) and any(
            getattr(eb, "fraction", 0.0) for eb in getattr(node, "energy_boundaries", []) or []
        ):
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

    for n, o in injections.items():
        Y = o.Y if o.Y is not None else list(cb.mole_to_mass(cb.species.dry_air()))
        X = cb.mass_to_mole(list(Y))
        e_in += float(o.m_dot) * float(cb.h_mass(float(o.Tt), X))
        thermal += abs(float(o.m_dot)) * float(cb.cp_mass(float(o.Tt), X)) * float(o.Tt)
        imbalance[n] += float(o.m_dot)

    q_user = sum(
        eb.Q
        for nid, node in nodes.items()
        for eb in getattr(node, "energy_boundaries", []) or []
        if eb.id != f"_wall_{nid}"
    )
    q_user += sum(v for k, v in result.items() if isinstance(k, str) and k.endswith(".Q_fraction"))
    q_wall_out = sum(
        v for k, v in result.items() if isinstance(k, str) and k.endswith(".Q_wall_out")
    )
    q_withheld = sum(
        v for k, v in result.items() if isinstance(k, str) and k.endswith(".Q_withheld")
    )
    mass = sum(abs(v) for n, v in imbalance.items() if n not in bnd)
    h_span = max((abs(v) for v in h.values()), default=0.0)
    scale = max(abs(e_in), abs(e_out), abs(q_user), thermal, 1.0)
    return {
        "energy": e_out + q_wall_out + q_withheld - e_in - q_user,
        "energy_tol": 1e-7 * scale + mass * (h_span + 1.0),
        "mass": mass,
        "scale": scale,
    }
