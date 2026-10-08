"""The state relay of a single-source, multi-sink element (branching junctions).

A node's mixed temperature depends on the flows INTO it when it carries heat
(EnergyBoundary): T = T_in + Q/(m cp). For a single-source element the stream
into a sink node is ``flow_at_node(sink)``, so its derivative must be taken
at the sink too. The relay took it at the SOURCE node, so a branching tee's
straight sink saw d/d(m_com) only and lost the -1 on m_branch: that Jacobian
entry was 0 against a finite difference of -1.9.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize._numdiff import approx_derivative

import combaero as cb
from combaero.network import (
    ChannelElement,
    EnergyBoundary,
    FlowNetwork,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    TeeJunctionElement,
)

Y = cb.mole_to_mass(cb.species.dry_air())


def _max_rel_jacobian_error(net: FlowNetwork) -> float:
    s = NetworkSolver(net)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    _, jac = s._residuals_and_jacobian(x)
    fd = approx_derivative(
        lambda v: s._residuals(v), x, method="3-point", abs_step=np.maximum(np.abs(x) * 1e-7, 1e-10)
    )
    mask = np.abs(fd) > 1e-6
    return float((np.abs(jac.toarray() - fd)[mask] / np.abs(fd[mask])).max())


def test_branching_tee_into_a_heated_node() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_common", Pt=2.1e5, Tt=400.0, Y=Y))
    net.add_node(PlenumNode("s"))
    net.add_node(PressureBoundary("pb_out", Pt=2.06e5, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_branch", Pt=2.06e5, Tt=400.0, Y=Y))
    net.nodes["s"].add_energy_boundary(EnergyBoundary("q", Q=5000.0))
    net.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pb_common",
            straight_node="s",
            branch_node="pb_branch",
            theta=math.pi / 2,
            F_C=0.01,
            psi=1.0,
            tee_type="branching",
        )
    )
    net.add_element(ChannelElement("c", "s", "pb_out", length=0.5, diameter=0.1128))
    assert _max_rel_jacobian_error(net) < 1e-3
