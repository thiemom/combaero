"""A discharge-hole Cd that follows its own flow is solved consistently (#582).

Evaluated at the solver's iterate, the row m - Cd(Re(m)) K was non-monotone
in m (IdelchikThick Cd rises 0.37 -> 0.86 between hole Re 25 and 1000), and
|F| had a local minimum at the table's Re floor that trapped Newton on the
GUI-default effusion liner. The element now solves m* = Cd(Re(m*)) K itself;
the row is m - m*(x).
"""

from __future__ import annotations

import numpy as np
import pytest
from test_effusion_wall import _pb, _plate

from combaero.network import FlowNetwork, NetworkSolver


def _single_plate() -> tuple[NetworkSolver, FlowNetwork]:
    g = FlowNetwork()
    g.add_node(_pb("cool", 1.03e5, 600.0))
    g.add_node(_pb("ch", 1.0e5, 1500.0))
    g.add_element(_plate())
    s = NetworkSolver(g)
    s.solve(auto_retry=False)  # sets up the unknowns and the propagation
    return s, g


def test_the_row_is_linear_in_the_element_flow() -> None:
    """Slope exactly 1 and the same m* at every iterate, including the old
    fold (Re 25 to ~1000 per hole) and below the table floor."""
    s, _ = _single_plate()
    rows = []
    for m in (1e-6, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2):
        r, J = s._residuals_and_jacobian(np.array([m]))
        rows.append(m - r[0])
        assert J.toarray()[0, 0] == pytest.approx(1.0, abs=1e-12)
    assert max(rows) - min(rows) <= 1e-12 * max(rows)


def test_the_flow_it_returns_is_self_consistent() -> None:
    """m* = Cd(Re(m*)) K: the C++ flow at Cd evaluated AT m* returns m*."""
    s, g = _single_plate()
    el = g.elements["eff"]
    r, _ = s._residuals_and_jacobian(np.array([1e-3]))
    m_star = 1e-3 - r[0]
    st_in = s._get_node_state(g.nodes["cool"], np.array([m_star]))
    st_out = s._get_node_state(g.nodes["ch"], np.array([m_star]))
    st_in.m_dot = m_star
    cd = el._effective_Cd(st_in, st_out)
    m_calc = el._flow_and_partials(st_in, st_out, cd)[0]
    assert m_calc == pytest.approx(m_star, rel=1e-10)


def test_single_plate_converges_in_one_newton_step() -> None:
    """The exactly linear row reaches |F| ~ 0 at once; MINPACK's progress
    test then misfires at the root, which the solver now accepts."""
    s, _ = _single_plate()
    sol = s.solve(auto_retry=False)
    assert sol["__success__"], sol.get("__message__")
    assert sol["__final_norm__"] < 1e-12
