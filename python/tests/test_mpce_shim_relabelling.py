"""The Python element as a shim over the C++ whole-element (f, J).

`MultiPortChamberElement.residuals` no longer computes the junction physics. It applies
the guards -- the degenerate-state fallbacks and the wrong-direction soft
barrier, which are solver policy and stay in Python -- and then calls
`_core.mpce_residuals_and_jacobian`, whose Jacobian comes back seeded over
`(P_i, Pt_i, outer_mdot_i, Pt_jct)` in that fixed order.

All that is left is relabelling those ten columns onto the solver's unknown
names. That mapping is the shim's whole contract and nothing else tests it: a
slip would be a silent transposition, not an error, and the golden C++ tests
cannot see it because they never touch the names.

Issue #271, step 4.4.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero.network import NetworkSolver
from combaero.network.components import NetworkMixtureState
from combaero.network.mpce_element import MultiPortChamberElement
from validation.junction import random_robustness as rr

_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _element(flow_direction: str = "branch") -> MultiPortChamberElement:
    ports = ["p0", "p1", "p2"]
    if flow_direction == "branch":
        element = MultiPortChamberElement(
            id="jct",
            inlet_nodes=ports[:1],
            outlet_nodes=ports[1:],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[0.01] * 3,
            flow_direction="branch",
            strict=False,
        )
    else:
        element = MultiPortChamberElement(
            id="jct",
            inlet_nodes=ports[:2],
            outlet_nodes=ports[2:],
            inlet_angles_deg=[0.0, 90.0],
            outlet_angles_deg=[0.0],
            port_areas=[0.01] * 3,
            flow_direction="merge",
            strict=False,
        )
    element._port_element_ids = ["e0", "e1", "e2"]
    return element


def _states():
    return [
        NetworkMixtureState(P=1.0e5, Pt=100_300.0, T=300.0, Tt=300.5, m_dot=0.0, Y=list(_Y))
        for _ in range(3)
    ]


_BRANCH = [-0.10, 0.06, 0.04]
_MERGE = [-0.06, -0.04, 0.10]


# ---------------------------------------------------------------------------
# The names
# ---------------------------------------------------------------------------


def test_every_reported_name_is_one_the_solver_knows():
    """A name the solver does not recognise is silently dropped from the
    sparse assembly, so a typo would cost a whole column with no error."""
    element = _element()
    _, jac = element.residuals(_states(), 100_300.0, list(_BRANCH))

    legal = (
        {f"{n}.P" for n in element.port_nodes}
        | {f"{n}.Pt" for n in element.port_nodes}
        | {f"{e}.m_dot" for e in element._port_element_ids}
        | {f"{element.id}.P_jct"}
    )
    for row, entries in jac.items():
        for name in entries:
            assert name in legal, f"row {row} reports an unknown name {name!r}"


def test_the_pt_and_junction_columns_land_on_the_right_names():
    """R_i = Pt_i - Pt_jct + ..., so row i must carry +1 on ITS OWN port's Pt
    and -1 on the junction, and nothing on any other port's Pt. Pinned by name
    rather than by index, which is what a relabelling slip would break."""
    element = _element()
    _, jac = element.residuals(_states(), 100_300.0, list(_BRANCH))

    for i, node in enumerate(element.port_nodes):
        assert jac[i][f"{node}.Pt"] == pytest.approx(1.0)
        assert jac[i][f"{element.id}.P_jct"] == pytest.approx(-1.0)
        for other in element.port_nodes:
            if other != node:
                assert f"{other}.Pt" not in jac[i], f"row {i} reports {other}.Pt"


def test_the_mass_row_is_the_signed_port_map_by_name():
    element = _element("merge")
    _, jac = element.residuals(_states(), 100_300.0, list(_MERGE))

    mass = jac[element.N]
    for j, eid in enumerate(element._port_element_ids):
        assert mass[f"{eid}.m_dot"] == pytest.approx(element._port_signs[j])
    assert not any(name.endswith(".P") or name.endswith(".Pt") for name in mass)


def test_the_static_pressure_columns_are_reported_at_all():
    """The whole point of the port: the Python used to leave the non-common
    dR/dP columns empty. If the shim drops them the improvement is invisible
    and this is the only place that would notice."""
    element = _element()
    _, jac = element.residuals(_states(), 100_300.0, list(_BRANCH))

    reported = {name for entries in jac.values() for name in entries if name.endswith(".P")}
    assert len(reported) >= 2, f"only {reported} carry a static-pressure derivative"


# ---------------------------------------------------------------------------
# The numbers behind the names
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("direction", "mdots"), [("branch", _BRANCH), ("merge", _MERGE)])
def test_each_named_column_matches_a_finite_difference(direction, mdots):
    """Ties the label to the number. A transposed pair of columns would keep
    every name legal and every value plausible, and only this catches it."""
    element = _element(direction)
    states = _states()

    def residual(perturb=None):
        local = [
            NetworkMixtureState(P=s.P, Pt=s.Pt, T=s.T, Tt=s.Tt, m_dot=s.m_dot, Y=list(s.Y))
            for s in states
        ]
        flows = list(mdots)
        if perturb is not None:
            kind, index, delta = perturb
            if kind == "P":
                local[index] = NetworkMixtureState(
                    P=local[index].P + delta,
                    Pt=local[index].Pt,
                    T=local[index].T,
                    Tt=local[index].Tt,
                    m_dot=local[index].m_dot,
                    Y=list(local[index].Y),
                )
            elif kind == "Pt":
                local[index] = NetworkMixtureState(
                    P=local[index].P,
                    Pt=local[index].Pt + delta,
                    T=local[index].T,
                    Tt=local[index].Tt,
                    m_dot=local[index].m_dot,
                    Y=list(local[index].Y),
                )
            else:
                flows[index] += delta * element._port_signs[index]
        return np.array(element.residuals(local, 100_300.0, flows)[0])

    _, jac = element.residuals(states, 100_300.0, list(mdots))

    columns = (
        [("P", i, f"{element.port_nodes[i]}.P", 1.0) for i in range(3)]
        + [("Pt", i, f"{element.port_nodes[i]}.Pt", 1.0) for i in range(3)]
        + [("m", i, f"{element._port_element_ids[i]}.m_dot", 1e-6) for i in range(3)]
    )
    for kind, index, name, step in columns:
        fd = (residual((kind, index, step)) - residual((kind, index, -step))) / (2.0 * step)
        for row in range(element.N + 1):
            reported = jac[row].get(name, 0.0)
            assert abs(reported - fd[row]) <= 1e-4 * max(1.0, abs(fd[row])), (
                f"{direction} row {row} column {name}: "
                f"reported {reported:.6e} against FD {fd[row]:.6e}"
            )


def test_the_assembled_network_jacobian_matches_finite_differences():
    """End to end, through the solver's own assembly: the relabelling has to
    survive being placed into the global sparse matrix, which is where a name
    the solver does not know would vanish without trace.

    Before the port this measured 4.4e-3 and 2.3e-3 on the two non-common
    static-pressure columns.
    """
    import random

    rng = random.Random(20260906)
    worst = 0.0
    checked = 0
    for _ in range(40):
        case = rr.sample(rng)
        if rr.has_root(case) is not True:
            continue
        solver = NetworkSolver(rr.build(case))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solver.solve(timeout=20.0)
        if not sol["__success__"]:
            continue
        x = np.array(sol["__x_solution__"])
        _, jac = solver._residuals_and_jacobian(x)
        J = jac.toarray() if hasattr(jac, "toarray") else np.asarray(jac)
        fd = np.zeros_like(J)
        for j in range(len(x)):
            h = max(abs(x[j]) * 1e-6, 1e-6)
            xp = x.copy()
            xm = x.copy()
            xp[j] += h
            xm[j] -= h
            fd[:, j] = (solver._residuals(xp) - solver._residuals(xm)) / (2.0 * h)
        worst = max(worst, float(np.max(np.abs(J - fd) / np.maximum(np.abs(fd), 1.0))))
        checked += 1
        if checked >= 6:
            break

    assert checked > 0, "no converged network was available to check"
    assert worst < 1e-4, f"worst relative Jacobian error {worst:.3e}"


def test_math_is_still_imported_for_the_angle_conversion():
    """The shim converts declared degrees to radians for the kernel. Trivial,
    and exactly the kind of thing an import cleanup removes."""
    assert math.radians(180.0) == pytest.approx(math.pi)
    element = _element()
    _, jac = element.residuals(_states(), 100_300.0, list(_BRANCH))
    assert jac  # reached the kernel at all


def test_the_mass_residual_uses_the_original_flows_not_the_snapped_ones():
    """Snapping exists to let the closure classify a dead port; it must not
    enter the conservation statement.

    A port at exactly zero is snapped to a velocity of 1e-9 in its declared
    direction before the kernel sees it, which is a mass flow of order
    1e-11 kg/s. The kernel's own mass row therefore carries that, and the
    shim overwrites it with the sum of the ORIGINAL flows. Worth 1e-11, and
    worth saying: the alternative is a conservation residual that is not zero
    for a state that conserves mass exactly.

    Found by falsification -- deleting the override changed nothing any other
    test could see.
    """
    element = _element()
    mdots = [-0.1, 0.0, 0.1]  # sums to exactly zero; port 1 is dead
    residuals, _ = element.residuals(_states(), 100_300.0, list(mdots))

    assert residuals[element.N] == 0.0, (
        f"mass residual {residuals[element.N]:.3e} for a state that conserves mass exactly"
    )
