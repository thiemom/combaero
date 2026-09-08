"""`ConstantKTeeElement`'s residual and Jacobian, audited before porting it.

#271 covers this element alongside `MPCEv2Element`, and the obvious next step
after #309 was to port it to C++ the same way. Measuring first said otherwise,
and the measurements are pinned here so the decision can be revisited on
evidence rather than memory.

**Its Jacobian is already complete.** `MPCEv2Element`'s was not -- its `K` came
from the Mynard closure and so depended on every port's velocity, while its
hand-derived `dR/dP` column covered only the common port, leaving two columns
absent and wrong by 4.4e-3 and 2.3e-3. Whole-element seeding in C++ supplied
them, which is what made that port worth doing.

Nothing of the sort applies here. `K` is a fixed per-port constant, so `q_dyn`
depends on the common port alone and every other column is genuinely zero. The
entire residual is eight lines of arithmetic with no closure call. A C++ port
would buy consistency of pattern, not correctness and not meaningful speed.

**One thing the audit did turn up**: the `.T` column. It is correct, and no
node type declares a temperature unknown, so the solver drops it. See
`test_the_temperature_column_is_correct_but_currently_unattachable`.

Issue #271.
"""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb
from combaero.network.components import NetworkMixtureState
from combaero.network.mpce_v2_element import ConstantKTeeElement

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_PT_JCT = 2.05e5
_BASE_P = [2.00e5, 1.97e5, 1.95e5]
_BASE_PT = [2.06e5, 2.02e5, 2.00e5]
_BASE_T = [320.0, 319.0, 318.0]
_OUTER = [0.90, 0.55, 0.35]


def _element(flow_direction: str = "branch") -> ConstantKTeeElement:
    if flow_direction == "branch":
        element = ConstantKTeeElement(
            id="jct",
            inlet_nodes=["p0"],
            outlet_nodes=["p1", "p2"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[0.01] * 3,
            flow_direction="branch",
            K_ports={1: 0.35, 2: 1.20},
        )
    else:
        element = ConstantKTeeElement(
            id="jct",
            inlet_nodes=["p0", "p1"],
            outlet_nodes=["p2"],
            inlet_angles_deg=[0.0, 90.0],
            outlet_angles_deg=[0.0],
            port_areas=[0.01] * 3,
            flow_direction="merge",
            K_ports={0: 0.35, 1: 1.20},
        )
    element._port_element_ids = ["e0", "e1", "e2"]
    return element


def _states(P, Pt, T):
    return [
        NetworkMixtureState(P=P[i], Pt=Pt[i], T=T[i], Tt=T[i], m_dot=0.0, Y=list(_Y))
        for i in range(3)
    ]


def _residual(element, P=None, Pt=None, T=None, outer=None) -> np.ndarray:
    P = P or _BASE_P
    Pt = Pt or _BASE_PT
    T = T or _BASE_T
    outer = outer or _OUTER
    mdots = [s * m for s, m in zip(element._port_signs, outer, strict=True)]
    return np.array(element.residuals(_states(P, Pt, T), _PT_JCT, mdots)[0])


def _jacobian(element):
    mdots = [s * m for s, m in zip(element._port_signs, _OUTER, strict=True)]
    return element.residuals(_states(_BASE_P, _BASE_PT, _BASE_T), _PT_JCT, mdots)[1]


def _columns(element):
    return (
        [(f"{element.port_nodes[i]}.P", "P", i, max(_BASE_P[i] * 1e-6, 1.0)) for i in range(3)]
        + [(f"{element.port_nodes[i]}.Pt", "Pt", i, 1.0) for i in range(3)]
        + [(f"{element.port_nodes[i]}.T", "T", i, 1e-4) for i in range(3)]
        + [(f"{element._port_element_ids[i]}.m_dot", "m", i, 1e-7) for i in range(3)]
    )


def _finite_difference(element, kind, index, step) -> np.ndarray:
    def bumped(delta):
        kwargs = {}
        if kind == "P":
            values = list(_BASE_P)
            values[index] += delta
            kwargs["P"] = values
        elif kind == "Pt":
            values = list(_BASE_PT)
            values[index] += delta
            kwargs["Pt"] = values
        elif kind == "T":
            values = list(_BASE_T)
            values[index] += delta
            kwargs["T"] = values
        else:
            values = list(_OUTER)
            values[index] += delta
            kwargs["outer"] = values
        return _residual(element, **kwargs)

    return (bumped(step) - bumped(-step)) / (2.0 * step)


# ---------------------------------------------------------------------------
# The audit result: nothing is missing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("direction", ["branch", "merge"])
def test_every_jacobian_entry_matches_a_finite_difference(direction):
    """The measurement that says this element does not need the #309
    treatment. Every column, both flow directions, against a central
    difference of its own residual."""
    element = _element(direction)
    jac = _jacobian(element)

    for name, kind, index, step in _columns(element):
        fd = _finite_difference(element, kind, index, step)
        for row in range(4):
            reported = jac[row].get(name, 0.0)
            assert abs(reported - fd[row]) <= 1e-5 * max(1.0, abs(fd[row])), (
                f"{direction} row {row} column {name}: "
                f"reported {reported:.6e} against FD {fd[row]:.6e}"
            )


@pytest.mark.parametrize("direction", ["branch", "merge"])
def test_only_the_common_port_enters_the_loss_term(direction):
    """The structural reason the Jacobian is already complete, stated
    directly: K is a fixed constant, so q_dyn depends on the common port
    alone and the other ports' P and mdot columns are genuinely zero -- not
    missing. This is exactly what was NOT true of MPCEv2Element, whose K came
    from the closure and moved with every port's velocity.
    """
    element = _element(direction)
    jac = _jacobian(element)
    common = element._common_port_index()

    for row in range(3):
        for i in range(3):
            if i == common:
                continue
            for suffix in (".P", ".T"):
                name = f"{element.port_nodes[i]}{suffix}"
                assert jac[row].get(name, 0.0) == 0.0, f"row {row} depends on {name}"
            mdot = f"{element._port_element_ids[i]}.m_dot"
            assert jac[row].get(mdot, 0.0) == 0.0, f"row {row} depends on {mdot}"


def test_the_audit_is_not_vacuous():
    """The common port's columns must actually be non-zero, or the test above
    would pass on an element whose Jacobian is empty."""
    element = _element()
    jac = _jacobian(element)
    common = element._common_port_index()

    live = [jac[row].get(f"{element._port_element_ids[common]}.m_dot", 0.0) for row in range(3)]
    assert any(abs(v) > 0.0 for v in live), "no row depends on the common port's flow"


# ---------------------------------------------------------------------------
# The one thing the audit turned up
# ---------------------------------------------------------------------------


def test_the_temperature_column_is_correct_but_currently_unattachable():
    """`dq/dT` is computed and written into the row, and it is right -- the
    test above finite-differences it. But no node class declares a `.T`
    unknown, so `NetworkSolver` never has that name and the entry is dropped
    by the sparse assembly.

    Kept rather than deleted: the derivative is correct, and removing it would
    make the Jacobian silently incomplete the day a node type starts solving
    temperature. Pinned here so the situation is visible instead of being
    rediscovered.
    """
    import inspect

    from combaero.network import components

    declares_temperature = []
    for name in dir(components):
        obj = getattr(components, name)
        if not (inspect.isclass(obj) and hasattr(obj, "unknowns")):
            continue
        try:
            source = inspect.getsource(obj.unknowns)
        except (OSError, TypeError):
            continue
        if '.T"' in source:
            declares_temperature.append(name)

    element = _element()
    jac = _jacobian(element)
    common = element._common_port_index()
    t_name = f"{element.port_nodes[common]}.T"
    assert any(t_name in jac[row] for row in range(3)), "the T column is no longer reported"

    assert not declares_temperature, (
        f"{declares_temperature} now declare a temperature unknown, so the "
        f"ConstantKTee T column is live -- confirm it is being consumed and "
        f"update this test's rationale"
    )


# ---------------------------------------------------------------------------
# The residual itself, which a Jacobian check cannot see
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("direction", ["branch", "merge"])
def test_the_residual_matches_the_documented_closed_form(direction):
    """Comparing an analytic Jacobian against a finite difference of its own
    residual proves the derivative, and NOTHING about the residual: both move
    together when the residual is wrong.

    Found by falsification -- swapping `_common_port_index` for the opposite
    flow direction changed the residual and every test above still passed.

    So the residual is checked against the form the class docstring states,

        Pt_i - Pt_jct + sign * K_i * q_dyn_com

    with `q_dyn_com` recomputed here from the COMMON port identified
    independently, out of the declared inlet/outlet lists rather than out of
    `_common_port_index()`.
    """
    element = _element(direction)
    # Independent: the common leg is the one with a single port on its side.
    common = 0 if direction == "branch" else 2
    assert len(element.inlet_nodes) == (1 if direction == "branch" else 2)
    sign = 1.0 if direction == "branch" else -1.0

    states = _states(_BASE_P, _BASE_PT, _BASE_T)
    rho = float(states[common].density())
    area = float(element.port_areas[common])
    m_common = element._port_signs[common] * _OUTER[common]
    q_dyn = m_common * m_common / (2.0 * rho * area * area)

    expected = []
    for i in range(3):
        K_i = 0.0 if i == common else float(element.K_ports.get(i, 0.0))
        expected.append(_BASE_PT[i] - _PT_JCT + sign * K_i * q_dyn)
    expected.append(sum(element._port_signs[i] * _OUTER[i] for i in range(3)))

    actual = _residual(element)
    for row in range(4):
        assert actual[row] == pytest.approx(expected[row], rel=1e-12, abs=1e-9), (
            f"{direction} row {row}: {actual[row]:.6f} against {expected[row]:.6f}"
        )


@pytest.mark.parametrize(("direction", "expected"), [("branch", 0), ("merge", 2)])
def test_the_common_port_is_the_single_sided_leg(direction, expected):
    """Stated directly, because it is the one input the residual form takes
    that a derivative check cannot reach."""
    assert _element(direction)._common_port_index() == expected
