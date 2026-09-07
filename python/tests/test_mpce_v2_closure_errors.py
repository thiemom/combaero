"""What happens when the Mynard closure raises inside MPCEv2Element.

Three call sites wrapped ``junction_loss_coefficient`` in ``except Exception``.
The residual site returned a lossless residual with an EMPTY Jacobian for any
exception at all -- which, measured across 2073 scorecard records and the
full suite, never caught a legitimate degenerate state (the guard before the
call handles those) and did catch a plumbing error, disguising it as a
plausible-looking lossless junction (issue #271, step 1.2).

These tests pin the replacement rule: a failure becomes a loud, named one
with the cause chained; anything else propagates untouched.

The rule outlived its original site. Since the whole-element (f, J) moved to
C++ (#271 step 4.4), the residual path no longer calls the Python closure at
all -- it calls ``_core.mpce_v2_residuals_and_jacobian``, which reports a
state it cannot classify as ``valid=False`` rather than by raising. The
residual-path tests below therefore inject at the KERNEL, and the rule they
pin is unchanged: a refusal the guards should have caught first is a named
RuntimeError, never a lossless residual.

The DIAGNOSTICS path still calls the Python closure, so those tests still
inject there.
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero.network import mpce_v2_element as v2
from combaero.network.components import NetworkMixtureState
from combaero.network.mpce_v2_element import MPCEv2Element

_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _element(flow_direction: str = "branch") -> MPCEv2Element:
    element = MPCEv2Element.__new__(MPCEv2Element)
    element.id = "jct"
    element.N = 3
    element.port_nodes = ["p0", "p1", "p2"]
    element.port_areas = [0.01, 0.01, 0.01]
    element.port_angles_deg = [180.0, 0.0, 90.0]
    element.flow_direction = flow_direction
    element._port_signs = [-1.0, 1.0, 1.0] if flow_direction == "branch" else [-1.0, -1.0, 1.0]
    element._port_element_ids = ["e0", "e1", "e2"]
    element.strict = False
    element.joining_etransfer_alpha = 0.2
    element.jacobian_method = "sympy"
    element.penalty_alpha = 0.0
    element.eta_scale = 1.0
    return element


def _states():
    return [
        NetworkMixtureState(P=1.0e5, Pt=100_300.0, T=300.0, Tt=300.5, m_dot=0.0, Y=_Y)
        for _ in range(3)
    ]


_BRANCH_MDOTS = [-0.10, 0.06, 0.04]  # one supplier -> sympy Jacobian path
_MERGE_MDOTS = [-0.06, -0.04, 0.10]  # two suppliers -> FD Jacobian path


def _raising(exc: BaseException, after_calls: int = 0):
    """A closure stand-in that raises ``exc`` on call number ``after_calls``."""
    calls = {"n": 0}
    real = v2.junction_loss_coefficient

    def fake(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] > after_calls:
            raise exc
        return real(*args, **kwargs)

    return fake


# ---------------------------------------------------------------------------
# Residual path
# ---------------------------------------------------------------------------


class _FakeKernelResult:
    """Stand-in for MpceResidualJacobian with a chosen validity."""

    def __init__(self, valid: bool):
        self.valid = valid
        self.residual = [0.0] * 4
        self.jacobian = [[0.0] * 10 for _ in range(4)]
        self.common_port = 0
        self.k_term_sign = 1.0
        self.k_per_port = [0.0] * 3


def test_programming_error_in_the_kernel_propagates(monkeypatch):
    """The 1.2 failure mode, pinned: an AttributeError must surface as itself.

    The old handler turned this into a lossless residual with jac = {} and
    19 tests failed with KeyError instead of the real error.
    """

    def boom(*args, **kwargs):
        raise AttributeError("plumbing")

    monkeypatch.setattr(v2._core, "mpce_v2_residuals_and_jacobian", boom)

    with pytest.raises(AttributeError, match="plumbing"):
        _element().residuals(_states(), 100_300.0, _BRANCH_MDOTS)


def test_a_kernel_refusal_becomes_a_named_failure(monkeypatch):
    """The kernel refuses a state that is not a junction in any regime. The
    guards own every such state, so a refusal here means the two disagree --
    which must be loud, and must name the element and the flows."""
    monkeypatch.setattr(
        v2._core,
        "mpce_v2_residuals_and_jacobian",
        lambda *a, **k: _FakeKernelResult(valid=False),
    )

    with pytest.raises(RuntimeError) as info:
        _element().residuals(_states(), 100_300.0, _BRANCH_MDOTS)

    message = str(info.value)
    assert "MPCEv2Element 'jct'" in message
    assert "port_mdots" in message


def test_residual_never_returns_an_empty_jacobian_on_a_refusal(monkeypatch):
    """The specific defect: a residual with jac = {} looks like a solved row."""
    monkeypatch.setattr(
        v2._core,
        "mpce_v2_residuals_and_jacobian",
        lambda *a, **k: _FakeKernelResult(valid=False),
    )

    with pytest.raises(RuntimeError):
        _element().residuals(_states(), 100_300.0, _BRANCH_MDOTS)
    # Reaching here without a raise means the old fallback is back.


# ---------------------------------------------------------------------------
# What replaced the FD fallback
# ---------------------------------------------------------------------------


def test_no_jacobian_column_is_silently_zero(monkeypatch):
    """The FD fallback's failure mode was a silently zero Jacobian column, and
    two tests used to pin its error handling. The fallback is gone -- the
    kernel differentiates every unknown -- so the intent is pinned directly
    instead: every entry the element reports must match a finite difference of
    its own residual, and no column that should move may be missing.

    Joining flow, which is the case the FD fallback used to serve.
    """
    element = _element("merge")
    states = _states()

    def residual_vector(mdots):
        return list(element.residuals(states, 100_300.0, list(mdots))[0])

    base = residual_vector(_MERGE_MDOTS)
    _, jac = element.residuals(states, 100_300.0, list(_MERGE_MDOTS))

    for j in range(3):
        step = max(abs(_MERGE_MDOTS[j]) * 1e-6, 1e-9)
        hi = list(_MERGE_MDOTS)
        lo = list(_MERGE_MDOTS)
        hi[j] += step
        lo[j] -= step
        fd = [
            (a - b) / (2.0 * step)
            for a, b in zip(residual_vector(hi), residual_vector(lo), strict=True)
        ]
        # port_mdots are junction convention; the reported column is on the
        # OUTER element's variable, hence the sign map.
        name = f"{element._port_element_ids[j]}.m_dot"
        for row in range(4):
            reported = jac[row].get(name, 0.0) * element._port_signs[j]
            assert abs(reported - fd[row]) <= 1e-4 * max(1.0, abs(fd[row])), (
                f"row {row}, column {name}: reported {reported:.6e} against FD {fd[row]:.6e}"
            )
    assert base  # the residual itself is real, not an empty list


# ---------------------------------------------------------------------------
# Diagnostics (post-solve): report, never silently omit
# ---------------------------------------------------------------------------


def test_diagnostics_annotate_a_closure_error_instead_of_dropping_K(monkeypatch):
    monkeypatch.setattr(v2, "junction_loss_coefficient", _raising(ValueError("mask")))

    diag = _element().diagnostics(_states(), 100_300.0, _BRANCH_MDOTS)

    assert diag["closure_error"] == "ValueError"
    assert not any(k.endswith("_K") for k in diag)


def test_diagnostics_let_programming_errors_through(monkeypatch):
    monkeypatch.setattr(v2, "junction_loss_coefficient", _raising(KeyError("plumbing")))

    with pytest.raises(KeyError):
        _element().diagnostics(_states(), 100_300.0, _BRANCH_MDOTS)


# ---------------------------------------------------------------------------
# The legitimate degenerate case still takes the guarded path, untouched
# ---------------------------------------------------------------------------


def test_degenerate_split_is_handled_before_the_closure_is_called(monkeypatch):
    """All ports flowing the same way is caught ahead of the closure call.

    It is a wrong-direction state for at least one port, so it must reach the
    soft barrier -- with a Jacobian (issue #271 step 2) -- and never call the
    closure; the closure stand-in raising here would prove the guard was
    bypassed. It once returned a lossless residual with jac = {}.
    """
    monkeypatch.setattr(
        v2, "junction_loss_coefficient", _raising(IndexError("should not be called"))
    )

    residuals, jac = _element().residuals(_states(), 100_300.0, [-0.1, -0.05, -0.05])

    assert len(residuals) == 4
    assert all(i in jac for i in range(4)), "the fallback must carry its Jacobian"
