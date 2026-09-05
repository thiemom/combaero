"""What happens when the Mynard closure raises inside MPCEv2Element.

Three call sites wrapped ``junction_loss_coefficient`` in ``except Exception``.
The residual site returned a lossless residual with an EMPTY Jacobian for any
exception at all -- which, measured across 2073 scorecard records and the
full suite, never caught a legitimate degenerate state (the guard before the
call handles those) and did catch a plumbing error, disguising it as a
plausible-looking lossless junction (issue #271, step 1.2).

These tests pin the replacement rule: the two exceptions the closure can
raise become a loud, named failure with the cause chained; anything else
propagates untouched. Errors are injected by monkeypatching the closure so
the tests do not depend on finding a real degenerate state.
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


def test_programming_error_in_the_closure_propagates(monkeypatch):
    """The 1.2 failure mode, pinned: an AttributeError must surface as itself.

    The old handler turned this into a lossless residual with jac = {} and
    19 tests failed with KeyError instead of the real error.
    """
    monkeypatch.setattr(v2, "junction_loss_coefficient", _raising(AttributeError("plumbing")))

    with pytest.raises(AttributeError, match="plumbing"):
        _element().residuals(_states(), 100_300.0, _BRANCH_MDOTS)


@pytest.mark.parametrize("exc_type", [IndexError, ValueError])
def test_degenerate_split_error_becomes_a_named_failure(monkeypatch, exc_type):
    """The closure's own two exceptions become a RuntimeError naming the element
    and the state, with the original chained -- not a lossless junction."""
    monkeypatch.setattr(v2, "junction_loss_coefficient", _raising(exc_type("mask")))

    with pytest.raises(RuntimeError) as info:
        _element().residuals(_states(), 100_300.0, _BRANCH_MDOTS)

    message = str(info.value)
    assert "MPCEv2Element 'jct'" in message
    assert "port_mdots" in message
    assert isinstance(info.value.__cause__, exc_type)


def test_residual_never_returns_an_empty_jacobian_on_error(monkeypatch):
    """The specific defect: a residual with jac = {} looks like a solved row."""
    monkeypatch.setattr(v2, "junction_loss_coefficient", _raising(IndexError("mask")))

    with pytest.raises(RuntimeError):
        _element().residuals(_states(), 100_300.0, _BRANCH_MDOTS)
    # If we get here without a raise, the old fallback is back.


# ---------------------------------------------------------------------------
# FD-fallback path (joining flow, two suppliers)
# ---------------------------------------------------------------------------


def test_fd_perturbation_error_becomes_a_named_failure(monkeypatch):
    """First call (the residual's own) succeeds; a perturbed call raises.

    Leaving that Jacobian column zero was a silent wrong derivative.
    """
    monkeypatch.setattr(
        v2, "junction_loss_coefficient", _raising(IndexError("mask"), after_calls=1)
    )

    with pytest.raises(RuntimeError, match="FD perturbation"):
        _element("merge").residuals(_states(), 100_300.0, _MERGE_MDOTS)


def test_fd_perturbation_programming_error_propagates(monkeypatch):
    monkeypatch.setattr(
        v2, "junction_loss_coefficient", _raising(TypeError("plumbing"), after_calls=1)
    )

    with pytest.raises(TypeError, match="plumbing"):
        _element("merge").residuals(_states(), 100_300.0, _MERGE_MDOTS)


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
    """All ports flowing the same way is caught by the guard ahead of the try.

    It must reach the continuity fallback and never call the closure -- the
    closure stand-in raising here would prove the guard was bypassed.
    """
    monkeypatch.setattr(
        v2, "junction_loss_coefficient", _raising(IndexError("should not be called"))
    )

    residuals, jac = _element().residuals(_states(), 100_300.0, [-0.1, -0.05, -0.05])

    assert len(residuals) == 4
    assert jac == {}
