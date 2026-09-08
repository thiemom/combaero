"""`MultiPortChamberElement`'s own model is deprecated; its machinery is not.

v0.5.0 is the last release carrying it. The retirement is decided on
measurement, not preference: on the same Bassett separating cells the v1 model
scores 0.5260 against `MPCEv2Element`'s 0.0564 and converges on 77 of 105
points against 94. It is also energetically inconsistent for joining flow and
sign-symmetric, so it admits mirror roots. Nothing in the package or the GUI
instantiates it (issue #271).

The distinction these tests exist to protect is that **only the model goes**.
The class also owns the topology and port machinery that `MPCEv2Element` and
`ConstantKTeeElement` inherit unchanged -- 13 of its 17 public members -- and
that is not deprecated. `MPCEv2Element.diagnostics` in particular delegates
here through `super()`, so a warning guarded on `isinstance` would fire on
every user of the element that REPLACES this one. It is guarded on the exact
type instead.
"""

from __future__ import annotations

import warnings

import pytest

import combaero as cb
from combaero.network.components import (
    BorderCarnotLossElement,
    MultiPortChamberElement,
    NetworkMixtureState,
)
from combaero.network.mpce_v2_element import ConstantKTeeElement, MPCEv2Element

_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _states():
    return [
        NetworkMixtureState(P=1.0e5, Pt=100_300.0, T=300.0, Tt=300.5, m_dot=0.0, Y=list(_Y))
        for _ in range(3)
    ]


def _v1() -> MultiPortChamberElement:
    element = MultiPortChamberElement(
        id="jct",
        inlet_nodes=["p0"],
        outlet_nodes=["p1", "p2"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0],
        port_areas=[0.01] * 3,
    )
    element._port_element_ids = ["e0", "e1", "e2"]
    return element


def _v2(cls=MPCEv2Element):
    kwargs = {
        "id": "jct",
        "inlet_nodes": ["p0"],
        "outlet_nodes": ["p1", "p2"],
        "inlet_angles_deg": [0.0],
        "outlet_angles_deg": [0.0, 90.0],
        "port_areas": [0.01] * 3,
        "flow_direction": "branch",
        "strict": False,
    }
    if cls is ConstantKTeeElement:
        kwargs["K_ports"] = {1: 0.35, 2: 1.20}
    element = cls(**kwargs)
    element._port_element_ids = ["e0", "e1", "e2"]
    return element


_MDOTS = [-0.10, 0.06, 0.04]


# ---------------------------------------------------------------------------
# The model warns
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["residuals", "diagnostics", "verify_solution_consistent"])
def test_each_v1_model_method_warns(method):
    element = _v1()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if method == "residuals":
            element.residuals(_states(), 1.0e5, list(_MDOTS))
        elif method == "diagnostics":
            element.diagnostics(_states(), 1.0e5, list(_MDOTS))
        else:
            element.verify_solution_consistent({})

    deprecations = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert deprecations, f"{method} did not warn"
    text = str(deprecations[0].message)
    assert method in text
    assert "0.6.0" in text, "the warning must say WHEN it goes"
    assert "MPCEv2Element" in text, "the warning must say what to use instead"


def test_it_warns_once_per_element_not_once_per_call():
    """`residuals` is called thousands of times in a Newton loop. A warning
    per call would be unusable, so it is guarded per instance."""
    element = _v1()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for _ in range(20):
            element.residuals(_states(), 1.0e5, list(_MDOTS))

    deprecations = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert len(deprecations) == 1, f"{len(deprecations)} warnings from 20 calls"


def test_a_second_element_warns_again():
    """Per instance, not per process: a second network must still be told."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _v1().residuals(_states(), 1.0e5, list(_MDOTS))
        _v1().residuals(_states(), 1.0e5, list(_MDOTS))

    deprecations = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert len(deprecations) == 2


# ---------------------------------------------------------------------------
# The machinery does not
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cls", [MPCEv2Element, ConstantKTeeElement])
def test_the_replacement_elements_do_not_warn(cls):
    """The load-bearing case. `MPCEv2Element.diagnostics` calls
    `super().diagnostics()`, so a warning guarded on `isinstance` would fire
    on every user of the element that replaces this one."""
    element = _v2(cls)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        element.residuals(_states(), 100_300.0, list(_MDOTS))
        element.diagnostics(_states(), 100_300.0, list(_MDOTS))
        element.verify_solution_consistent({})

    deprecations = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert not deprecations, f"{cls.__name__} warned: {[str(w.message) for w in deprecations]}"


def test_the_subclass_really_does_reach_v1_diagnostics():
    """Otherwise the test above passes for the wrong reason -- it would prove
    nothing if the delegation had quietly gone away."""
    import inspect

    assert "super()" in inspect.getsource(MPCEv2Element.diagnostics)


def test_constructing_v1_does_not_warn():
    """Only USING the model warns. Constructing one is how the machinery gets
    exercised, and a warning at construction would fire for subclasses too."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _v1()
        _v2()

    assert not [w for w in caught if issubclass(w.category, DeprecationWarning)]


def test_the_companion_loss_element_is_not_deprecated():
    """`BorderCarnotLossElement` pairs with v1 in the Tier-1 design but is a
    separate element and is not going anywhere with it."""
    assert BorderCarnotLossElement is not None
    assert not hasattr(BorderCarnotLossElement, "_v1_model_deprecation_warned")
