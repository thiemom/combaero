"""`MultiPortChamberBase`'s own junction model is gone; its machinery is not.

Deprecated in 0.5.0, removed here. The retirement was decided on measurement:
on the same Bassett separating cells the v1 model scored 0.5260 against
`MultiPortChamberElement`'s 0.0564 and converged on 77 of 105 points against 94. It was
also energetically inconsistent for joining flow and sign-symmetric, so it
admitted mirror roots. Nothing in the package or the GUI instantiated it
(issue #271).

**Only the model went.** The class still owns the topology and port machinery
that `MultiPortChamberElement` and `ConstantKTeeElement` inherit -- port ordering, the
inlet/outlet sign map, area and angle resolution, the connecting-element
wiring checks -- and `diagnostics`, which reports the port map and computes no
physics. `MultiPortChamberElement.diagnostics` calls that one through `super()`.

These tests exist because "remove the model" and "remove the class" are one
careless step apart, and the second would take both shipped junction elements
with it.
"""

from __future__ import annotations

import inspect

import pytest

from combaero.network.components import (
    BorderCarnotLossElement,
    MultiPortChamberBase,
)
from combaero.network.mpce_element import ConstantKTeeElement, MultiPortChamberElement

# ---------------------------------------------------------------------------
# The model is gone
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["residuals", "verify_solution_consistent"])
def test_the_v1_model_methods_are_removed(method):
    assert method not in MultiPortChamberBase.__dict__, (
        f"{method} is back on MultiPortChamberBase; it was removed in 0.6.0"
    )


def test_the_class_is_abstract_and_cannot_be_instantiated():
    """Without a residual it is a base, not an element. Instantiating it used
    to give a working junction, so this is the visible break for anyone who
    had one."""
    with pytest.raises(TypeError, match="abstract"):
        MultiPortChamberBase(
            id="jct",
            inlet_nodes=["a"],
            outlet_nodes=["b", "c"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[0.01] * 3,
        )


def test_the_cpp_kernel_behind_it_is_gone_too():
    """`multi_port_chamber_residuals_and_jacobian` had exactly one caller, the
    method removed above. Leaving it bound would ship a public entry point
    nothing reaches."""
    from combaero import _core, _solver_tools

    assert not hasattr(_core, "multi_port_chamber_residuals_and_jacobian")
    assert not hasattr(_solver_tools, "multi_port_chamber_residuals_and_jacobian")
    assert "multi_port_chamber_residuals_and_jacobian" not in _solver_tools.__all__


# ---------------------------------------------------------------------------
# The machinery is not
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cls", [MultiPortChamberElement, ConstantKTeeElement])
def test_the_shipped_junctions_still_build_and_carry_the_inherited_machinery(cls):
    kwargs = {
        "id": "jct",
        "inlet_nodes": ["a"],
        "outlet_nodes": ["b", "c"],
        "inlet_angles_deg": [0.0],
        "outlet_angles_deg": [0.0, 90.0],
        "port_areas": [0.01] * 3,
        "flow_direction": "branch",
        "strict": False,
    }
    if cls is ConstantKTeeElement:
        kwargs["K_ports"] = {1: 0.35, 2: 1.20}
    element = cls(**kwargs)

    assert isinstance(element, MultiPortChamberBase)
    assert element.N == 3
    assert list(element._port_signs) == [-1.0, 1.0, 1.0]
    assert list(element.port_nodes) == ["a", "b", "c"]


def test_diagnostics_survived_as_machinery():
    """It reports the port map and computes no physics, so it belongs with the
    machinery rather than the model -- and `MultiPortChamberElement` reaches it through
    `super()`, so removing it breaks the element that replaced v1."""
    assert "diagnostics" in MultiPortChamberBase.__dict__
    assert "super()" in inspect.getsource(MultiPortChamberElement.diagnostics)


def test_the_border_carnot_loss_element_is_untouched():
    """It pairs with a junction in the Tier-1 design but is a separate element
    and did not leave with v1. Its own coverage lives in test_area_inference
    and test_network_scenarios."""
    assert BorderCarnotLossElement is not None
    assert "residuals" in BorderCarnotLossElement.__dict__


def test_the_border_carnot_cpp_path_survives():
    """`include/multi_port_chamber.h` was NOT deleted with the model: it also
    holds `border_carnot_L` and `HAGER_FRACTION`, which
    `border_carnot_loss_residual_and_jacobian` is built on. Deleting the
    header wholesale broke the build once while doing this, which is why the
    check is here rather than in a comment.

    `border_carnot_L` itself is C++-internal and not bound, so the bound
    function that depends on it is what gets exercised.
    """
    from combaero import _core, _solver_tools

    assert hasattr(_core, "border_carnot_loss_residual_and_jacobian")
    assert hasattr(_solver_tools, "border_carnot_loss_residual_and_jacobian")
