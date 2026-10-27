"""Every name deprecated in 0.7.0 is gone in 0.8.0.

0.7.0 made each deprecated name warn and name its removal release; 0.8.0
removes them. This is the inventory: each removed function, alias, element,
keyword and strategy must fail the way a missing name fails (AttributeError,
ImportError, TypeError, ValueError), not keep working silently, and the
replacements must stay.
"""

from __future__ import annotations

import pytest

import combaero as cb
import combaero.network as net
from combaero import _core, _solver_tools
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
    VortexElement,
    pressure_loss,
)


@pytest.mark.parametrize(
    "name",
    ["standard_dry_air_composition", "smr_wgs_equilibrium", "smr_wgs_equilibrium_adiabatic"],
)
def test_the_removed_functions_are_gone(name: str) -> None:
    assert name not in cb.__all__
    with pytest.raises(AttributeError):
        getattr(cb, name)
    with pytest.raises(AttributeError):
        getattr(_core, name)


def test_the_replacement_functions_remain() -> None:
    assert len(cb.species.dry_air()) == cb.num_species()
    assert callable(cb.reforming_equilibrium)
    assert callable(cb.reforming_equilibrium_adiabatic)


@pytest.mark.parametrize("alias", ["ConstantLoss", "LinearThetaLoss"])
@pytest.mark.parametrize("module", [net, pressure_loss], ids=["network", "pressure_loss"])
def test_the_pressure_loss_aliases_are_gone(module: object, alias: str) -> None:
    with pytest.raises(AttributeError):
        getattr(module, alias)
    assert net.ConstantFractionLoss is pressure_loss.ConstantFractionLoss
    assert net.LinearThetaFractionLoss is pressure_loss.LinearThetaFractionLoss


def test_the_border_carnot_loss_element_is_gone() -> None:
    assert "BorderCarnotLossElement" not in net.__all__
    with pytest.raises(AttributeError):
        net.BorderCarnotLossElement  # noqa: B018
    with pytest.raises(ImportError):
        from combaero.network.components import BorderCarnotLossElement  # noqa: F401
    for module in (_core, _solver_tools):
        assert not hasattr(module, "border_carnot_loss_residual_and_jacobian")
    assert not hasattr(_core, "BorderCarnotLossResult")


@pytest.mark.parametrize(
    "kwargs",
    [{"omega_rpm": 3000.0}, {"r_c": 0.02}, {"n": 2.0}],
    ids=["omega_rpm", "r_c", "n"],
)
def test_the_removed_vortex_keywords_raise(kwargs: dict[str, float]) -> None:
    with pytest.raises(TypeError):
        VortexElement("v", "a", "b", **kwargs)


def test_the_vortex_omega_rpm_attribute_is_gone() -> None:
    el = VortexElement("v", "a", "b", shaft_speed_rpm=3000.0)
    assert el.shaft_speed_rpm == 3000.0
    with pytest.raises(AttributeError):
        el.omega_rpm  # noqa: B018


def test_the_vatistas_functions_remain() -> None:
    """Only the element's r_c/n went; the free functions are unchanged."""
    assert callable(cb.vatistas_v0_bar)


def test_incompressible_warmstart_raises_and_names_the_valid_strategies() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=200000.0, Tt=300.0))
    g.add_node(PressureBoundary("out", Pt=100000.0, Tt=300.0))
    g.add_element(OrificeElement("orf", "in", "out", 0.6, 0.05, correlation="fixed"))
    with pytest.raises(ValueError, match="outletref_warmstart") as info:
        NetworkSolver(g).solve(init_strategy="incompressible_warmstart", timeout=10.0)
    assert "'incompressible_warmstart'" in str(info.value)  # the rejected value
    for valid in ("default", "analytical_pt_prop", "homotopy", "continuation"):
        assert valid in str(info.value)


def test_the_warmstart_maxfev_keyword_went_with_it() -> None:
    with pytest.raises(TypeError):
        NetworkSolver(FlowNetwork()).solve(warmstart_maxfev=150)  # type: ignore[call-arg]


def test_an_unknown_network_attribute_still_raises() -> None:
    with pytest.raises(AttributeError):
        net.NoSuchElement  # noqa: B018


def test_a_gui_file_saved_with_incompressible_warmstart_still_loads() -> None:
    """The library refuses the strategy; a saved GUI file (or a headless
    payload) gets its outlet-referenced successor, as the frontend maps it."""
    from gui.backend.schemas import SolverSettings

    s = SolverSettings(init_strategy="incompressible_warmstart")  # type: ignore[arg-type]
    assert s.init_strategy == "outletref_warmstart"
