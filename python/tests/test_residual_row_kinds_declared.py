"""RESIDUAL-ROW-SCALING: every element declares its rows' units; there is no
default.

A silent "mdot" default misfiled Pa rows 1e5-1e9 too heavily three times
(MPCE impulse rows, the compressible channel #481, four pressure-balance
elements #614). ``declared_row_kinds`` now raises at network validation
unless the units are declared on the class that defines ``residuals`` or a
subclass of it.
"""

from __future__ import annotations

import pytest

import combaero as cb
import combaero.network.ejector_element  # noqa: F401 -- register the subclasses
import combaero.network.ideal_tee_element  # noqa: F401
import combaero.network.mpce_element  # noqa: F401
from combaero.network import FlowNetwork, OrificeElement, PlenumNode, PressureBoundary
from combaero.network import components as comp

Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _subclasses(c: type) -> list[type]:
    out = []
    for s in c.__subclasses__():
        out += [s, *_subclasses(s)]
    return out


SHIPPED = sorted(
    {
        c
        for c in _subclasses(comp.NetworkElement)
        if c.__module__.startswith("combaero.") and c is not comp.MultiPortChamberBase
    },
    key=lambda c: c.__name__,
)


@pytest.mark.parametrize("cls", SHIPPED, ids=lambda c: c.__name__)
def test_every_shipped_element_declares_its_row_units_where_its_residual_is(cls: type) -> None:
    junction = issubclass(cls, comp.MultiPortChamberBase)
    attr = "row_scale_kinds" if junction else "residual_scale_kind"
    res_owner = comp._defining_class(cls, "residuals")
    kind_owner = comp._defining_class(cls, attr)
    assert kind_owner is not comp.NetworkElement, f"{cls.__name__} declares no {attr}"
    assert issubclass(kind_owner, res_owner), (
        f"{cls.__name__}: residuals on {res_owner.__name__}, {attr} on {kind_owner.__name__}"
    )


def _net(element: comp.NetworkElement) -> FlowNetwork:
    g = FlowNetwork()
    for name, Pt in (("a", 1.2e5), ("b", 1.0e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 300.0, Y
        g.add_node(bnd)
    g.add_node(PlenumNode("p"))
    g.add_element(element)
    g.add_element(OrificeElement("o2", "p", "b", Cd=0.8, diameter=0.02, correlation="fixed"))
    return g


def test_a_subclass_that_rewrites_its_residual_must_redeclare() -> None:
    """Inheriting the parent's units unseen is exactly the failure mode."""

    class _NewForm(OrificeElement):
        def residuals(self, state_in, state_out, flows=None):  # noqa: ANN001, ANN201
            return super().residuals(state_in, state_out, flows)

    g = _net(_NewForm("o1", "a", "p", Cd=0.8, diameter=0.02, correlation="fixed"))
    with pytest.raises(TypeError, match="_NewForm.*declare the units"):
        g.validate()

    class _Declared(_NewForm):
        residual_scale_kind = "mdot"

    _net(_Declared("o1", "a", "p", Cd=0.8, diameter=0.02, correlation="fixed")).validate()


def test_a_kind_must_be_p_or_mdot() -> None:
    class _Typo(OrificeElement):
        residual_scale_kind = "pressure"

    g = _net(_Typo("o1", "a", "p", Cd=0.8, diameter=0.02, correlation="fixed"))
    with pytest.raises(TypeError, match="'p' or 'mdot'"):
        g.validate()
