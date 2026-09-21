"""The GUI must reject removed surface types, not silently substitute smooth.

Enhanced-surface correlations were removed in 0.7.0 (issue #339) because they
could not be traced to their cited sources. A saved network asking for a ribbed
channel and getting an unribbed one back would be a wrong answer presented as a
right one, so `map_surface_model` raises instead of falling through.

Ribbed returned in 0.8.0 on a provenanced correlation set; impingement (both
the jet array and the single-jet form) returned in 0.9.0 on Florschuetz,
Truman and Metzger (1981) with a real crossflow term (#337). Dimpled and
pin-fin remain deferred indefinitely (#339). This file pins three things:
that ribbed and both impingement forms map again, that the still-deferred two
refuse, and that every field survives the JSON round trip.
"""

from __future__ import annotations

import pytest

from gui.backend.graph_builder import map_surface_model


class _Surface:
    """Minimal stand-in: the schema types for these were removed with them."""

    def __init__(self, type_: str) -> None:
        self.type = type_


def test_smooth_still_maps() -> None:
    from combaero.network.components import SmoothModel

    assert isinstance(map_surface_model(_Surface("smooth")), SmoothModel)


def test_ribbed_maps_again_and_tolerates_an_older_saved_network() -> None:
    """Ribbed returned in 0.8.0 on a provenanced correlation set.

    A network saved before the newer fields existed carries only `type`, so the
    mapping falls back to the schema defaults rather than raising an
    AttributeError the user cannot act on.
    """
    from combaero.network.components import RibbedModel

    m = map_surface_model(_Surface("ribbed"))
    assert isinstance(m, RibbedModel)
    assert m.n_ribbed_walls == 2
    assert m.smooth_wall_Nu_multiplier == 1.0


def test_impingement_maps_again_and_tolerates_an_older_saved_network() -> None:
    """Impingement (jet array) returned in 0.9.0 on Florschuetz, Truman and
    Metzger (1981) with a real crossflow term -- see issue #337."""
    from combaero.network.components import ImpingementModel

    m = map_surface_model(_Surface("impingement"))
    assert isinstance(m, ImpingementModel)
    assert m.row == 1
    assert pytest.approx(0.79) == m.C_D


def test_single_jet_impingement_maps_again() -> None:
    import combaero as cb
    from combaero.network.components import SingleJetImpingementModel

    m = map_surface_model(_Surface("single_jet_impingement"))
    assert isinstance(m, SingleJetImpingementModel)
    assert m.bc == cb.ImpingementThermalBC.ConstantHeatFlux
    assert pytest.approx(7.75) == m.L_D


@pytest.mark.parametrize("surface", ["dimpled", "pin_fin"])
def test_deferred_surfaces_report_removal(surface: str) -> None:
    with pytest.raises(ValueError, match="removed in 0.7.0") as exc:
        map_surface_model(_Surface(surface))
    assert "#339" in str(exc.value)
    assert "temporarily" not in str(exc.value)


def test_unknown_surface_is_rejected_too() -> None:
    """A type nobody recognises must not quietly become smooth either."""
    with pytest.raises(ValueError):
        map_surface_model(_Surface("not_a_real_surface"))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


# ---------------------------------------------------------------------------
# The ribbed round trip, restored in 0.8.0
# ---------------------------------------------------------------------------


def test_ribbed_survives_the_json_to_element_round_trip() -> None:
    """A saved network's JSON must reach the element with every field intact.

    The schema discriminates the surface union, so `map_surface_model` receives
    a typed `RibbedModelData` rather than a dict. Pinned because the fields are
    easy to lose one at a time -- `p_e` was called `pitch_to_height` before the
    rebuild, and a silent default would look like a working network giving
    quietly wrong answers.
    """
    from combaero.network.components import RibbedModel
    from gui.backend.schemas import ChannelData

    data = ChannelData(
        length=0.6,
        diameter=0.025,
        surface={
            "type": "ribbed",
            "e_D": 0.06,
            "p_e": 12.0,
            "alpha_deg": 60.0,
            "W_H": 2.0,
            "n_ribbed_walls": 4,
            "smooth_wall_Nu_multiplier": 1.67,
        },
    )
    model = map_surface_model(data.surface)
    assert isinstance(model, RibbedModel)
    assert model.e_D == 0.06
    assert model.p_e == 12.0
    assert model.alpha_deg == 60.0
    assert model.W_H == 2.0
    assert model.n_ribbed_walls == 4
    assert model.smooth_wall_Nu_multiplier == 1.67


def test_a_channel_with_no_surface_still_defaults_to_smooth() -> None:
    from combaero.network.components import SmoothModel
    from gui.backend.schemas import ChannelData

    data = ChannelData(length=0.6, diameter=0.025)
    assert isinstance(map_surface_model(data.surface), SmoothModel)


# ---------------------------------------------------------------------------
# Impingement, restored in 0.9.0
# ---------------------------------------------------------------------------


def test_impingement_survives_the_json_to_element_round_trip() -> None:
    """Same rationale as the ribbed round trip above: every field must reach
    the element intact, not silently default -- xn_d/yn_d/z_d/row/C_D are all
    new, easy to typo or drop one at a time."""
    from combaero.network.components import ImpingementModel
    from gui.backend.schemas import ChannelData

    data = ChannelData(
        length=0.6,
        diameter=0.025,
        surface={
            "type": "impingement",
            "d_jet": 0.0025,
            "xn_d": 10.0,
            "yn_d": 5.0,
            "z_d": 1.5,
            "row": 4,
            "C_D": 0.82,
        },
    )
    model = map_surface_model(data.surface)
    assert isinstance(model, ImpingementModel)
    assert model.d_jet == 0.0025
    assert model.xn_d == 10.0
    assert model.yn_d == 5.0
    assert model.z_d == 1.5
    assert model.row == 4
    assert pytest.approx(0.82) == model.C_D


def test_single_jet_impingement_survives_the_json_to_element_round_trip() -> None:
    import combaero as cb
    from combaero.network.components import SingleJetImpingementModel
    from gui.backend.schemas import ChannelData

    data = ChannelData(
        length=0.6,
        diameter=0.025,
        surface={
            "type": "single_jet_impingement",
            "bc": "constant_wall_temperature",
            "d_jet": 0.004,
            "L_D": 6.0,
            "R_D": 3.0,
        },
    )
    model = map_surface_model(data.surface)
    assert isinstance(model, SingleJetImpingementModel)
    assert model.bc == cb.ImpingementThermalBC.ConstantWallTemperature
    assert model.d_jet == 0.004
    assert model.L_D == 6.0
    assert model.R_D == 3.0
