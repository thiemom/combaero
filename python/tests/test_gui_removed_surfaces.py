"""The GUI must reject removed surface types, not silently substitute smooth.

Enhanced-surface correlations were removed in 0.7.0 (issue #339) because they
could not be traced to their cited sources. A saved network asking for a ribbed
channel and getting an unribbed one back would be a wrong answer presented as a
right one, so `map_surface_model` raises instead of falling through.

The two messages differ on purpose. Ribbed returns once a provenanced
correlation lands (#334); the other three are deferred indefinitely. Telling a
user with a ribbed network that their surface is gone for good would be false.
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


def test_ribbed_says_it_is_coming_back() -> None:
    with pytest.raises(ValueError, match="temporarily unavailable") as exc:
        map_surface_model(_Surface("ribbed"))
    assert "#334" in str(exc.value)


@pytest.mark.parametrize("surface", ["dimpled", "pin_fin", "impingement"])
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
