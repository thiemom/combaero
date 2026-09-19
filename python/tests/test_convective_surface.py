"""Test ConvectiveSurface dataclass and model subclasses."""

import pytest

import combaero as cb
from combaero.heat_transfer import (
    ConvectiveSurface,
    SmoothModel,
)


def test_smooth_model_defaults():
    """Verify SmoothModel default values."""
    model = SmoothModel()
    assert model.correlation == "gnielinski"
    assert model.mu_ratio == 1.0
    assert model.roughness == 0.0


def test_convective_surface_defaults():
    """Verify ConvectiveSurface default values."""
    surface = ConvectiveSurface()
    assert surface.area == 0.0
    assert isinstance(surface.model, SmoothModel)
    assert surface.heating is None
    assert surface.Nu_multiplier == 1.0
    assert surface.f_multiplier == 1.0


def test_convective_surface_disabled():
    """Verify that area=0 disables heat transfer."""
    surface = ConvectiveSurface(area=0.0)
    result = surface.htc_and_T(
        T=700.0,
        P=2e5,
        X=cb.species.dry_air(),
        velocity=50.0,
        diameter=0.02,
        length=0.5,
    )
    assert result is None


def test_convective_surface_smooth_model():
    """Test ConvectiveSurface with SmoothModel."""
    surface = ConvectiveSurface(area=0.1, model=SmoothModel(correlation="gnielinski"), heating=True)

    result = surface.htc_and_T(
        T=700.0,
        P=2e5,
        X=cb.species.dry_air(),
        velocity=50.0,
        diameter=0.02,
        length=0.5,
        T_hot=1000.0,
    )

    assert result.h > 0.0
    assert result.T_aw > 700.0  # Should be higher due to recovery
    assert surface.area == 0.1


def test_convective_surface_multipliers():
    """Test that Nu_multiplier and f_multiplier are applied."""
    surface_base = ConvectiveSurface(area=0.1, model=SmoothModel(), heating=True)

    surface_mult = ConvectiveSurface(area=0.1, model=SmoothModel(), heating=True, Nu_multiplier=1.2)

    result_base = surface_base.htc_and_T(
        T=700.0,
        P=2e5,
        X=cb.species.dry_air(),
        velocity=50.0,
        diameter=0.02,
        length=0.5,
        T_hot=1000.0,
    )

    result_mult = surface_mult.htc_and_T(
        T=700.0,
        P=2e5,
        X=cb.species.dry_air(),
        velocity=50.0,
        diameter=0.02,
        length=0.5,
        T_hot=1000.0,
    )

    # HTC should be scaled by Nu_multiplier
    assert abs(result_mult.h / result_base.h - 1.2) < 1e-6


def test_convective_surface_heating_auto_detect():
    """Test that heating parameter is properly handled."""
    surface = ConvectiveSurface(
        area=0.1, model=SmoothModel(correlation="dittus_boelter"), heating=None
    )

    # When heating=None, it defaults to True in the implementation
    result = surface.htc_and_T(
        T=700.0,
        P=2e5,
        X=cb.species.dry_air(),
        velocity=50.0,
        diameter=0.02,
        length=0.5,
        T_hot=1000.0,
    )

    assert result.h > 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
