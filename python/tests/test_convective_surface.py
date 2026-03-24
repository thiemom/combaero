"""Test ConvectiveSurface dataclass and model subclasses."""

import pytest

import combaero as cb
from combaero.heat_transfer import (
    ConvectiveSurface,
    DimpledModel,
    ImpingementModel,
    PinFinModel,
    RibbedModel,
    SmoothModel,
)


def test_smooth_model_defaults():
    """Verify SmoothModel default values."""
    model = SmoothModel()
    assert model.correlation == "gnielinski"
    assert model.mu_ratio == 1.0
    assert model.roughness == 0.0


def test_ribbed_model_defaults():
    """Verify RibbedModel default values."""
    model = RibbedModel()
    assert model.e_D == 0.0
    assert model.pitch_to_height == 0.0
    assert model.alpha_deg == 90.0


def test_dimpled_model_defaults():
    """Verify DimpledModel default values."""
    model = DimpledModel()
    assert model.d_Dh == 0.0
    assert model.h_d == 0.0
    assert model.S_d == 0.0


def test_pin_fin_model_defaults():
    """Verify PinFinModel default values."""
    model = PinFinModel()
    assert model.pin_diameter == 0.0
    assert model.channel_height == 0.0
    assert model.S_D == 2.0
    assert model.X_D == 2.0
    assert model.N_rows == 1
    assert model.is_staggered is True


def test_impingement_model_defaults():
    """Verify ImpingementModel default values."""
    model = ImpingementModel()
    assert model.d_jet == 0.0
    assert model.z_D == 0.0
    assert model.x_D == 0.0
    assert model.y_D == 0.0
    assert model.A_target == 0.0
    assert model.Cd_jet == 0.8


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
        T_wall=1000.0,
    )

    assert result.h > 0.0
    assert result.T_aw > 700.0  # Should be higher due to recovery
    assert surface.area == 0.1


def test_convective_surface_ribbed_model():
    """Test ConvectiveSurface with RibbedModel."""
    surface = ConvectiveSurface(
        area=0.15,
        model=RibbedModel(e_D=0.05, pitch_to_height=10.0, alpha_deg=60.0),
        heating=True,
    )

    result = surface.htc_and_T(
        T=700.0,
        P=2.5e5,
        X=cb.species.dry_air(),
        velocity=60.0,
        diameter=0.025,
        length=0.6,
        T_wall=1100.0,
    )

    assert result.h > 0.0
    assert result.T_aw > 700.0
    assert surface.area == 0.15


def test_convective_surface_dimpled_model():
    """Test ConvectiveSurface with DimpledModel."""
    surface = ConvectiveSurface(
        area=0.12, model=DimpledModel(d_Dh=0.2, h_d=0.2, S_d=2.0), heating=True
    )

    result = surface.htc_and_T(
        T=650.0,
        P=2e5,
        X=cb.species.dry_air(),
        velocity=45.0,
        diameter=0.022,
        length=0.55,
        T_wall=950.0,
    )

    assert result.h > 0.0
    assert result.T_aw > 650.0
    assert surface.area == 0.12


def test_convective_surface_pin_fin_model():
    """Test ConvectiveSurface with PinFinModel."""
    surface = ConvectiveSurface(
        area=0.08,
        model=PinFinModel(
            pin_diameter=0.003,
            channel_height=0.006,
            S_D=2.5,
            X_D=2.5,
            N_rows=5,
            is_staggered=True,
        ),
        heating=True,
    )

    result = surface.htc_and_T(
        T=600.0,
        P=2e5,
        X=cb.species.dry_air(),
        velocity=40.0,
        diameter=0.003,  # Use pin diameter
        length=0.015,  # 5 rows * 2.5 * 0.003 = 0.0375m streamwise
        T_wall=900.0,
    )

    assert result.h > 0.0
    assert result.T_aw > 600.0
    assert surface.area == 0.08


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
        T_wall=1000.0,
    )

    result_mult = surface_mult.htc_and_T(
        T=700.0,
        P=2e5,
        X=cb.species.dry_air(),
        velocity=50.0,
        diameter=0.02,
        length=0.5,
        T_wall=1000.0,
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
        T_wall=1000.0,
    )

    assert result.h > 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
