import numpy as np
import pytest

import combaero as cb
from combaero import HumidAir, dewpoint, humidity_ratio, wet_bulb_temperature


def test_humidair_state_object():
    air = HumidAir()

    T = 300.0
    P = 101325.0
    RH = 0.5

    # Test fluent setter
    assert air.set_TP_RH(T, P, RH) is air

    # Test properties
    assert air.state.T == T
    assert air.state.P == P
    assert air.rh == RH

    # Test read-only getters (with some tolerance due to float conversions)
    assert abs(air.dewpoint - dewpoint(T, P, RH)) < 1e-6
    assert abs(air.wet_bulb - wet_bulb_temperature(T, P, RH)) < 1e-6
    assert abs(air.humidity_ratio - humidity_ratio(T, P, RH)) < 1e-6

    # Note: h_mass is tested fully in C++, we just verify it exists and is a float
    assert isinstance(air.h_mass, float)

    # Test alternate fluent setters
    w = air.humidity_ratio
    air2 = HumidAir()
    assert air2.set_TP_omega(T, P, w) is air2
    assert abs(air2.rh - RH) < 1e-6

    tdp = air.dewpoint
    air3 = HumidAir()
    assert air3.set_TP_dewpoint(T, P, tdp) is air3
    assert abs(air3.rh - RH) < 1e-6

    # Test tuple property getter
    tp_rh = air.TP_RH
    assert tp_rh[0] == T
    assert tp_rh[1] == P
    assert tp_rh[2] == RH

    # Test tuple property setter
    air.TP_RH = (310.0, 100000.0, 0.8)
    assert air.state.T == 310.0
    assert air.state.P == 100000.0
    assert air.rh == 0.8

    # Test TP_dewpoint tuple property setter
    air.TP_dewpoint = (310.0, 100000.0, 300.0)
    assert air.state.T == 310.0
    assert air.state.P == 100000.0
    assert abs(air.dewpoint - 300.0) < 1e-6


def test_humidair_rh_bounds():
    air = HumidAir()
    with pytest.raises(ValueError):
        air.set_TP_RH(300.0, 101325.0, 1.5)  # RH > 1.0

    with pytest.raises(ValueError):
        air.set_TP_RH(300.0, 101325.0, -0.5)  # RH < 0.0


def test_humid_air_rh0_equals_dry_air():
    T, P = 288.15, 101325.0
    humid_rh0 = cb.species.humid_air_mass(T, P, 0.0)
    dry = cb.species.dry_air_mass()
    assert np.allclose(humid_rh0, dry, atol=1e-10), (
        f"humid_air_mass(RH=0) differs from dry_air_mass: max delta = {np.max(np.abs(humid_rh0 - dry))}"
    )


# Reference wet-bulb temperatures from CoolProp HAPropsSI (ASHRAE RP-1485),
# evaluated once offline; the reference is not a project dependency. The
# simplified psychrometric equation used here (constant cp_air, no vapour cp
# term, no enhancement factor) stays within 0.5 K of it on a liquid wet bulb.
@pytest.mark.parametrize(
    ("T", "P", "RH", "Twb_ref"),
    [
        (355.0, 50000.0, 0.0, 290.819),
        (360.0, 50000.0, 0.5, 343.189),
        (370.0, 50000.0, 0.05, 310.174),
        (375.0, 101325.0, 0.0, 304.352),
        (375.0, 101325.0, 0.5, 356.645),
        (380.0, 101325.0, 0.5, 361.021),
        (300.0, 101325.0, 0.5, 292.555),
    ],
)
def test_wet_bulb_converges_when_T_exceeds_boiling_point(
    T: float, P: float, RH: float, Twb_ref: float
) -> None:
    # Above the boiling point at P the saturation humidity ratio is undefined
    # at Twb = T, which is where the old Newton start and clamp put the iterate.
    Twb = wet_bulb_temperature(T, P, RH)
    assert abs(Twb - Twb_ref) < 0.5


def test_dewpoint_below_minus_100C_is_not_pinned_to_the_clamp() -> None:
    # saturation_vapor_pressure accepts T down to 163.15 K; a dewpoint below
    # 173.15 K used to be clamped there and returned as converged.
    T, P, RH = 175.0, 101325.0, 0.1
    Tdp = dewpoint(T, P, RH)
    p_w = RH * cb.saturation_vapor_pressure(T)
    assert cb.saturation_vapor_pressure(Tdp) == pytest.approx(p_w, rel=1e-6)


def test_dewpoint_outside_svp_range_raises() -> None:
    # A dewpoint below the saturation-pressure extrapolation limit cannot be
    # represented; it must refuse, not return the limit.
    with pytest.raises((ValueError, IndexError, RuntimeError)):
        dewpoint(170.0, 101325.0, 1e-4)
