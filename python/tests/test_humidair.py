import pytest

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
