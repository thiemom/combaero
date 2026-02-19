"""Tests for compressible flow functions."""

import numpy as np

import combaero as cb


def test_nozzle_flow_subsonic() -> None:
    """Subsonic nozzle flow (P_back > P_critical)."""
    X = cb.standard_dry_air_composition()
    T0 = 500.0  # K (high enough to stay above 300K after expansion)
    P0 = 200000.0  # Pa
    P_back = 150000.0  # Pa (subsonic)
    A_eff = 0.001  # m^2

    sol = cb.nozzle_flow(T0, P0, P_back, A_eff, X)

    assert sol.M < 1.0  # Subsonic
    assert not sol.choked
    assert sol.mdot > 0
    assert sol.outlet.T < T0  # Temperature drops
    assert sol.outlet.P < P0  # Pressure drops


def test_nozzle_flow_choked() -> None:
    """Choked nozzle flow (P_back < P_critical)."""
    X = cb.standard_dry_air_composition()
    T0 = 600.0  # K (high enough for large expansion)
    P0 = 200000.0  # Pa
    P_back = 50000.0  # Pa (choked)
    A_eff = 0.001  # m^2

    sol = cb.nozzle_flow(T0, P0, P_back, A_eff, X)

    # Flow should be choked (M=1 at throat)
    assert sol.choked
    # Mass flow should be positive
    assert sol.mdot > 0


def test_critical_pressure_ratio() -> None:
    """Critical pressure ratio P*/P0."""
    X = cb.standard_dry_air_composition()
    T0 = 500.0  # K
    P0 = 100000.0  # Pa

    ratio = cb.critical_pressure_ratio(T0, P0, X)

    # For variable-gamma gas, ratio is higher than ideal constant-gamma
    # Actual value ~0.726 for air at 300K
    assert 0.5 < ratio < 0.8


def test_mach_from_pressure_ratio() -> None:
    """Mach number from isentropic pressure ratio."""
    X = cb.standard_dry_air_composition()
    T0 = 500.0  # K
    P0 = 500000.0  # Pa
    P = 80000.0  # Pa

    M = cb.mach_from_pressure_ratio(T0, P0, P, X)

    assert M > 0
    assert M < 1.0  # Subsonic for P/P0 = 0.8


def test_solve_A_eff_from_mdot() -> None:
    """Inverse problem: find area from mass flow."""
    X = cb.standard_dry_air_composition()
    T0 = 500.0
    P0 = 200000.0
    P_back = 150000.0
    A_eff_orig = 0.001

    # Get mass flow for known area
    sol = cb.nozzle_flow(T0, P0, P_back, A_eff_orig, X)
    mdot = sol.mdot

    # Solve inverse
    A_eff_solved = cb.solve_A_eff_from_mdot(T0, P0, P_back, mdot, X)

    assert np.isclose(A_eff_solved, A_eff_orig, rtol=1e-4)


def test_fanno_pipe_basic() -> None:
    """Basic Fanno flow through a pipe."""
    X = cb.standard_dry_air_composition()
    T_in = 400.0  # K
    P_in = 200000.0  # Pa
    u_in = 50.0  # m/s (subsonic)
    L = 1.0  # m
    D = 0.05  # m
    f = 0.02  # Darcy friction factor

    sol = cb.fanno_pipe(T_in, P_in, u_in, L, D, f, X)

    # Pressure should drop due to friction
    assert P_in > sol.outlet.P
    # Velocity should increase (subsonic Fanno)
    assert sol.outlet.T > 0
    assert sol.mdot > 0


def test_fanno_pipe_re_in_populated() -> None:
    """fanno_pipe should populate Re_in and f_avg."""
    X = cb.standard_dry_air_composition()
    T_in = 400.0
    P_in = 200000.0
    u_in = 50.0
    L = 1.0
    D = 0.05
    f = 0.02

    sol = cb.fanno_pipe(T_in, P_in, u_in, L, D, f, X)

    assert sol.Re_in > 0
    # For constant-f overload, f_avg == f
    assert np.isclose(sol.f_avg, f)
    # f field should match input
    assert np.isclose(sol.f, f)


def test_fanno_pipe_profile_has_f_re() -> None:
    """Profile stations should carry f and Re fields."""
    X = cb.standard_dry_air_composition()
    sol = cb.fanno_pipe(400.0, 200000.0, 50.0, 1.0, 0.05, 0.02, X, store_profile=True)

    assert len(sol.profile) > 0
    for st in sol.profile:
        assert st.f >= 0.0
        assert st.Re >= 0.0


def test_fanno_pipe_rough_basic() -> None:
    """fanno_pipe_rough: pressure drops, f_avg and Re_in populated."""
    X = cb.standard_dry_air_composition()
    T_in = 400.0
    P_in = 200000.0
    u_in = 50.0
    L = 1.0
    D = 0.05
    roughness = 0.0  # smooth pipe

    sol = cb.fanno_pipe_rough(T_in, P_in, u_in, L, D, roughness, X)

    assert P_in > sol.outlet.P
    assert sol.mdot > 0
    assert sol.Re_in > 0
    assert sol.f_avg > 0


def test_fanno_pipe_rough_smooth_matches_constant_f() -> None:
    """Smooth pipe (roughness=0): rough variant should be close to constant-f variant
    when f is pre-computed from inlet Re."""
    X = cb.standard_dry_air_composition()
    T_in = 400.0
    P_in = 200000.0
    u_in = 50.0
    L = 1.0
    D = 0.05

    # Compute inlet f from Re using Haaland
    Re_in = cb.reynolds(T_in, P_in, X, u_in, D)
    f_in = cb.friction_haaland(Re_in, 0.0)

    sol_const = cb.fanno_pipe(T_in, P_in, u_in, L, D, f_in, X, n_steps=200)
    sol_rough = cb.fanno_pipe_rough(T_in, P_in, u_in, L, D, 0.0, X, n_steps=200)

    # Outlet pressure should be close (within 1% — Re varies slightly along pipe)
    assert np.isclose(sol_const.outlet.P, sol_rough.outlet.P, rtol=0.01)
    # f_avg should be close to f_in for smooth pipe
    assert np.isclose(sol_rough.f_avg, f_in, rtol=0.05)


def test_fanno_pipe_rough_increases_dP() -> None:
    """Rough pipe should have higher pressure drop than smooth pipe."""
    X = cb.standard_dry_air_composition()
    T_in = 400.0
    P_in = 200000.0
    u_in = 50.0
    L = 1.0
    D = 0.05

    sol_smooth = cb.fanno_pipe_rough(T_in, P_in, u_in, L, D, 0.0, X)
    sol_rough = cb.fanno_pipe_rough(T_in, P_in, u_in, L, D, 1e-4, X)

    dP_smooth = P_in - sol_smooth.outlet.P
    dP_rough = P_in - sol_rough.outlet.P
    assert dP_rough > dP_smooth


def test_fanno_pipe_rough_profile() -> None:
    """Profile stations should carry local f and Re."""
    X = cb.standard_dry_air_composition()
    sol = cb.fanno_pipe_rough(400.0, 200000.0, 50.0, 1.0, 0.05, 1e-4, X, store_profile=True)

    assert len(sol.profile) > 0
    for st in sol.profile:
        assert st.f > 0.0
        assert st.Re > 0.0


def test_fanno_pipe_rough_correlations() -> None:
    """All supported correlations should give similar results for smooth pipe."""
    X = cb.standard_dry_air_composition()
    T_in, P_in, u_in, L, D = 400.0, 200000.0, 50.0, 1.0, 0.05

    sol_h = cb.fanno_pipe_rough(T_in, P_in, u_in, L, D, 0.0, X, correlation="haaland")
    sol_s = cb.fanno_pipe_rough(T_in, P_in, u_in, L, D, 0.0, X, correlation="serghides")
    sol_c = cb.fanno_pipe_rough(T_in, P_in, u_in, L, D, 0.0, X, correlation="colebrook")

    # All three should give very close outlet pressures
    assert np.isclose(sol_h.outlet.P, sol_s.outlet.P, rtol=0.005)
    assert np.isclose(sol_h.outlet.P, sol_c.outlet.P, rtol=0.005)


def test_fanno_pipe_rough_choked() -> None:
    """Long rough pipe should choke: use a very long pipe so choking is certain."""
    X = cb.standard_dry_air_composition()
    T_in = 400.0
    P_in = 200000.0
    u_in = 250.0  # high subsonic velocity — short L* expected
    D = 0.05
    roughness = 1e-4

    # With high velocity and roughness, L* is short; 1000 m is certainly beyond it
    sol = cb.fanno_pipe_rough(T_in, P_in, u_in, 1000.0, D, roughness, X)

    assert sol.choked
    assert sol.L_choke < 1000.0  # choking happened before the end


def test_fanno_pipe_rough_invalid_correlation() -> None:
    """Unknown correlation should raise ValueError."""
    import pytest

    X = cb.standard_dry_air_composition()
    with pytest.raises(ValueError, match="unknown correlation"):
        cb.fanno_pipe_rough(400.0, 200000.0, 50.0, 1.0, 0.05, 0.0, X, correlation="invalid")


def test_fanno_max_length() -> None:
    """Maximum pipe length before choking."""
    X = cb.standard_dry_air_composition()
    T_in = 400.0
    P_in = 200000.0
    u_in = 100.0  # m/s
    D = 0.05
    f = 0.02

    L_star = cb.fanno_max_length(T_in, P_in, u_in, D, f, X)

    assert L_star > 0


def test_friction_haaland() -> None:
    """Haaland friction factor."""
    Re = 50000.0
    e_D = 0.001

    f = cb.friction_haaland(Re, e_D)

    # Should be in reasonable range for turbulent flow
    assert 0.01 < f < 0.1


def test_friction_serghides() -> None:
    """Serghides friction factor."""
    Re = 50000.0
    e_D = 0.001

    f = cb.friction_serghides(Re, e_D)

    assert 0.01 < f < 0.1


def test_friction_colebrook() -> None:
    """Colebrook friction factor."""
    Re = 50000.0
    e_D = 0.001

    f = cb.friction_colebrook(Re, e_D)

    assert 0.01 < f < 0.1


def test_friction_consistency() -> None:
    """All friction correlations should be close."""
    Re = 100000.0
    e_D = 0.0001

    f_haaland = cb.friction_haaland(Re, e_D)
    f_serghides = cb.friction_serghides(Re, e_D)
    f_colebrook = cb.friction_colebrook(Re, e_D)

    # Serghides should be very close to Colebrook (<0.3%)
    assert np.isclose(f_serghides, f_colebrook, rtol=0.003)
    # Haaland within ~3%
    assert np.isclose(f_haaland, f_colebrook, rtol=0.03)


def test_de_laval_nozzle_supersonic() -> None:
    """De Laval nozzle: choked flow with supersonic exit."""
    X = cb.standard_dry_air_composition()

    # Nozzle geometry
    A_inlet = 0.002  # m^2
    A_throat = 0.001  # m^2
    A_exit = 0.0015  # m^2 (area ratio 1.5)
    x_throat = 0.05
    x_exit = 0.10

    T0 = 800.0  # K
    P0 = 500000.0  # Pa
    P_exit = 100000.0  # Pa (low enough for supersonic)

    sol = cb.nozzle_cd(
        T0, P0, P_exit, A_inlet, A_throat, A_exit, x_throat, x_exit, X, n_stations=50
    )

    # Should be choked
    assert sol.choked
    assert sol.mdot > 0

    # Find throat station
    throat_idx = min(range(len(sol.profile)), key=lambda i: abs(sol.profile[i].x - x_throat))
    throat = sol.profile[throat_idx]

    # Throat should be near sonic (M ≈ 1)
    assert 0.9 < throat.M < 1.1

    # Exit should be supersonic (M > 1)
    exit_station = sol.profile[-1]
    assert exit_station.M > 1.0


def test_de_laval_nozzle_subsonic() -> None:
    """De Laval nozzle: subsonic throughout (high back pressure)."""
    X = cb.standard_dry_air_composition()

    A_inlet = 0.002
    A_throat = 0.001
    A_exit = 0.0015
    x_throat = 0.05
    x_exit = 0.10

    T0 = 800.0
    P0 = 500000.0
    P_exit = 400000.0  # High back pressure -> subsonic

    sol = cb.nozzle_cd(
        T0, P0, P_exit, A_inlet, A_throat, A_exit, x_throat, x_exit, X, n_stations=50
    )

    # Should not be choked with high back pressure
    assert not sol.choked

    # Exit should be subsonic
    exit_station = sol.profile[-1]
    assert exit_station.M < 1.0


def test_nozzle_mass_conservation() -> None:
    """Mass flow should be constant along nozzle."""
    X = cb.standard_dry_air_composition()

    A_inlet = 0.002
    A_throat = 0.001
    A_exit = 0.0015

    T0 = 800.0
    P0 = 500000.0
    P_exit = 100000.0

    sol = cb.nozzle_cd(T0, P0, P_exit, A_inlet, A_throat, A_exit, 0.05, 0.10, X, n_stations=50)

    # Check mass flow at several stations: mdot = rho * u * A
    for st in sol.profile[::10]:
        mdot_local = st.rho * st.u * st.A
        assert np.isclose(mdot_local, sol.mdot, rtol=0.01)


def test_nozzle_thrust() -> None:
    """Rocket nozzle thrust calculation."""
    X = cb.standard_dry_air_composition()

    T0 = 2000.0  # K (hot combustion gas)
    P0 = 2000000.0  # Pa (20 bar chamber)
    P_amb = 101325.0  # Pa (sea level)

    A_inlet = 0.002
    A_throat = 0.0005
    A_exit = 0.0012

    sol = cb.nozzle_cd(T0, P0, P_amb, A_inlet, A_throat, A_exit, 0.05, 0.12, X)

    thrust = cb.nozzle_thrust(sol, P_amb)

    # Thrust should be positive
    assert thrust.thrust > 0

    # Isp should be reasonable (100-300 s for hot gas)
    assert 100 < thrust.specific_impulse < 300

    # Thrust coefficient typically 1.0-2.0
    assert 1.0 < thrust.thrust_coefficient < 2.0

    # Verify thrust formula: F = mdot*u_e + (P_e - P_amb)*A_e
    exit_station = sol.profile[-1]
    A_e = exit_station.A
    F_calc = thrust.mdot * thrust.u_exit + (thrust.P_exit - P_amb) * A_e
    assert np.isclose(thrust.thrust, F_calc, rtol=1e-6)


def test_nozzle_thrust_vacuum() -> None:
    """Thrust in vacuum (P_amb = 0)."""
    X = cb.standard_dry_air_composition()

    T0 = 2000.0
    P0 = 2000000.0
    P_amb_sea = 101325.0
    P_amb_vac = 0.0

    sol = cb.nozzle_cd(T0, P0, P_amb_sea, 0.002, 0.0005, 0.0012, 0.05, 0.12, X)

    thrust_sea = cb.nozzle_thrust(sol, P_amb_sea)
    thrust_vac = cb.nozzle_thrust(sol, P_amb_vac)

    # Vacuum thrust should be higher (no back pressure)
    assert thrust_vac.thrust > thrust_sea.thrust
    assert thrust_vac.specific_impulse > thrust_sea.specific_impulse
