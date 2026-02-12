"""Tests for mass-specific thermodynamic properties.

Verifies Python bindings for cp_mass, cv_mass, h_mass, s_mass.
These functions provide a symmetric interface to molar-basis properties.
Units: All mass-specific properties in [J/kg] or [J/(kg*K)].
"""

import combaero as cb


class TestMassSpecificThermodynamicProperties:
    """Test mass-specific thermodynamic property functions (Phase 2)."""

    def test_cp_mass_air(self):
        """Test mass-specific heat capacity for air."""
        T = 300.0  # K
        X_air = cb.standard_dry_air_composition()

        cp_molar = cb.cp(T, X_air)  # [J/(mol*K)]
        mw = cb.mwmix(X_air)  # [g/mol]
        cp_mass_expected = cp_molar / mw * 1000.0  # [J/(kg*K)]

        cp_mass_actual = cb.cp_mass(T, X_air)

        # Must match within machine precision
        assert abs(cp_mass_actual - cp_mass_expected) < 1e-12

        # Sanity check: air cp_mass ~ 1005 J/(kg*K) at 300K
        assert 1000.0 < cp_mass_actual < 1010.0

    def test_cv_mass_air(self):
        """Test mass-specific heat capacity at constant volume for air."""
        T = 300.0  # K
        X_air = cb.standard_dry_air_composition()

        cv_molar = cb.cv(T, X_air)  # [J/(mol*K)]
        mw = cb.mwmix(X_air)  # [g/mol]
        cv_mass_expected = cv_molar / mw * 1000.0  # [J/(kg*K)]

        cv_mass_actual = cb.cv_mass(T, X_air)

        # Must match within machine precision
        assert abs(cv_mass_actual - cv_mass_expected) < 1e-12

        # Sanity check: air cv_mass ~ 718 J/(kg*K) at 300K
        assert 710.0 < cv_mass_actual < 725.0

    def test_h_mass_air(self):
        """Test mass-specific enthalpy for air."""
        T = 500.0  # K
        X_air = cb.standard_dry_air_composition()

        h_molar = cb.h(T, X_air)  # [J/mol]
        mw = cb.mwmix(X_air)  # [g/mol]
        h_mass_expected = h_molar / mw * 1000.0  # [J/kg]

        h_mass_actual = cb.h_mass(T, X_air)

        # Must match within machine precision
        assert abs(h_mass_actual - h_mass_expected) < 1e-10

    def test_s_mass_air(self):
        """Test mass-specific entropy for air."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()

        s_molar = cb.s(T, X_air, P)  # [J/(mol*K)]
        mw = cb.mwmix(X_air)  # [g/mol]
        s_mass_expected = s_molar / mw * 1000.0  # [J/(kg*K)]

        s_mass_actual = cb.s_mass(T, X_air, P)

        # Must match within machine precision
        assert abs(s_mass_actual - s_mass_expected) < 1e-12

    def test_gamma_consistency(self):
        """Test that gamma = cp/cv is same for molar and mass basis."""
        T = 400.0  # K
        X_air = cb.standard_dry_air_composition()

        # Molar basis
        gamma_molar = cb.cp(T, X_air) / cb.cv(T, X_air)

        # Mass basis
        gamma_mass = cb.cp_mass(T, X_air) / cb.cv_mass(T, X_air)

        # Must be identical (mwmix cancels out)
        assert abs(gamma_molar - gamma_mass) < 1e-14

        # Sanity check: air gamma ~ 1.4
        assert 1.39 < gamma_molar < 1.41

    def test_cp_mass_high_temperature(self):
        """Test cp_mass at high temperature."""
        T = 1500.0  # K
        X_air = cb.standard_dry_air_composition()

        cp_molar = cb.cp(T, X_air)
        mw = cb.mwmix(X_air)
        cp_mass_expected = cp_molar / mw * 1000.0

        cp_mass_actual = cb.cp_mass(T, X_air)

        # Must match within machine precision
        assert abs(cp_mass_actual - cp_mass_expected) / cp_mass_expected < 1e-14

        # At high T, cp increases
        cp_mass_300K = cb.cp_mass(300.0, X_air)
        assert cp_mass_actual > cp_mass_300K

    def test_h_mass_temperature_dependence(self):
        """Test h_mass increases with temperature."""
        X_air = cb.standard_dry_air_composition()

        T1 = 300.0  # K
        T2 = 600.0  # K

        h1 = cb.h_mass(T1, X_air)
        h2 = cb.h_mass(T2, X_air)

        # Enthalpy should increase with temperature
        assert h2 > h1

        # Rough check: delta_h ~ cp_avg * delta_T
        cp_avg = (cb.cp_mass(T1, X_air) + cb.cp_mass(T2, X_air)) / 2.0
        delta_h_expected = cp_avg * (T2 - T1)
        delta_h_actual = h2 - h1

        # Should be close (within 5% due to cp temperature dependence)
        assert abs(delta_h_actual - delta_h_expected) / delta_h_expected < 0.05

    def test_s_mass_pressure_dependence(self):
        """Test s_mass decreases with pressure (at constant T)."""
        T = 300.0  # K
        X_air = cb.standard_dry_air_composition()

        P1 = 101325.0  # 1 atm
        P2 = 202650.0  # 2 atm

        s1 = cb.s_mass(T, X_air, P1)
        s2 = cb.s_mass(T, X_air, P2)

        # Entropy should decrease with increasing pressure
        assert s2 < s1

    def test_all_functions_together(self):
        """Test all 4 mass-specific functions work together."""
        T = 400.0  # K
        P = 200000.0  # Pa
        X_air = cb.standard_dry_air_composition()

        mw = cb.mwmix(X_air)

        # Test all 4 functions
        cp_mass = cb.cp_mass(T, X_air)
        cv_mass = cb.cv_mass(T, X_air)
        h_mass = cb.h_mass(T, X_air)
        s_mass = cb.s_mass(T, X_air, P)

        # Verify against molar basis
        assert abs(cp_mass - cb.cp(T, X_air) / mw * 1000.0) < 1e-12
        assert abs(cv_mass - cb.cv(T, X_air) / mw * 1000.0) < 1e-12
        assert abs(h_mass - cb.h(T, X_air) / mw * 1000.0) < 1e-10
        assert abs(s_mass - cb.s(T, X_air, P) / mw * 1000.0) < 1e-12

        # All should be positive
        assert cp_mass > 0
        assert cv_mass > 0
        assert h_mass > 0  # At 400K, h should be positive
        # s can be negative depending on reference state

    def test_units_are_correct(self):
        """Verify units are in J/kg and J/(kg*K) not J/mol."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()

        cp_mass = cb.cp_mass(T, X_air)
        cv_mass = cb.cv_mass(T, X_air)
        h_mass = cb.h_mass(T, X_air)
        s_mass = cb.s_mass(T, X_air, P)

        # Mass-specific values should be ~29x larger than molar for air
        # (air molar mass ~ 29 g/mol)
        mw_air = cb.mwmix(X_air)
        ratio = mw_air / 1000.0  # Convert g/mol to kg/mol

        cp_molar = cb.cp(T, X_air)
        assert abs(cp_mass * ratio - cp_molar) / cp_molar < 1e-12

        # Typical air values at 300K
        assert 1000 < cp_mass < 1010  # [J/(kg*K)]
        assert 710 < cv_mass < 725  # [J/(kg*K)]
        assert isinstance(h_mass, float)  # [J/kg]
        assert isinstance(s_mass, float)  # [J/(kg*K)]
