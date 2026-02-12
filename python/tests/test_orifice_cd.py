"""Tests for orifice Cd correlations."""

import math

import pytest

import combaero as cb


class TestOrificeGeometry:
    """Tests for OrificeGeometry struct."""

    def test_geometry_beta(self):
        geom = cb.OrificeGeometry()
        geom.d = 0.050
        geom.D = 0.100
        assert geom.beta == pytest.approx(0.5)

    def test_geometry_area(self):
        geom = cb.OrificeGeometry()
        geom.d = 0.050
        geom.D = 0.100
        expected = math.pi * 0.050**2 / 4
        assert geom.area == pytest.approx(expected)

    def test_geometry_valid(self):
        geom = cb.OrificeGeometry()
        geom.d = 0.050
        geom.D = 0.100
        assert geom.is_valid()

        # Invalid: d > D
        invalid = cb.OrificeGeometry()
        invalid.d = 0.15
        invalid.D = 0.10
        assert not invalid.is_valid()


class TestCdCorrelations:
    """Tests for Cd correlation functions."""

    @pytest.fixture
    def standard_geom(self):
        """Standard test geometry: 50mm orifice in 100mm pipe."""
        geom = cb.OrificeGeometry()
        geom.d = 0.050
        geom.D = 0.100
        return geom

    @pytest.fixture
    def standard_state(self):
        """Standard flow state."""
        state = cb.OrificeState()
        state.Re_D = 100000.0
        state.dP = 10000.0
        state.rho = 1.2
        state.mu = 1.8e-5
        return state

    def test_Cd_sharp_thin_plate_range(self, standard_geom, standard_state):
        """Cd should be in reasonable range for sharp-edged orifice."""
        Cd = cb.Cd_sharp_thin_plate(standard_geom, standard_state)
        assert 0.58 < Cd < 0.68

    def test_Cd_thick_plate_higher(self, standard_geom, standard_state):
        """Thick plate should have higher Cd than thin plate."""
        standard_geom.t = 0.010  # 10mm thickness
        Cd_thin = cb.Cd_sharp_thin_plate(standard_geom, standard_state)
        Cd_thick = cb.Cd_thick_plate(standard_geom, standard_state)
        assert Cd_thick > Cd_thin
        assert Cd_thick < Cd_thin * 1.3

    def test_Cd_rounded_entry_higher(self, standard_geom, standard_state):
        """Rounded entry should have higher Cd than sharp edge."""
        standard_geom.r = 0.010  # 10mm radius
        Cd_sharp = cb.Cd_sharp_thin_plate(standard_geom, standard_state)
        Cd_round = cb.Cd_rounded_entry(standard_geom, standard_state)
        assert Cd_round > Cd_sharp
        assert Cd_round < 1.0

    def test_Cd_orifice_auto_select(self, standard_geom, standard_state):
        """Auto-select should match thin plate for default geometry."""
        Cd_auto = cb.Cd_orifice(standard_geom, standard_state)
        Cd_thin = cb.Cd_sharp_thin_plate(standard_geom, standard_state)
        assert Cd_auto == pytest.approx(Cd_thin)


class TestOrificeFlowCalculations:
    """Tests for orifice flow calculations with Cd."""

    @pytest.fixture
    def geom(self):
        geom = cb.OrificeGeometry()
        geom.d = 0.050
        geom.D = 0.100
        return geom

    def test_mdot_calculation(self, geom):
        """Test mass flow calculation."""
        Cd = 0.61
        dP = 10000.0
        rho = 1.2
        mdot = cb.orifice_mdot_Cd(geom, Cd, dP, rho)
        expected = Cd * geom.area * math.sqrt(2 * rho * dP)
        assert mdot == pytest.approx(expected)

    def test_dP_round_trip(self, geom):
        """Test dP calculation round-trips with mdot."""
        Cd = 0.61
        dP = 10000.0
        rho = 1.2
        mdot = cb.orifice_mdot_Cd(geom, Cd, dP, rho)
        dP_calc = cb.orifice_dP_Cd(geom, Cd, mdot, rho)
        assert dP_calc == pytest.approx(dP)

    def test_Cd_from_measurement(self, geom):
        """Test Cd extraction from measurement."""
        Cd = 0.61
        dP = 10000.0
        rho = 1.2
        mdot = cb.orifice_mdot_Cd(geom, Cd, dP, rho)
        Cd_calc = cb.orifice_Cd_from_measurement(geom, mdot, dP, rho)
        assert Cd_calc == pytest.approx(Cd)


class TestUtilityFunctions:
    """Tests for utility functions."""

    def test_K_from_Cd_round_trip(self):
        """Test K <-> Cd conversion round-trips."""
        Cd = 0.61
        beta = 0.5
        K = cb.orifice_K_from_Cd(Cd, beta)
        assert K > 0
        Cd_back = cb.orifice_Cd_from_K(K, beta)
        assert Cd_back == pytest.approx(Cd)

    def test_thickness_correction(self):
        """Test thickness correction factor with Idelchik model."""
        Re_d = 1e5  # Typical Reynolds number
        
        # Thin plate: no correction
        assert cb.orifice_thickness_correction(0.01, 0.5, Re_d) == pytest.approx(1.0)

        # Moderate thickness: correction > 1 (reattachment benefit)
        corr1 = cb.orifice_thickness_correction(0.5, 0.5, Re_d)
        assert corr1 > 1.0
        
        # Peak around t/d ~ 1.0
        corr_peak = cb.orifice_thickness_correction(1.0, 0.5, Re_d)
        assert corr_peak > corr1
        
        # Large thickness: friction reduces Cd (non-monotonic behavior)
        corr2 = cb.orifice_thickness_correction(3.0, 0.5, Re_d)
        assert corr2 < corr_peak  # Falls at large t/d due to friction
