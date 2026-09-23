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
        """Standard test geometry: 50mm orifice in 100mm channel."""
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
        beta = geom.d / geom.D
        E = 1.0 / math.sqrt(1.0 - beta**4)
        expected = Cd * E * geom.area * math.sqrt(2 * rho * dP)
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
        # With beta=0.5, E=1.033, so the effective Cd is higher
        # The test geometry has D=0.100, d=0.050, so beta=0.5
        # This means the flow includes the velocity-of-approach factor
        # For a round-trip test, we need to use the same beta in both directions
        assert Cd_calc == pytest.approx(Cd)  # Should return the original Cd


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

        # Small thickness: reattachment benefit
        corr_small = cb.orifice_thickness_correction(0.2, 0.5, Re_d)
        assert corr_small > 1.0

        # Peak around t/d ~ 0.3 (calibrated to Idelchik data)
        corr_peak = cb.orifice_thickness_correction(0.3, 0.5, Re_d)
        assert corr_peak > corr_small

        # Large thickness: friction reduces Cd (non-monotonic behavior)
        corr_long = cb.orifice_thickness_correction(3.0, 0.5, Re_d)
        assert corr_long < corr_peak  # Falls at large t/d due to friction
        assert corr_long < 1.0  # Long-tube behavior: k_t < 1.0


class TestMcGreehanSchotsch1988:
    """McGreehan and Schotsch (1988) composite Cd.

    The C++ suite (tests/test_orifice.cpp) holds the full set of literature
    anchors. These cover the Python surface: the binding, its default, and the
    two behaviours a caller can most easily get wrong.
    """

    def test_bare_jet_plate_agrees_with_florschuetz_default(self):
        # A different paper, rig and decade: Florschuetz, Truman and Metzger
        # (1981) recommend C_D = 0.79 for a plenum-fed jet plate.
        cd = cb.mcgreehan_schotsch_1988_cd(1.0e4, 0.0, 1.0)
        assert cd == pytest.approx(cb.FLORSCHUETZ_1981_DEFAULT_CD, rel=0.01)

    def test_reproduces_the_papers_stated_baseline(self):
        # p.213 states a baseline sharp-edged Cd of 0.60 at Re = 3.2e4.
        assert cb.mcgreehan_schotsch_1988_cd(3.2e4, 0.0, 0.0) == pytest.approx(0.60, abs=5e-4)

    def test_crossflow_defaults_to_zero(self):
        assert cb.mcgreehan_schotsch_1988_cd(1.0e4, 0.0, 1.0) == cb.mcgreehan_schotsch_1988_cd(
            1.0e4, 0.0, 1.0, 0.0
        )

    def test_crossflow_is_not_monotonic(self):
        # Fig. 4 draws a rise to ~0.635 near U1/Vi ~ 0.1 before the decay.
        # Deliberate; see the extraction's item I3.
        base = cb.mcgreehan_schotsch_1988_cd(3.2e4, 0.0, 0.0, 0.0)
        peak = cb.mcgreehan_schotsch_1988_cd(3.2e4, 0.0, 0.0, 0.085)
        far = cb.mcgreehan_schotsch_1988_cd(3.2e4, 0.0, 0.0, 4.0)
        assert peak > base
        assert peak == pytest.approx(0.635, rel=0.02)
        assert far < 0.5 * base

    def test_held_at_the_reynolds_validity_floor(self):
        # Eq. (8) diverges below Re ~ 904, so it is held at its stated floor
        # rather than extrapolated.
        at_floor = cb.mcgreehan_schotsch_1988_cd(1.0e4, 0.0, 1.0)
        for re in (1.0e-3, 1.0, 9.0e2, 5.0e3):
            assert cb.mcgreehan_schotsch_1988_cd(re, 0.0, 1.0) == at_floor
        assert cb.mcgreehan_schotsch_1988_cd(1.0e5, 0.0, 1.0) < at_floor

    def test_stays_inside_florschuetz_measured_band(self):
        # Florschuetz Table 1 measures 0.73-0.85 across his configurations.
        for re in (5.0e3, 1.0e4, 3.0e4, 7.0e4):
            for t_over_d in (1.0, 1.5, 2.0, 3.0):
                cd = cb.mcgreehan_schotsch_1988_cd(re, 0.0, t_over_d)
                assert 0.73 <= cd <= 0.85, f"Re={re}, t/d={t_over_d} -> {cd}"

    def test_a_constant_cd_remains_available_for_the_jet_plate(self):
        # The correlation is opt-in: ImpingementModel.C_D stays a plain float,
        # so a caller can pin a measured or literature value instead.
        from combaero.network.components import ImpingementModel

        assert ImpingementModel().C_D == cb.FLORSCHUETZ_1981_DEFAULT_CD

        measured = ImpingementModel(C_D=0.82)
        assert measured.C_D == 0.82

        computed = ImpingementModel(C_D=cb.mcgreehan_schotsch_1988_cd(1.0e4, 0.0, 2.0))
        assert pytest.approx(0.8211, rel=1e-3) == computed.C_D
