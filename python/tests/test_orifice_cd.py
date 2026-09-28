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


class TestIdelchikWallOrifice:
    """Idelchik (1966) diagrams 4-17 and 4-18: a hole in a large wall.

    Ground truth is the source's own tabulated zeta, not the code's output.
    These replace three directional tests over Cd_thick_plate /
    Cd_rounded_entry / Cd_orifice, which asserted only that thick > thin and
    round > sharp and so passed straight over an 11% jump discontinuity at
    Re = 1e5 and a silent fallback to Stolz at r/d = 0.
    """

    def test_sharp_edge_reproduces_diagram_4_17(self) -> None:
        """zeta reaches 2.85 at the table's last knot, Re = 1e6.

        NOT at Re = 1e5: item 1's "Re >= 1e5 -> 2.85" is the coarse statement
        of the same curve the table gives finely, and the table says 2.60 at
        1e5. Treating it as a separate branch put a 4.5% step there.
        """
        hole = cb.DischargeHoleGeometry(d=1e-3)
        sel = cb.DischargeCdCorrelation.Idelchik1966Sharp
        cd = cb.discharge_cd(sel, hole, cb.DischargeHoleState(Re=1e6))
        assert cd == pytest.approx(1.0 / math.sqrt(2.85), rel=1e-9)
        # The table's own value at 1e5, which the old hard switch overrode.
        cd_1e5 = cb.discharge_cd(sel, hole, cb.DischargeHoleState(Re=1e5))
        assert cd_1e5 == pytest.approx(1.0 / math.sqrt(0.04 + 2.56), rel=1e-9)

    def test_rounded_reproduces_diagram_4_18c(self) -> None:
        # Read off the page rendered at 400 dpi, not the OCR layer.
        table = {
            0.01: 2.72,
            0.02: 2.56,
            0.03: 2.40,
            0.04: 2.27,
            0.06: 2.06,
            0.08: 1.88,
            0.12: 1.60,
            0.16: 1.38,
            0.20: 1.37,
        }
        flow = cb.DischargeHoleState(Re=2e5)
        for r_over_d, zeta in table.items():
            hole = cb.DischargeHoleGeometry(d=1e-3, r=r_over_d * 1e-3)
            cd = cb.discharge_cd(cb.DischargeCdCorrelation.Idelchik1966Rounded, hole, flow)
            assert cd == pytest.approx(1.0 / math.sqrt(zeta), rel=1e-9), f"r/d = {r_over_d}"

    def test_thick_hole_reproduces_diagram_4_18a(self) -> None:
        table = {
            0.0: 2.85,
            0.2: 2.72,
            0.4: 2.60,
            0.6: 2.34,
            0.8: 1.95,
            1.0: 1.76,
        }
        # At Re = 1e6 the low-Re form reduces exactly to item 1's
        # zeta' + lam*l/Dh, because k = 1/zeta_sharp (see the header note on
        # the printed 0.342).
        flow = cb.DischargeHoleState(Re=1e6)
        for l_over_d, zeta_prime in table.items():
            hole = cb.DischargeHoleGeometry(d=1e-3, L=l_over_d * 1e-3)
            cd = cb.discharge_cd(cb.DischargeCdCorrelation.Idelchik1966Thick, hole, flow)
            # zeta = zeta' + lam*l/Dh, and lam > 0, so Cd sits just below the
            # friction-free value rather than on it.
            cd_frictionless = 1.0 / math.sqrt(zeta_prime)
            assert cd <= cd_frictionless + 1e-12, f"l/d = {l_over_d}"
            assert cd > cd_frictionless * 0.97, f"l/d = {l_over_d}"

    def test_cd_is_continuous_across_the_1e5_boundary(self) -> None:
        """The removed rounded-entry code jumped 11.1% at Re_D = 1e5.

        0.881391 -> 0.979326 across 9.999e4 -> 1.0e5: a hard C0 break in a
        solver input. Idelchik's own low-Re branch has no such step.
        """
        hole = cb.DischargeHoleGeometry(d=1e-3)
        sel = cb.DischargeCdCorrelation.Idelchik1966Sharp
        below = cb.discharge_cd(hole=hole, correlation=sel, flow=cb.DischargeHoleState(Re=9.999e4))
        above = cb.discharge_cd(hole=hole, correlation=sel, flow=cb.DischargeHoleState(Re=1.0001e5))
        assert abs(above - below) < 1e-4

    def test_low_re_reaches_three_decades_below_mcgreehan(self) -> None:
        """Idelchik's table starts at Re = 25; McGreehan's floor is 1e4."""
        hole = cb.DischargeHoleGeometry(d=1e-3)
        sel = cb.DischargeCdCorrelation.Idelchik1966Sharp
        cd_creep = cb.discharge_cd(sel, hole, cb.DischargeHoleState(Re=25.0))
        cd_min = cb.discharge_cd(sel, hole, cb.DischargeHoleState(Re=400.0))
        cd_turb = cb.discharge_cd(sel, hole, cb.DischargeHoleState(Re=1e6))
        # zeta = 1.94 + 1.00 at Re = 25, 0.54 + 1.37 at 400, 0 + 2.85 at 1e6.
        assert cd_creep == pytest.approx(1.0 / math.sqrt(2.94), rel=1e-9)
        assert cd_min == pytest.approx(1.0 / math.sqrt(1.91), rel=1e-9)
        assert cd_turb == pytest.approx(1.0 / math.sqrt(2.85), rel=1e-9)
        # Cd is NOT monotone in Re: it peaks near Re = 400 where zeta is
        # smallest. Idelchik's own table, not an artifact.
        assert cd_min > cd_creep and cd_min > cd_turb

    def test_crossflow_derivative_is_exactly_zero(self) -> None:
        """Idelchik is plenum-to-plenum: no approach velocity exists.

        An honest absence, not a truncation -- so exactly 0.0, not small.
        """
        hole = cb.DischargeHoleGeometry(d=1e-3, r=1e-4)
        flow = cb.DischargeHoleState(Re=5e4)
        _, _, dcd_du = cb.discharge_cd_and_derivatives(
            cb.DischargeCdCorrelation.Idelchik1966Rounded, hole, flow
        )
        assert dcd_du == 0.0

    def test_agrees_with_mcgreehan_schotsch_cross_source(self) -> None:
        """CROSS-SOURCE accuracy, labelled as such per the validation policy.

        Idelchik (1966) against McGreehan and Schotsch (1988): independent
        sources, 22 years apart. This is an accuracy check, not fidelity --
        it says the two agree, not that either mirrors its own paper.
        """
        Re = 1e5
        sel = cb.DischargeCdCorrelation.Idelchik1966Thick
        for l_over_d in (0.0, 0.4, 1.0, 2.0, 4.0):
            hole = cb.DischargeHoleGeometry(d=1e-3, L=l_over_d * 1e-3)
            cd_i = cb.discharge_cd(sel, hole, cb.DischargeHoleState(Re=Re))
            cd_m = cb.mcgreehan_schotsch_1988_cd(Re, 0.0, l_over_d, 0.0)
            assert cd_m == pytest.approx(cd_i, rel=0.06), f"l/d = {l_over_d}"
        # The sharp-edged anchor is the tight one: 0.5923 vs 0.5926.
        assert cb.mcgreehan_schotsch_1988_cd(Re, 0.0, 0.0, 0.0) == pytest.approx(
            1.0 / math.sqrt(2.85), abs=5e-4
        )

    def test_lichtarowicz_refuses_rather_than_substituting(self) -> None:
        hole = cb.DischargeHoleGeometry(d=1e-3)
        flow = cb.DischargeHoleState(Re=1e5)
        with pytest.raises(ValueError, match="not implemented"):
            cb.discharge_cd(cb.DischargeCdCorrelation.Lichtarowicz1965, hole, flow)


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
        # eps = 0 is the paper exactly: the chain reproduces its stated
        # 0.60 to 5e-4, which is what this test has always pinned.
        base = cb.mcgreehan_schotsch_1988_cd_and_derivatives
        exact = cb.mcgreehan_schotsch_1988_cd(3.2e4, 0.0, 0.0, 0.0, eps=0.0)
        assert exact == pytest.approx(0.60, abs=5e-4)

        # The DEFAULT regularises the crossflow input (rv_smooth_eps), which
        # moves the value at U1/Vi = 0 by ~4e-4 -- 2% of the +/-0.02 scatter
        # the correlation's own data sits in. Stated, not silent.
        assert cb.mcgreehan_schotsch_1988_cd(3.2e4, 0.0, 0.0) == pytest.approx(0.60, abs=1.1e-3)
        assert base(3.2e4, 0.0, 0.0, 0.0)[0] == pytest.approx(0.60, abs=1.1e-3)

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
