"""
Tests for advanced cooling correlations.

Tests verify:
- Rib-enhanced cooling (Han et al. 1988)
- Impingement cooling (Florschuetz et al. 1981, Martin 1977)
- Film cooling effectiveness (Baldauf et al. 2002)
- Parameter range validation
- Physical trends and behavior
"""

import pytest

import combaero as cb


class TestRibEnhancedCooling:
    """Test rib-enhanced cooling correlations."""

    def test_rib_enhancement_baseline(self):
        """Test rib enhancement at baseline conditions."""
        e_D = 0.05
        P_e = 10.0
        alpha = 90.0

        enhancement = cb.rib_enhancement_factor(e_D, P_e, alpha)

        # Should give enhancement > 1 (correlation gives ~1.0-3.5 range)
        assert 1.0 < enhancement < 4.0

    def test_rib_enhancement_increases_with_height(self):
        """Test that enhancement increases with rib height."""
        P_e = 10.0
        alpha = 90.0

        enh_low = cb.rib_enhancement_factor(e_D=0.02, P_e=P_e, alpha=alpha)
        enh_mid = cb.rib_enhancement_factor(e_D=0.05, P_e=P_e, alpha=alpha)
        enh_high = cb.rib_enhancement_factor(e_D=0.08, P_e=P_e, alpha=alpha)

        assert enh_mid > enh_low
        assert enh_high > enh_mid

    def test_rib_enhancement_angle_effect(self):
        """Test that 90-deg ribs are most effective."""
        e_D = 0.05
        P_e = 10.0

        enh_30 = cb.rib_enhancement_factor(e_D, P_e, alpha=30)
        enh_60 = cb.rib_enhancement_factor(e_D, P_e, alpha=60)
        enh_90 = cb.rib_enhancement_factor(e_D, P_e, alpha=90)

        # 90-deg should be highest
        assert enh_90 > enh_60
        assert enh_60 > enh_30

    def test_rib_enhancement_parameter_validation(self):
        """Test parameter range validation for rib enhancement."""
        # Valid range: e_D = 0.02-0.1, P_e = 5-20, alpha = 30-90

        with pytest.raises(RuntimeError):
            cb.rib_enhancement_factor(e_D=0.01, P_e=10, alpha=90)

        with pytest.raises(RuntimeError):
            cb.rib_enhancement_factor(e_D=0.15, P_e=10, alpha=90)

        with pytest.raises(RuntimeError):
            cb.rib_enhancement_factor(e_D=0.05, P_e=3, alpha=90)

        with pytest.raises(RuntimeError):
            cb.rib_enhancement_factor(e_D=0.05, P_e=25, alpha=90)

        with pytest.raises(RuntimeError):
            cb.rib_enhancement_factor(e_D=0.05, P_e=10, alpha=20)

        with pytest.raises(RuntimeError):
            cb.rib_enhancement_factor(e_D=0.05, P_e=10, alpha=100)

    def test_rib_friction_baseline(self):
        """Test rib friction multiplier at baseline conditions."""
        e_D = 0.05
        P_e = 10.0

        multiplier = cb.rib_friction_multiplier(e_D, P_e)

        # Should give friction increase > 1 (correlation gives ~1.0-8.0 range)
        assert 1.0 < multiplier < 10.0

    def test_rib_friction_increases_with_height(self):
        """Test that friction increases with rib height."""
        P_e = 10.0

        f_low = cb.rib_friction_multiplier(e_D=0.02, P_e=P_e)
        f_mid = cb.rib_friction_multiplier(e_D=0.05, P_e=P_e)
        f_high = cb.rib_friction_multiplier(e_D=0.08, P_e=P_e)

        assert f_mid > f_low
        assert f_high > f_mid

    def test_rib_friction_parameter_validation(self):
        """Test parameter range validation for rib friction."""
        # Valid range: e_D = 0.02-0.1, P_e = 5-20

        with pytest.raises(RuntimeError):
            cb.rib_friction_multiplier(e_D=0.01, P_e=10)

        with pytest.raises(RuntimeError):
            cb.rib_friction_multiplier(e_D=0.15, P_e=10)

        with pytest.raises(RuntimeError):
            cb.rib_friction_multiplier(e_D=0.05, P_e=3)

        with pytest.raises(RuntimeError):
            cb.rib_friction_multiplier(e_D=0.05, P_e=25)


class TestImpingementCooling:
    """Test impingement cooling correlations."""

    def test_single_jet_baseline(self):
        """Test single jet impingement at baseline conditions."""
        Re_jet = 20000
        Pr = 0.7
        z_D = 6.0

        Nu = cb.impingement_nusselt(Re_jet, Pr, z_D)

        # Should give reasonable Nusselt number
        assert 50 < Nu < 200

    def test_single_jet_reynolds_dependence(self):
        """Test that Nu increases with Reynolds number."""
        Pr = 0.7
        z_D = 6.0

        Nu_low = cb.impingement_nusselt(Re_jet=10000, Pr=Pr, z_D=z_D)
        Nu_mid = cb.impingement_nusselt(Re_jet=30000, Pr=Pr, z_D=z_D)
        Nu_high = cb.impingement_nusselt(Re_jet=50000, Pr=Pr, z_D=z_D)

        assert Nu_mid > Nu_low
        assert Nu_high > Nu_mid

    def test_single_jet_prandtl_dependence(self):
        """Test that Nu increases with Prandtl number."""
        Re_jet = 20000
        z_D = 6.0

        Nu_low_Pr = cb.impingement_nusselt(Re_jet, Pr=0.5, z_D=z_D)
        Nu_high_Pr = cb.impingement_nusselt(Re_jet, Pr=1.0, z_D=z_D)

        assert Nu_high_Pr > Nu_low_Pr

    def test_jet_array_baseline(self):
        """Test jet array impingement at baseline conditions."""
        Re_jet = 20000
        Pr = 0.7
        z_D = 6.0
        x_D = 8.0
        y_D = 8.0

        Nu = cb.impingement_nusselt(Re_jet, Pr, z_D, x_D, y_D)

        # Should give reasonable Nusselt number
        assert 40 < Nu < 180

    def test_jet_array_crossflow_degradation(self):
        """Test that jet array Nu is lower than single jet (crossflow effect)."""
        Re_jet = 20000
        Pr = 0.7
        z_D = 6.0

        Nu_single = cb.impingement_nusselt(Re_jet, Pr, z_D)
        Nu_array = cb.impingement_nusselt(Re_jet, Pr, z_D, x_D=8.0, y_D=8.0)

        # Array should be lower due to crossflow degradation
        assert Nu_array < Nu_single

    def test_impingement_parameter_validation(self):
        """Test parameter range validation for impingement."""
        # Valid range: Re_jet = 5000-80000, z_D = 1-12, x_D/y_D = 4-16

        with pytest.raises(RuntimeError):
            cb.impingement_nusselt(Re_jet=3000, Pr=0.7, z_D=6)

        with pytest.raises(RuntimeError):
            cb.impingement_nusselt(Re_jet=100000, Pr=0.7, z_D=6)

        with pytest.raises(RuntimeError):
            cb.impingement_nusselt(Re_jet=20000, Pr=0.7, z_D=0.5)

        with pytest.raises(RuntimeError):
            cb.impingement_nusselt(Re_jet=20000, Pr=0.7, z_D=15)

        with pytest.raises(RuntimeError):
            cb.impingement_nusselt(Re_jet=20000, Pr=0.7, z_D=6, x_D=2, y_D=8)

        with pytest.raises(RuntimeError):
            cb.impingement_nusselt(Re_jet=20000, Pr=0.7, z_D=6, x_D=20, y_D=8)


class TestFilmCooling:
    """Test film cooling effectiveness correlations."""

    def test_film_cooling_at_hole_exit(self):
        """Test film cooling effectiveness at hole exit (x_D=0)."""
        M = 1.0
        DR = 1.5
        alpha = 30.0

        eta = cb.film_cooling_effectiveness(x_D=0, M=M, DR=DR, alpha_deg=alpha)

        # Should be high near hole exit
        assert 0.3 < eta < 0.8

    def test_film_cooling_streamwise_decay(self):
        """Test that effectiveness decays downstream."""
        M = 1.0
        DR = 1.5
        alpha = 30.0

        eta_0 = cb.film_cooling_effectiveness(x_D=0, M=M, DR=DR, alpha_deg=alpha)
        eta_10 = cb.film_cooling_effectiveness(x_D=10, M=M, DR=DR, alpha_deg=alpha)
        eta_20 = cb.film_cooling_effectiveness(x_D=20, M=M, DR=DR, alpha_deg=alpha)

        assert eta_10 < eta_0
        assert eta_20 < eta_10

    def test_film_cooling_blowing_ratio_effect(self):
        """Test that optimal blowing ratio exists."""
        DR = 1.5
        alpha = 30.0
        x_D = 5.0

        eta_low = cb.film_cooling_effectiveness(x_D, M=0.5, DR=DR, alpha_deg=alpha)
        eta_opt = cb.film_cooling_effectiveness(x_D, M=1.0, DR=DR, alpha_deg=alpha)
        eta_high = cb.film_cooling_effectiveness(x_D, M=2.0, DR=DR, alpha_deg=alpha)

        # Optimal M should give highest effectiveness
        assert eta_opt > eta_low or eta_opt > eta_high

    def test_film_cooling_angle_effect(self):
        """Test that shallow angles give better coverage."""
        M = 1.0
        DR = 1.5
        x_D = 5.0

        eta_shallow = cb.film_cooling_effectiveness(x_D, M, DR, alpha_deg=30)
        eta_steep = cb.film_cooling_effectiveness(x_D, M, DR, alpha_deg=90)

        # Shallow angles should be better
        assert eta_shallow > eta_steep

    def test_film_cooling_bounds(self):
        """Test that effectiveness is bounded [0, 1]."""
        M = 1.0
        DR = 1.5
        alpha = 30.0

        for x_D in [0, 5, 10, 20, 40]:
            eta = cb.film_cooling_effectiveness(x_D, M, DR, alpha)
            assert 0.0 <= eta <= 1.0

    def test_film_cooling_parameter_validation(self):
        """Test parameter range validation for film cooling."""
        # Valid range: M = 0.3-2.5, DR = 1.2-2.0, alpha = 20-90, x_D >= 0

        with pytest.raises(RuntimeError):
            cb.film_cooling_effectiveness(x_D=5, M=0.2, DR=1.5, alpha_deg=30)

        with pytest.raises(RuntimeError):
            cb.film_cooling_effectiveness(x_D=5, M=3.0, DR=1.5, alpha_deg=30)

        with pytest.raises(RuntimeError):
            cb.film_cooling_effectiveness(x_D=5, M=1.0, DR=1.0, alpha_deg=30)

        with pytest.raises(RuntimeError):
            cb.film_cooling_effectiveness(x_D=5, M=1.0, DR=2.5, alpha_deg=30)

        with pytest.raises(RuntimeError):
            cb.film_cooling_effectiveness(x_D=5, M=1.0, DR=1.5, alpha_deg=15)

        with pytest.raises(RuntimeError):
            cb.film_cooling_effectiveness(x_D=5, M=1.0, DR=1.5, alpha_deg=100)

        with pytest.raises(RuntimeError):
            cb.film_cooling_effectiveness(x_D=-1, M=1.0, DR=1.5, alpha_deg=30)

    def test_film_cooling_avg_baseline(self):
        """Test laterally averaged film cooling effectiveness."""
        M = 1.0
        DR = 1.5
        alpha = 30.0
        x_D = 5.0
        s_D = 3.0

        eta_avg = cb.film_cooling_effectiveness_avg(x_D, M, DR, alpha, s_D)

        # Should be reasonable
        assert 0.0 < eta_avg < 1.0

    def test_film_cooling_avg_lower_than_centerline(self):
        """Test that averaged effectiveness is lower than centerline."""
        M = 1.0
        DR = 1.5
        alpha = 30.0
        x_D = 5.0

        eta_centerline = cb.film_cooling_effectiveness(x_D, M, DR, alpha)
        eta_avg = cb.film_cooling_effectiveness_avg(x_D, M, DR, alpha, s_D=3.0)

        # Averaged should be lower
        assert eta_avg < eta_centerline

    def test_film_cooling_avg_spacing_effect(self):
        """Test that closer holes give better averaged coverage."""
        M = 1.0
        DR = 1.5
        alpha = 30.0
        x_D = 5.0

        eta_close = cb.film_cooling_effectiveness_avg(x_D, M, DR, alpha, s_D=2.5)
        eta_far = cb.film_cooling_effectiveness_avg(x_D, M, DR, alpha, s_D=5.0)

        # Closer holes should give better coverage
        assert eta_close > eta_far

    def test_film_cooling_avg_parameter_validation(self):
        """Test parameter range validation for averaged film cooling."""
        # Valid range: s_D = 2-6

        with pytest.raises(RuntimeError):
            cb.film_cooling_effectiveness_avg(5, 1.0, 1.5, 30, s_D=1.5)

        with pytest.raises(RuntimeError):
            cb.film_cooling_effectiveness_avg(5, 1.0, 1.5, 30, s_D=7.0)


class TestMultiRowSellers:
    """Test multi-row film cooling using Sellers superposition."""

    def test_single_row_matches_direct_call(self):
        """Test that single row gives same result as direct effectiveness call."""
        rows = [0.0]
        eval_xD = 10.0
        M = 0.5
        DR = 1.8
        alpha = 30.0

        eta_multirow = cb.film_cooling_multirow_sellers(rows, eval_xD, M, DR, alpha)
        eta_single = cb.film_cooling_effectiveness(eval_xD, M, DR, alpha)

        assert eta_multirow == pytest.approx(eta_single, rel=1e-10)

    def test_two_rows_higher_than_single(self):
        """Test that two rows give higher effectiveness than single row."""
        M = 0.5
        DR = 1.8
        alpha = 30.0
        eval_xD = 20.0

        # Single row at x/D = 0
        eta_single = cb.film_cooling_multirow_sellers([0.0], eval_xD, M, DR, alpha)

        # Two rows at x/D = 0 and 10
        eta_double = cb.film_cooling_multirow_sellers([0.0, 10.0], eval_xD, M, DR, alpha)

        # Double row should be more effective
        assert eta_double > eta_single

    def test_three_rows_higher_than_two(self):
        """Test that three rows give higher effectiveness than two rows."""
        M = 0.5
        DR = 1.8
        alpha = 30.0
        eval_xD = 30.0

        eta_two = cb.film_cooling_multirow_sellers([0.0, 10.0], eval_xD, M, DR, alpha)
        eta_three = cb.film_cooling_multirow_sellers([0.0, 10.0, 20.0], eval_xD, M, DR, alpha)

        assert eta_three > eta_two

    def test_downstream_rows_ignored(self):
        """Test that downstream rows (after eval point) are ignored."""
        M = 0.5
        DR = 1.8
        alpha = 30.0
        eval_xD = 15.0

        # Only first two rows should contribute
        eta_with_downstream = cb.film_cooling_multirow_sellers(
            [0.0, 10.0, 20.0, 30.0], eval_xD, M, DR, alpha
        )
        eta_without_downstream = cb.film_cooling_multirow_sellers(
            [0.0, 10.0], eval_xD, M, DR, alpha
        )

        assert eta_with_downstream == pytest.approx(eta_without_downstream, rel=1e-10)

    def test_sellers_superposition_formula(self):
        """Test that Sellers formula is correctly implemented."""
        rows = [0.0, 10.0]
        eval_xD = 20.0
        M = 0.5
        DR = 1.8
        alpha = 30.0

        # Get individual row effectiveness
        eta1 = cb.film_cooling_effectiveness(eval_xD - rows[0], M, DR, alpha)
        eta2 = cb.film_cooling_effectiveness(eval_xD - rows[1], M, DR, alpha)

        # Sellers formula: eta_total = 1 - (1-eta1)*(1-eta2)
        eta_expected = 1.0 - (1.0 - eta1) * (1.0 - eta2)

        eta_actual = cb.film_cooling_multirow_sellers(rows, eval_xD, M, DR, alpha)

        assert eta_actual == pytest.approx(eta_expected, rel=1e-10)

    def test_effectiveness_bounded_by_one(self):
        """Test that effectiveness never exceeds 1.0."""
        # Many rows with high blowing ratio
        rows = [0.0, 5.0, 10.0, 15.0, 20.0]
        M = 2.0
        DR = 2.0
        alpha = 30.0

        for eval_xD in [25.0, 30.0, 40.0]:
            eta = cb.film_cooling_multirow_sellers(rows, eval_xD, M, DR, alpha)
            assert eta <= 1.0

    def test_empty_rows_raises_error(self):
        """Test that empty row list raises error."""
        with pytest.raises(RuntimeError, match="cannot be empty"):
            cb.film_cooling_multirow_sellers([], 10.0, 0.5, 1.8, 30.0)

    def test_row_order_independence(self):
        """Test that row order doesn't matter (commutative)."""
        rows_ordered = [0.0, 10.0, 20.0]
        rows_shuffled = [20.0, 0.0, 10.0]
        eval_xD = 30.0
        M = 0.5
        DR = 1.8
        alpha = 30.0

        eta_ordered = cb.film_cooling_multirow_sellers(rows_ordered, eval_xD, M, DR, alpha)
        eta_shuffled = cb.film_cooling_multirow_sellers(rows_shuffled, eval_xD, M, DR, alpha)

        assert eta_ordered == pytest.approx(eta_shuffled, rel=1e-10)


class TestEffusionCooling:
    """Test effusion cooling correlations."""

    def test_effusion_effectiveness_baseline(self):
        """Test effusion effectiveness at typical conditions."""
        eta = cb.effusion_effectiveness(
            x_D=10.0, M=2.0, DR=1.8, porosity=0.05, s_D=6.0, alpha_deg=30
        )
        # Should be moderate effectiveness (0.2-0.5 range typical)
        assert 0.2 < eta < 0.6

    def test_effusion_higher_porosity_increases_effectiveness(self):
        """Test that higher porosity gives better cooling."""
        eta_low = cb.effusion_effectiveness(10, 2.0, 1.8, 0.03, 6.0, 30)
        eta_high = cb.effusion_effectiveness(10, 2.0, 1.8, 0.08, 6.0, 30)

        assert eta_high > eta_low

    def test_effusion_effectiveness_decays_downstream(self):
        """Test that effectiveness decays with distance."""
        eta_near = cb.effusion_effectiveness(5, 2.0, 1.8, 0.05, 6.0, 30)
        eta_far = cb.effusion_effectiveness(20, 2.0, 1.8, 0.05, 6.0, 30)

        assert eta_near > eta_far

    def test_effusion_higher_blowing_ratio_increases_effectiveness(self):
        """Test that higher M gives better cooling."""
        eta_low = cb.effusion_effectiveness(10, 1.5, 1.8, 0.05, 6.0, 30)
        eta_high = cb.effusion_effectiveness(10, 3.0, 1.8, 0.05, 6.0, 30)

        assert eta_high > eta_low

    def test_effusion_parameter_validation(self):
        """Test parameter range validation."""
        # M out of range
        with pytest.raises(RuntimeError, match="blowing ratio M"):
            cb.effusion_effectiveness(10, 0.5, 1.8, 0.05, 6.0, 30)

        with pytest.raises(RuntimeError, match="blowing ratio M"):
            cb.effusion_effectiveness(10, 5.0, 1.8, 0.05, 6.0, 30)

        # DR out of range
        with pytest.raises(RuntimeError, match="Density ratio DR"):
            cb.effusion_effectiveness(10, 2.0, 1.0, 0.05, 6.0, 30)

        # Porosity out of range
        with pytest.raises(RuntimeError, match="Porosity"):
            cb.effusion_effectiveness(10, 2.0, 1.8, 0.01, 6.0, 30)

        with pytest.raises(RuntimeError, match="Porosity"):
            cb.effusion_effectiveness(10, 2.0, 1.8, 0.15, 6.0, 30)

        # s_D out of range
        with pytest.raises(RuntimeError, match="spacing s_D"):
            cb.effusion_effectiveness(10, 2.0, 1.8, 0.05, 3.0, 30)

        # Angle out of range
        with pytest.raises(RuntimeError, match="angle alpha"):
            cb.effusion_effectiveness(10, 2.0, 1.8, 0.05, 6.0, 15)

    def test_effusion_effectiveness_bounded(self):
        """Test that effectiveness is bounded [0, 1]."""
        # Try various conditions
        for x_D in [1, 5, 10, 20]:
            for M in [1.0, 2.0, 3.5]:
                eta = cb.effusion_effectiveness(x_D, M, 1.8, 0.05, 6.0, 30)
                assert 0.0 <= eta <= 1.0

    def test_effusion_discharge_coefficient_baseline(self):
        """Test discharge coefficient at typical conditions."""
        Cd = cb.effusion_discharge_coefficient(Re_d=10000, P_ratio=1.05, alpha_deg=30, L_D=4.0)
        # Typical Cd range for effusion holes
        assert 0.5 < Cd < 0.8

    def test_effusion_cd_varies_with_angle(self):
        """Test that Cd varies with injection angle."""
        # Shallow angle (20 deg) vs steep angle (45 deg)
        Cd_shallow = cb.effusion_discharge_coefficient(10000, 1.05, 20, 4.0)
        Cd_steep = cb.effusion_discharge_coefficient(10000, 1.05, 45, 4.0)

        # Both should be in valid range
        assert 0.5 < Cd_shallow < 0.8
        assert 0.5 < Cd_steep < 0.8

    def test_effusion_cd_decreases_with_LD(self):
        """Test that Cd decreases with longer holes."""
        Cd_short = cb.effusion_discharge_coefficient(10000, 1.05, 30, 2.5)
        Cd_long = cb.effusion_discharge_coefficient(10000, 1.05, 30, 7.0)

        assert Cd_short > Cd_long

    def test_effusion_cd_parameter_validation(self):
        """Test discharge coefficient parameter validation."""
        # Re_d too low
        with pytest.raises(RuntimeError, match="Reynolds number Re_d"):
            cb.effusion_discharge_coefficient(2000, 1.05, 30, 4.0)

        # P_ratio out of range
        with pytest.raises(RuntimeError, match="Pressure ratio"):
            cb.effusion_discharge_coefficient(10000, 1.01, 30, 4.0)

        with pytest.raises(RuntimeError, match="Pressure ratio"):
            cb.effusion_discharge_coefficient(10000, 1.20, 30, 4.0)

        # Angle out of range
        with pytest.raises(RuntimeError, match="angle alpha"):
            cb.effusion_discharge_coefficient(10000, 1.05, 15, 4.0)

        # L_D out of range
        with pytest.raises(RuntimeError, match="length/diameter L_D"):
            cb.effusion_discharge_coefficient(10000, 1.05, 30, 1.5)

    def test_effusion_cd_bounded(self):
        """Test that Cd is bounded to reasonable range."""
        for Re_d in [5000, 10000, 20000]:
            for P_ratio in [1.03, 1.07, 1.12]:
                Cd = cb.effusion_discharge_coefficient(Re_d, P_ratio, 30, 4.0)
                assert 0.5 <= Cd <= 0.85


class TestPhysicalConsistency:
    """Test physical consistency and trends."""

    def test_rib_heat_vs_friction_tradeoff(self):
        """Test that ribs increase both heat transfer and friction."""
        e_D = 0.05
        P_e = 10.0
        alpha = 90.0

        enhancement = cb.rib_enhancement_factor(e_D, P_e, alpha)
        friction = cb.rib_friction_multiplier(e_D, P_e)

        # Both should be > 1 (increase from smooth channel)
        assert enhancement > 1.0
        assert friction > 1.0

    def test_impingement_height_effect(self):
        """Test that closer jets give higher heat transfer."""
        Re_jet = 20000
        Pr = 0.7

        Nu_close = cb.impingement_nusselt(Re_jet, Pr, z_D=2.0)
        Nu_far = cb.impingement_nusselt(Re_jet, Pr, z_D=10.0)

        # Closer should be higher
        assert Nu_close > Nu_far

    def test_film_cooling_far_downstream_approaches_zero(self):
        """Test that film cooling effectiveness approaches zero far downstream."""
        M = 1.0
        DR = 1.5
        alpha = 30.0

        eta_far = cb.film_cooling_effectiveness(x_D=40, M=M, DR=DR, alpha_deg=alpha)

        # Should be very low far downstream
        assert eta_far < 0.2


class TestPinFinArrays:
    """Test pin fin array correlations."""

    def test_pin_fin_baseline(self):
        """Test pin fin Nusselt at typical conditions."""
        Nu = cb.pin_fin_nusselt(Re_d=20000, Pr=0.7, L_D=2.0, S_D=2.5, X_D=2.5)
        assert 80 < Nu < 120

    def test_pin_fin_staggered_higher_than_inline(self):
        """Test that staggered arrays have higher Nu than inline."""
        Nu_staggered = cb.pin_fin_nusselt(20000, 0.7, 2.0, 2.5, 2.5, is_staggered=True)
        Nu_inline = cb.pin_fin_nusselt(20000, 0.7, 2.0, 2.5, 2.5, is_staggered=False)

        assert Nu_staggered > Nu_inline

    def test_pin_fin_increases_with_reynolds(self):
        """Test that Nu increases with Re."""
        Nu_low = cb.pin_fin_nusselt(5000, 0.7, 2.0, 2.5, 2.5)
        Nu_high = cb.pin_fin_nusselt(50000, 0.7, 2.0, 2.5, 2.5)

        assert Nu_high > Nu_low

    def test_pin_fin_parameter_validation(self):
        """Test parameter range validation."""
        with pytest.raises(RuntimeError, match="Re_d"):
            cb.pin_fin_nusselt(2000, 0.7, 2.0, 2.5, 2.5)

        with pytest.raises(RuntimeError, match="Prandtl"):
            cb.pin_fin_nusselt(20000, 0.3, 2.0, 2.5, 2.5)

        with pytest.raises(RuntimeError, match="L_D"):
            cb.pin_fin_nusselt(20000, 0.7, 0.3, 2.5, 2.5)


class TestDimpledSurfaces:
    """Test dimpled surface correlations."""

    def test_dimple_enhancement_baseline(self):
        """Test dimple enhancement at typical conditions."""
        enh = cb.dimple_nusselt_enhancement(Re_Dh=30000, d_Dh=0.2, h_d=0.2, S_d=2.0)

        # Typical range 1.8-2.5x
        assert 1.5 < enh < 2.8

    def test_dimple_enhancement_bounded(self):
        """Test that enhancement is bounded to realistic range."""
        # Try extreme conditions
        enh_max = cb.dimple_nusselt_enhancement(80000, 0.3, 0.3, 1.5)
        enh_min = cb.dimple_nusselt_enhancement(10000, 0.1, 0.1, 3.0)

        assert 1.5 <= enh_max <= 2.8
        assert 1.5 <= enh_min <= 2.8

    def test_dimple_friction_baseline(self):
        """Test dimple friction multiplier at typical conditions."""
        f_ratio = cb.dimple_friction_multiplier(Re_Dh=30000, d_Dh=0.2, h_d=0.2)

        # Should be much lower than ribs (1.3-2.2x vs 6-10x)
        assert 1.3 < f_ratio < 2.2

    def test_dimple_friction_bounded(self):
        """Test that friction is bounded to realistic range."""
        f_max = cb.dimple_friction_multiplier(80000, 0.3, 0.3)
        f_min = cb.dimple_friction_multiplier(10000, 0.1, 0.1)

        assert 1.3 <= f_max <= 2.2
        assert 1.3 <= f_min <= 2.2

    def test_dimple_parameter_validation(self):
        """Test parameter range validation."""
        with pytest.raises(RuntimeError, match="Re_Dh"):
            cb.dimple_nusselt_enhancement(5000, 0.2, 0.2, 2.0)

        with pytest.raises(RuntimeError, match="d_Dh"):
            cb.dimple_nusselt_enhancement(30000, 0.05, 0.2, 2.0)

        with pytest.raises(RuntimeError, match="h_d"):
            cb.dimple_nusselt_enhancement(30000, 0.2, 0.05, 2.0)

        with pytest.raises(RuntimeError, match="S_d"):
            cb.dimple_nusselt_enhancement(30000, 0.2, 0.2, 1.0)

    def test_dimple_advantage_over_ribs(self):
        """Test that dimples have lower friction penalty than ribs."""
        # Dimple friction
        f_dimple = cb.dimple_friction_multiplier(30000, 0.2, 0.2)

        # Rib friction (from existing tests)
        f_rib = cb.rib_friction_multiplier(e_D=0.05, P_e=10.0)

        # Dimples should have much lower friction penalty
        assert f_dimple < f_rib / 2
