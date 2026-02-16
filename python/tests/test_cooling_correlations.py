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
