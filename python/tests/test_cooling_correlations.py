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


class TestThermalPerformanceFactor:
    """Validate Webb & Eckert (1972) thermal-hydraulic performance factor.

    Cross-checked against Singh & Ekkad tabulated eta values.
    """

    @pytest.mark.parametrize(
        "Nu_ratio, f_ratio, expected_eta, tol",
        [
            (2.45, 6.50, 1.31, 0.03),
            (2.15, 6.35, 1.16, 0.03),
            (1.95, 6.30, 1.05, 0.03),
            (1.75, 6.25, 0.94, 0.03),
            (1.55, 6.20, 0.84, 0.03),
        ],
    )
    def test_eta_vs_table(self, Nu_ratio, f_ratio, expected_eta, tol):
        """eta = Nu_ratio / f_ratio^(1/3) matches Singh & Ekkad table within 3%."""
        eta = cb.thermal_performance_factor(Nu_ratio, f_ratio)
        assert abs(eta - expected_eta) / expected_eta < tol

    def test_smooth_surface_gives_one(self):
        """Smooth surface (Nu_ratio=1, f_ratio=1) gives eta=1."""
        eta = cb.thermal_performance_factor(Nu_ratio=1.0, f_ratio=1.0)
        assert abs(eta - 1.0) < 1e-10

    def test_eta_greater_than_one_for_good_surface(self):
        """A surface with high Nu gain and modest friction gives eta > 1."""
        eta = cb.thermal_performance_factor(Nu_ratio=2.0, f_ratio=4.0)
        assert eta > 1.0

    def test_eta_less_than_one_for_bad_surface(self):
        """A surface with modest Nu gain and high friction gives eta < 1."""
        eta = cb.thermal_performance_factor(Nu_ratio=1.2, f_ratio=8.0)
        assert eta < 1.0


class TestCoolingWorkflowHelpers:
    """Test cooling workflow convenience helper functions."""

    def test_adiabatic_wall_temperature_matches_definition(self):
        """T_aw should satisfy eta definition exactly."""
        t_hot = 1800.0
        t_coolant = 800.0
        eta = 0.35

        t_aw = cb.adiabatic_wall_temperature(t_hot, t_coolant, eta)
        eta_back = (t_hot - t_aw) / (t_hot - t_coolant)

        assert t_coolant <= t_aw <= t_hot
        assert eta_back == pytest.approx(eta, rel=1e-12)

    def test_cooled_wall_heat_flux_reduces_with_effectiveness(self):
        """Film/effusion effectiveness should lower wall heat flux."""
        q_uncooled = cb.cooled_wall_heat_flux(
            T_hot=1850.0,
            T_coolant=760.0,
            h_hot=1100.0,
            h_coolant=600.0,
            eta=0.0,
            t_wall=0.0025,
            k_wall=18.0,
        )
        q_cooled = cb.cooled_wall_heat_flux(
            T_hot=1850.0,
            T_coolant=760.0,
            h_hot=1100.0,
            h_coolant=600.0,
            eta=0.45,
            t_wall=0.0025,
            k_wall=18.0,
        )

        assert q_uncooled > 0.0
        assert q_cooled > 0.0
        assert q_cooled < q_uncooled

    def test_cooling_workflow_helper_validation(self):
        """Helper functions should validate parameters and fail clearly."""
        with pytest.raises(RuntimeError, match="eta"):
            cb.adiabatic_wall_temperature(1700.0, 700.0, 1.2)

        with pytest.raises(RuntimeError, match="Heat transfer coefficients"):
            cb.cooled_wall_heat_flux(1700.0, 700.0, -10.0, 500.0, 0.2, 0.002, 20.0)
