"""
Analytical unit tests for the Vatistas n-vortex model.

Tests cover:
  - Physical limits (r=0, r>>r_c, r<<r_c, peak location)
  - Closed-form specialisations (n=1, n=2)
  - Jacobian FD verification (100% coverage, rtol=1e-6)
  - Input validation (ValueError for n<1, Gamma<=0, r_c<=0)

No digitised data required here; see test_vatistas_validation.py for that.
"""

import math

import numpy as np
import pytest

import combaero as cb


def central_diff(func, x, h=None):
    if h is None:
        h = max(1e-6, abs(x) * 1e-6)
    return (func(x + h) - func(x - h)) / (2.0 * h)


GAMMA = 0.5  # m^2/s
R_C = 0.05  # m
RHO = 1.2  # kg/m^3


# ---------------------------------------------------------------------------
# Shape function -- physical limits
# ---------------------------------------------------------------------------


class TestV0Bar:
    def test_zero_at_origin(self):
        for n in [1, 2, 3, 4]:
            assert cb.vatistas_v0_bar(0.0, n) == 0.0

    @pytest.mark.parametrize("n", [1, 2, 3, 4])
    def test_peak_at_r_bar_one(self, n):
        expected_peak = 2.0 ** (-1.0 / n)
        assert math.isclose(cb.vatistas_v0_bar(1.0, n), expected_peak, rel_tol=1e-12)

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_free_vortex_limit(self, n):
        r_bar = 1000.0
        v = cb.vatistas_v0_bar(r_bar, n)
        free = 1.0 / r_bar
        assert math.isclose(v, free, rel_tol=1e-3)

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_solid_body_limit(self, n):
        r_bar = 0.001
        v = cb.vatistas_v0_bar(r_bar, n)
        solid = r_bar
        assert math.isclose(v, solid, rel_tol=1e-3)

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_dv0_drbar_positive_inside_core(self, n):
        assert cb.vatistas_dv0_bar_drbar(0.5, n) > 0.0

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_dv0_drbar_zero_at_peak(self, n):
        assert abs(cb.vatistas_dv0_bar_drbar(1.0, n)) < 1e-12

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_dv0_drbar_negative_outside_core(self, n):
        assert cb.vatistas_dv0_bar_drbar(2.0, n) < 0.0


# ---------------------------------------------------------------------------
# Pressure integral -- closed forms
# ---------------------------------------------------------------------------


class TestPressureIntegral:
    @pytest.mark.parametrize("r_bar", [0.5, 1.0, 2.0, 5.0, 10.0])
    def test_n1_closed_form(self, r_bar):
        expected = r_bar**2 / (2.0 * (1.0 + r_bar**2))
        assert math.isclose(cb.vatistas_pressure_integral(r_bar, 1.0), expected, rel_tol=1e-12)

    @pytest.mark.parametrize("r_bar", [0.5, 1.0, 2.0, 5.0, 10.0])
    def test_n2_closed_form(self, r_bar):
        expected = math.atan(r_bar**2) / 2.0
        assert math.isclose(cb.vatistas_pressure_integral(r_bar, 2.0), expected, rel_tol=1e-12)

    def test_zero_at_origin(self):
        for n in [1, 2, 3]:
            assert cb.vatistas_pressure_integral(0.0, n) == 0.0

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_monotone_increasing(self, n):
        r_vals = [0.5, 1.0, 2.0, 5.0]
        vals = [cb.vatistas_pressure_integral(r, n) for r in r_vals]
        assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))

    @pytest.mark.parametrize("r_bar", [0.5, 1.0, 2.0])
    def test_n3_matches_n2_limit_roughness(self, r_bar):
        # n=3 and n=2 give different values; just check n=3 > 0 and < n=1 asymptote
        I3 = cb.vatistas_pressure_integral(r_bar, 3.0)
        assert I3 > 0.0
        assert I3 < 1.0  # integral is bounded

    def test_d_pressure_integral_drbar_analytical(self):
        # FTC: d(I)/d(r_bar) == V0_bar(r_bar)^2 / r_bar  (exact identity, bypasses
        # quadrature noise that would contaminate a central-difference check)
        for n in [1, 2, 3]:
            for r_bar in [0.5, 1.0, 2.0, 5.0]:
                analytic = cb.vatistas_d_pressure_integral_drbar(r_bar, n)
                ftc = cb.vatistas_v0_bar(r_bar, n) ** 2 / r_bar
                assert math.isclose(analytic, ftc, rel_tol=1e-12), (
                    f"d_pressure_integral FTC mismatch at n={n}, r_bar={r_bar}: "
                    f"analytic={analytic}, ftc={ftc}"
                )


# ---------------------------------------------------------------------------
# Radial velocity
# ---------------------------------------------------------------------------


class TestVrBar:
    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_negative_sign(self, n):
        # Inward flow: Vr_bar < 0 everywhere (except r_bar=0)
        for r_bar in [0.5, 1.0, 2.0, 5.0]:
            assert cb.vatistas_vr_bar(r_bar, n) < 0.0

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_dvr_bar_drbar_fd(self, n):
        for r_bar in [0.5, 1.5, 3.0]:
            analytic = cb.vatistas_dvr_bar_drbar(r_bar, n)
            fd = central_diff(lambda rb, _n=n: cb.vatistas_vr_bar(rb, _n), r_bar)
            assert math.isclose(analytic, fd, rel_tol=1e-6), (
                f"dvr_bar FD mismatch at n={n}, r_bar={r_bar}: analytic={analytic}, fd={fd}"
            )


# ---------------------------------------------------------------------------
# Dimensional functions and Jacobians
# ---------------------------------------------------------------------------


class TestDimensional:
    def test_v_theta_peak_at_rc(self):
        for n in [1, 2, 3]:
            scale = GAMMA / (2.0 * math.pi * R_C)
            expected_peak = scale * 2.0 ** (-1.0 / n)
            assert math.isclose(
                cb.vatistas_v_theta(R_C, GAMMA, R_C, n), expected_peak, rel_tol=1e-12
            )

    def test_v_theta_zero_at_origin(self):
        assert cb.vatistas_v_theta(0.0, GAMMA, R_C, 2.0) == 0.0

    def test_delta_p_zero_at_origin(self):
        for n in [1, 2, 3]:
            assert cb.vatistas_delta_p(0.0, RHO, GAMMA, R_C, n) == 0.0

    def test_delta_p_monotone(self):
        r_vals = [R_C * 0.5, R_C, R_C * 2, R_C * 5]
        for n in [1, 2, 3]:
            vals = [cb.vatistas_delta_p(r, RHO, GAMMA, R_C, n) for r in r_vals]
            assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))

    @pytest.mark.parametrize("n", [1, 2, 3])
    @pytest.mark.parametrize("r_bar", [0.1, 0.5, 1.0, 2.0, 5.0])
    def test_dv_theta_dr_fd(self, n, r_bar):
        r = R_C * r_bar
        analytic_tuple = cb.vatistas_v_theta_and_jacobians(r, GAMMA, R_C, n)
        dV_dr_analytic = analytic_tuple[1]
        fd = central_diff(lambda ri, _n=n: cb.vatistas_v_theta(ri, GAMMA, R_C, _n), r)
        assert math.isclose(dV_dr_analytic, fd, rel_tol=1e-6, abs_tol=2e-8), (
            f"dV/dr mismatch n={n}, r_bar={r_bar}: {dV_dr_analytic} vs {fd}"
        )

    @pytest.mark.parametrize("n", [1, 2, 3])
    @pytest.mark.parametrize("r_bar", [0.5, 1.0, 2.0])
    def test_dv_theta_dGamma_fd(self, n, r_bar):
        r = R_C * r_bar
        analytic_tuple = cb.vatistas_v_theta_and_jacobians(r, GAMMA, R_C, n)
        dV_dGamma_analytic = analytic_tuple[2]
        fd = central_diff(lambda G, _n=n: cb.vatistas_v_theta(r, G, R_C, _n), GAMMA)
        assert math.isclose(dV_dGamma_analytic, fd, rel_tol=1e-6)

    @pytest.mark.parametrize("n", [1, 2, 3])
    @pytest.mark.parametrize("r_bar", [0.5, 1.0, 2.0])
    def test_dv_theta_drc_fd(self, n, r_bar):
        r = R_C * r_bar
        analytic_tuple = cb.vatistas_v_theta_and_jacobians(r, GAMMA, R_C, n)
        dV_drc_analytic = analytic_tuple[3]
        fd = central_diff(lambda rc, _n=n: cb.vatistas_v_theta(r, GAMMA, rc, _n), R_C)
        assert math.isclose(dV_drc_analytic, fd, rel_tol=1e-6)

    @pytest.mark.parametrize("n", [1, 2, 3])
    @pytest.mark.parametrize("r_bar", [0.5, 1.0, 2.0])
    def test_d_delta_p_dr_fd(self, n, r_bar):
        r = R_C * r_bar
        analytic_tuple = cb.vatistas_delta_p_and_jacobians(r, RHO, GAMMA, R_C, n)
        d_dP_dr_analytic = analytic_tuple[1]
        fd = central_diff(lambda ri, _n=n: cb.vatistas_delta_p(ri, RHO, GAMMA, R_C, _n), r)
        assert math.isclose(d_dP_dr_analytic, fd, rel_tol=1e-6), (
            f"d_delta_P/dr mismatch n={n}, r_bar={r_bar}: {d_dP_dr_analytic} vs {fd}"
        )

    @pytest.mark.parametrize("n", [1, 2, 3])
    @pytest.mark.parametrize("r_bar", [0.5, 1.0, 2.0])
    def test_d_delta_p_dGamma_fd(self, n, r_bar):
        r = R_C * r_bar
        analytic_tuple = cb.vatistas_delta_p_and_jacobians(r, RHO, GAMMA, R_C, n)
        d_dP_dGamma_analytic = analytic_tuple[2]
        fd = central_diff(lambda G, _n=n: cb.vatistas_delta_p(r, RHO, G, R_C, _n), GAMMA)
        assert math.isclose(d_dP_dGamma_analytic, fd, rel_tol=1e-6)

    @pytest.mark.parametrize("n", [1, 2, 3])
    @pytest.mark.parametrize("r_bar", [0.5, 1.0, 2.0])
    def test_d_delta_p_drc_fd(self, n, r_bar):
        r = R_C * r_bar
        analytic_tuple = cb.vatistas_delta_p_and_jacobians(r, RHO, GAMMA, R_C, n)
        d_dP_drc_analytic = analytic_tuple[3]
        fd = central_diff(lambda rc, _n=n: cb.vatistas_delta_p(r, RHO, GAMMA, rc, _n), R_C)
        assert math.isclose(d_dP_drc_analytic, fd, rel_tol=1e-6)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


class TestValidation:
    def test_n_below_one_raises(self):
        with pytest.raises((ValueError, RuntimeError)):
            cb.vatistas_v_theta(0.05, 1.0, 0.05, 0.5)

    def test_gamma_zero_raises(self):
        with pytest.raises((ValueError, RuntimeError)):
            cb.vatistas_v_theta(0.05, 0.0, 0.05, 2.0)

    def test_gamma_negative_raises(self):
        with pytest.raises((ValueError, RuntimeError)):
            cb.vatistas_v_theta(0.05, -1.0, 0.05, 2.0)

    def test_rc_zero_raises(self):
        with pytest.raises((ValueError, RuntimeError)):
            cb.vatistas_v_theta(0.05, 1.0, 0.0, 2.0)

    def test_vortex_class_n_below_one_raises(self):
        with pytest.raises(ValueError):
            cb.VatistasVortex(1.0, 0.05, n=0.5)

    def test_vortex_class_gamma_invalid_raises(self):
        with pytest.raises(ValueError):
            cb.VatistasVortex(-1.0, 0.05)

    def test_vortex_class_rc_invalid_raises(self):
        with pytest.raises(ValueError):
            cb.VatistasVortex(1.0, 0.0)


# ---------------------------------------------------------------------------
# VatistasVortex Python class
# ---------------------------------------------------------------------------


class TestVatistasVortexClass:
    def test_v_theta_vectorised(self):
        vortex = cb.VatistasVortex(GAMMA, R_C)
        r = np.array([0.0, R_C * 0.5, R_C, R_C * 2.0])
        v = vortex.V_theta(r)
        assert v.shape == (4,)
        assert v[0] == 0.0

    def test_v_theta_scalar(self):
        vortex = cb.VatistasVortex(GAMMA, R_C)
        v = vortex.V_theta(R_C)
        assert isinstance(v, float)

    def test_v_theta_max_matches_v_theta_at_rc(self):
        vortex = cb.VatistasVortex(GAMMA, R_C)
        assert math.isclose(vortex.V_theta_max(), float(vortex.V_theta(R_C)), rel_tol=1e-12)

    def test_delta_p_bar_n2_matches_eq11(self):
        vortex = cb.VatistasVortex(GAMMA, R_C, n=2.0)
        for r_bar in [0.5, 1.0, 2.0, 5.0, 10.0]:
            r = R_C * r_bar
            expected = (2.0 / math.pi) * math.atan(r_bar**2)
            # delta_P_bar is normalised by I_inf ~ pi/4, so
            # delta_P_bar = I(r_bar) / (pi/4) = arctan(r_bar^2)/2 / (pi/4)
            # = (2/pi)*arctan(r_bar^2)  -- matches paper Eq.11
            assert math.isclose(vortex.delta_P_bar(r), expected, rel_tol=1e-3), (
                f"delta_P_bar mismatch at r_bar={r_bar}: {vortex.delta_P_bar(r)} vs {expected}"
            )
