"""
Tests for materials database and thermal conductivity functions.

Tests verify:
- Material thermal conductivity k(T) values against published data
- Temperature range validation
- TBC aging model behavior
- Database access functions
"""

import math

import pytest

import combaero as cb


class TestSuperalloys:
    """Test superalloy thermal conductivity functions."""

    def test_inconel718_ambient(self):
        """Test Inconel 718 at ambient temperature."""
        T = 300.0  # K (27degC)
        k = cb.k_inconel718(T)

        # At 300K (27degC): k ~= 11.4 + 0.0146*27 ~= 11.79 W/(m*K)
        assert k == pytest.approx(11.79, abs=0.1)

    def test_inconel718_high_temp(self):
        """Test Inconel 718 at high temperature."""
        T = 1000.0  # K (727degC)
        k = cb.k_inconel718(T)

        # At 1000K (727C): k approx 11.4 + 0.0146*727 approx 22.0 W/(m*K)
        assert k == pytest.approx(22.0, abs=0.5)

    def test_inconel718_temperature_range(self):
        """Test Inconel 718 valid temperature range."""
        # Valid range: 300-1200 K
        assert cb.k_inconel718(300) > 0
        assert cb.k_inconel718(1200) > 0

        # Out of range should raise error
        with pytest.raises(RuntimeError):
            cb.k_inconel718(250)

        with pytest.raises(RuntimeError):
            cb.k_inconel718(1300)

    def test_haynes230_ambient(self):
        """Test Haynes 230 at ambient temperature."""
        T = 300.0  # K
        k = cb.k_haynes230(T)

        # Brochure H-3000O: 8.9 W/m K at room temperature (#761)
        assert k == pytest.approx(8.94, abs=0.01)

    def test_haynes230_high_temp(self):
        """Test Haynes 230 at high temperature."""
        T = 1200.0  # K (927C)
        k = cb.k_haynes230(T)

        # Brochure H-3000O: 26.4 W/m K at 900 C, 28.4 at 1000 C (linear)
        assert k == pytest.approx(8.4 + 0.020 * (1200.0 - 273.15), rel=1e-12)

    def test_haynes230_temperature_range(self):
        """Test Haynes 230 valid temperature range."""
        # Valid range: 300-1273 K, where the brochure table ends (1000 C)
        assert cb.k_haynes230(300) > 0
        assert cb.k_haynes230(1000.0 + 273.15) == pytest.approx(28.4, rel=1e-12)

        with pytest.raises(RuntimeError):
            cb.k_haynes230(250)

        with pytest.raises(RuntimeError):
            cb.k_haynes230(1274.15)


class TestStructuralAlloys:
    """Test structural alloy thermal conductivity functions."""

    def test_ss316_ambient(self):
        """Test Stainless Steel 316 at ambient temperature."""
        T = 300.0  # K
        k = cb.k_stainless_steel_316(T)

        # At 300K: k approx 13.5 + 0.015*27 - 2.1e-6*27^2 approx 13.9 W/(m*K)
        assert k == pytest.approx(13.9, abs=0.2)

    def test_ss316_high_temp(self):
        """Test SS316 at high temperature."""
        T = 800.0  # K (527C)
        k = cb.k_stainless_steel_316(T)

        # Should be in reasonable range for austenitic stainless
        assert 18 < k < 25

    def test_ss316_temperature_range(self):
        """Test SS316 valid temperature range."""
        # Valid range: 300-1200 K
        assert cb.k_stainless_steel_316(300) > 0
        assert cb.k_stainless_steel_316(1200) > 0

        with pytest.raises(RuntimeError):
            cb.k_stainless_steel_316(250)

    def test_aluminum_6061_ambient(self):
        """Test Aluminum 6061 at ambient temperature."""
        T = 300.0  # K
        k = cb.k_aluminum_6061(T)

        # At 300K: k approx 167.0 - 0.012*27 approx 166.7 W/(m*K)
        assert k == pytest.approx(166.7, abs=1.0)

    def test_aluminum_6061_elevated_temp(self):
        """Test Al6061 at elevated temperature."""
        T = 500.0  # K (227C)
        k = cb.k_aluminum_6061(T)

        # Should decrease with temperature
        k_ambient = cb.k_aluminum_6061(300)
        assert k < k_ambient

    def test_aluminum_6061_temperature_range(self):
        """Test Al6061 valid temperature range (limited by T6 temper)."""
        # Valid range: 200-600 K
        assert cb.k_aluminum_6061(200) > 0
        assert cb.k_aluminum_6061(600) > 0

        with pytest.raises(RuntimeError):
            cb.k_aluminum_6061(150)

        with pytest.raises(RuntimeError):
            cb.k_aluminum_6061(700)


class TestThermalBarrierCoatings:
    """YSZ: Zhu & Miller 2000 (NASA/TM-2000-210237), plasma-sprayed (#761)."""

    def test_as_sprayed_is_about_one_and_flat_in_temperature(self):
        """Zhu & Miller: as-sprayed plasma coating 'about 1 W/m-K'. The old
        uncited fit rose 0.81 -> 1.39 from 300 to 1589 K."""
        for T in (300.0, 800.0, 1200.0, 1700.0):
            assert cb.k_tbc_ysz(T, hours=0) == pytest.approx(1.0, rel=1e-12)

    def test_the_sintering_is_zhu_miller_eq_1(self):
        """k = k0 + (k_inf - k0) A exp(-Q/RT)[1 - exp(-t/tau)], tau = tau0
        exp(Q_tau/RT) with t in s (Eq. 1, p. 7; k_inf = 1.74 inferred)."""
        R = 8.31446261815324
        for T, h in ((1263.15, 4.0), (1373.15, 12.0), (1593.15, 2.0), (1500.0, 100.0)):
            tau = 572.5 * math.exp(41710.0 / (R * T))
            prog = 102.2 * math.exp(-68228.0 / (R * T)) * (1.0 - math.exp(-h * 3600.0 / tau))
            assert cb.k_tbc_ysz(T, hours=h) == pytest.approx(1.0 + 0.74 * prog, rel=1e-12)

    def test_the_plateaus_match_figure_4(self):
        """Fig. 4 plateaus (read off the figure, unconfirmed): about 1.09 /
        1.17 / 1.44 W/m K at 990 / 1100 / 1320 C surface after 32 h."""
        for T_C, k_fig in ((990.0, 1.09), (1100.0, 1.17), (1320.0, 1.44)):
            assert cb.k_tbc_ysz(T_C + 273.15, hours=32) == pytest.approx(k_fig, abs=0.03)

    def test_sintering_rises_then_saturates(self):
        T = 1593.15
        k = [cb.k_tbc_ysz(T, hours=h) for h in (0.0, 1.0, 4.0, 8.0, 32.0, 1000.0)]
        assert all(b >= a for a, b in zip(k, k[1:], strict=False))
        assert k[-1] == pytest.approx(k[-2], rel=1e-3)  # saturated within a day
        assert k[-1] < 1.74  # below the fully sintered limit

    def test_hotter_sinters_more(self):
        assert cb.k_tbc_ysz(1500.0, hours=100) > cb.k_tbc_ysz(1300.0, hours=100)

    def test_cold_exposure_is_negligible_without_a_jump(self):
        """No 1073 K gate any more: the Arrhenius factor makes it small, and
        k is continuous through 1073 K."""
        assert cb.k_tbc_ysz(773.0, hours=1000) == pytest.approx(1.0, rel=5e-3)
        lo, hi = cb.k_tbc_ysz(1072.9, hours=100), cb.k_tbc_ysz(1073.1, hours=100)
        assert hi - lo < 1e-3

    def test_ebpvd_is_as_deposited_only(self):
        """EB-PVD keeps its as-deposited fit (uncited) and gets no sintering:
        Zhu & Miller's kinetics are for plasma-sprayed coatings."""
        T = 1300.0
        assert cb.k_tbc_ysz(T, hours=0, is_ebpvd=True) == pytest.approx(1.5 + 0.0002 * (T - 273.15))
        assert cb.k_tbc_ysz(T, hours=100, is_ebpvd=True) == cb.k_tbc_ysz(T, hours=0, is_ebpvd=True)

    def test_ysz_temperature_range(self):
        assert cb.k_tbc_ysz(300) > 0
        assert cb.k_tbc_ysz(1700) > 0
        with pytest.raises(RuntimeError):
            cb.k_tbc_ysz(250)
        with pytest.raises(RuntimeError):
            cb.k_tbc_ysz(1800)

    def test_ysz_negative_hours_treated_as_zero(self):
        assert cb.k_tbc_ysz(1000.0, hours=-100) == cb.k_tbc_ysz(1000.0, hours=0)


class TestMaterialDatabase:
    """Test that list_materials returns expected materials."""

    materials = cb.list_materials()

    # Should be a list
    assert isinstance(materials, list)

    # Should contain our known materials
    assert "inconel718" in materials
    assert "haynes230" in materials
    assert "ss316" in materials
    assert "al6061" in materials
    assert "ysz" in materials

    # Should have at least 5 materials
    assert len(materials) >= 5


class TestTemperatureTrends:
    """Test that thermal conductivity trends are physically reasonable."""

    def test_superalloys_increase_with_temperature(self):
        """Test that superalloy k increases with T (electron dominance)."""
        # Inconel 718
        k_300 = cb.k_inconel718(300)
        k_600 = cb.k_inconel718(600)
        k_900 = cb.k_inconel718(900)

        assert k_600 > k_300
        assert k_900 > k_600

    def test_aluminum_decreases_with_temperature(self):
        """Test that aluminum k decreases with T (phonon scattering)."""
        k_300 = cb.k_aluminum_6061(300)
        k_400 = cb.k_aluminum_6061(400)
        k_500 = cb.k_aluminum_6061(500)

        assert k_400 < k_300
        assert k_500 < k_400

    def test_tbc_as_sprayed_is_flat_in_temperature(self):
        """Plasma-sprayed YSZ as sprayed: about 1 W/m K at every T (Zhu &
        Miller 2000); the literature has it roughly flat, 0.8-1.1."""
        assert cb.k_tbc_ysz(500, hours=0) == cb.k_tbc_ysz(1500, hours=0)


class TestMagnitudeRanges:
    """Test that thermal conductivity values are in reasonable ranges."""

    def test_superalloy_magnitude(self):
        """Test superalloy k is in typical range (10-30 W/(m*K))."""
        T = 800.0

        k_inconel = cb.k_inconel718(T)
        k_haynes = cb.k_haynes230(T)

        assert 10 < k_inconel < 30
        assert 10 < k_haynes < 30

    def test_stainless_steel_magnitude(self):
        """Test stainless steel k is in typical range (13-25 W/(m*K))."""
        T = 600.0
        k = cb.k_stainless_steel_316(T)

        assert 13 < k < 25

    def test_aluminum_magnitude(self):
        """Test aluminum k is in typical range (150-180 W/(m*K))."""
        T = 300.0
        k = cb.k_aluminum_6061(T)

        assert 150 < k < 180

    def test_tbc_magnitude(self):
        """Test TBC k is in typical range (0.8-2.0 W/(m*K))."""
        T = 1200.0

        k_fresh = cb.k_tbc_ysz(T, hours=0)
        k_aged = cb.k_tbc_ysz(T, hours=10000)

        assert 0.8 < k_fresh < 2.0
        assert 0.8 < k_aged < 2.0
