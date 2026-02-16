"""
Tests for materials database and thermal conductivity functions.

Tests verify:
- Material thermal conductivity k(T) values against published data
- Temperature range validation
- TBC aging model behavior
- Database access functions
"""

import pytest

import combaero as cb


class TestSuperalloys:
    """Test superalloy thermal conductivity functions."""

    def test_inconel718_ambient(self):
        """Test Inconel 718 at ambient temperature."""
        T = 300.0  # K (27°C)
        k = cb.k_inconel718(T)

        # At 300K (27°C): k ≈ 11.4 + 0.0146*27 ≈ 11.79 W/(m·K)
        assert k == pytest.approx(11.79, abs=0.1)

    def test_inconel718_high_temp(self):
        """Test Inconel 718 at high temperature."""
        T = 1000.0  # K (727°C)
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

        # At 300K: k approx 11.8 + 0.0158*27 + 1.1e-6*27^2 approx 12.23 W/(m*K)
        assert k == pytest.approx(12.23, abs=0.1)

    def test_haynes230_high_temp(self):
        """Test Haynes 230 at high temperature."""
        T = 1200.0  # K (927C)
        k = cb.k_haynes230(T)

        # Quadratic fit should give reasonable value
        assert 20 < k < 35

    def test_haynes230_temperature_range(self):
        """Test Haynes 230 valid temperature range."""
        # Valid range: 300-1400 K
        assert cb.k_haynes230(300) > 0
        assert cb.k_haynes230(1400) > 0

        with pytest.raises(RuntimeError):
            cb.k_haynes230(250)

        with pytest.raises(RuntimeError):
            cb.k_haynes230(1500)


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
    """Test thermal barrier coating thermal conductivity with aging."""

    def test_ysz_as_sprayed_ambient(self):
        """Test YSZ as-sprayed at ambient temperature."""
        T = 300.0  # K
        k = cb.k_tbc_ysz(T, hours=0)

        # At 300K: k approx 0.8 + 0.00045*27 approx 0.81 W/(m*K)
        assert k == pytest.approx(0.81, abs=0.05)

    def test_ysz_as_sprayed_high_temp(self):
        """Test YSZ as-sprayed at high temperature."""
        T = 1500.0  # K (1227C)
        k = cb.k_tbc_ysz(T, hours=0)

        # At 1500K: k approx 0.8 + 0.00045*1227 approx 1.35 W/(m*K)
        assert k == pytest.approx(1.35, abs=0.1)

    def test_ysz_sintering_increases_conductivity(self):
        """Test that sintering increases thermal conductivity."""
        T = 1500.0  # K

        k_fresh = cb.k_tbc_ysz(T, hours=0)
        k_aged_100h = cb.k_tbc_ysz(T, hours=100)
        k_aged_1000h = cb.k_tbc_ysz(T, hours=1000)
        k_aged_10000h = cb.k_tbc_ysz(T, hours=10000)

        # Conductivity should increase with aging
        assert k_aged_100h > k_fresh
        assert k_aged_1000h > k_aged_100h
        assert k_aged_10000h > k_aged_1000h

        # Should approach fully sintered value (~1.85 W/(m*K))
        assert k_aged_10000h < 2.0

    def test_ysz_sintering_temperature_dependence(self):
        """Test that sintering rate increases with temperature."""
        hours = 1000

        # Higher temperature -> faster sintering -> higher k increase
        k_low_T = cb.k_tbc_ysz(T=1000, hours=hours)
        k_high_T = cb.k_tbc_ysz(T=1500, hours=hours)

        k_low_T_fresh = cb.k_tbc_ysz(T=1000, hours=0)
        k_high_T_fresh = cb.k_tbc_ysz(T=1500, hours=0)

        # Relative increase should be larger at high temperature
        increase_low_T = (k_low_T - k_low_T_fresh) / k_low_T_fresh
        increase_high_T = (k_high_T - k_high_T_fresh) / k_high_T_fresh

        assert increase_high_T > increase_low_T

    def test_ysz_temperature_range(self):
        """Test YSZ valid temperature range."""
        # Valid range: 300-1700 K
        assert cb.k_tbc_ysz(300) > 0
        assert cb.k_tbc_ysz(1700) > 0

        with pytest.raises(RuntimeError):
            cb.k_tbc_ysz(250)

        with pytest.raises(RuntimeError):
            cb.k_tbc_ysz(1800)

    def test_ysz_negative_hours_treated_as_zero(self):
        """Test that negative hours are treated as as-sprayed."""
        T = 1000.0

        k_zero = cb.k_tbc_ysz(T, hours=0)
        k_negative = cb.k_tbc_ysz(T, hours=-100)

        assert k_negative == k_zero


class TestTBCSintering:
    """Test TBC sintering model enhancements."""

    def test_tbc_no_sintering_at_low_temperature(self):
        """Verify no sintering occurs below 1073 K (800 C)."""
        # At 773 K (500 C), sintering should not occur
        k_cold_fresh = cb.k_tbc_ysz(T=773, hours=0)
        k_cold_aged = cb.k_tbc_ysz(T=773, hours=1000)

        # Should be identical (no sintering)
        assert k_cold_fresh == pytest.approx(k_cold_aged, rel=1e-10)

    def test_tbc_sintering_at_high_temperature(self):
        """Verify sintering occurs above 1073 K."""
        # At 1373 K (1100 C), sintering should be significant
        k_hot_fresh = cb.k_tbc_ysz(T=1373, hours=0)
        k_hot_aged = cb.k_tbc_ysz(T=1373, hours=200)

        # Aged coating should have higher k due to pore closure
        assert k_hot_aged > k_hot_fresh

        # Should be bounded by fully dense YSZ limit
        assert k_hot_aged <= 1.85

    def test_tbc_sintering_increases_with_time(self):
        """Verify k increases monotonically with exposure time."""
        T = 1500  # K
        k_0h = cb.k_tbc_ysz(T, hours=0)
        k_100h = cb.k_tbc_ysz(T, hours=100)
        k_500h = cb.k_tbc_ysz(T, hours=500)
        k_1000h = cb.k_tbc_ysz(T, hours=1000)

        # Should increase monotonically
        assert k_100h > k_0h
        assert k_500h > k_100h
        assert k_1000h > k_500h

    def test_tbc_ebpvd_higher_than_aps(self):
        """Verify EB-PVD has higher initial k than APS."""
        T = 800  # K - use moderate temperature

        # As-sprayed comparison
        k_aps = cb.k_tbc_ysz(T, hours=0, is_ebpvd=False)
        k_ebpvd = cb.k_tbc_ysz(T, hours=0, is_ebpvd=True)

        # EB-PVD should have higher initial k (columnar vs splat)
        assert k_ebpvd > k_aps

        # Typical values at 800K
        assert 0.8 < k_aps < 1.1  # APS range
        assert 1.5 < k_ebpvd < 1.7  # EB-PVD range

    def test_tbc_ebpvd_also_sinters(self):
        """Verify EB-PVD also undergoes sintering at high T."""
        T = 1500  # K

        k_ebpvd_fresh = cb.k_tbc_ysz(T, hours=0, is_ebpvd=True)
        k_ebpvd_aged = cb.k_tbc_ysz(T, hours=500, is_ebpvd=True)

        # EB-PVD should also sinter (though starts higher)
        assert k_ebpvd_aged > k_ebpvd_fresh

    def test_tbc_sintering_rate_increases_with_temperature(self):
        """Verify Arrhenius behavior - faster sintering at higher T."""
        # Compare sintering progress at two temperatures
        hours = 100

        k_low_fresh = cb.k_tbc_ysz(T=1273, hours=0)
        k_low_aged = cb.k_tbc_ysz(T=1273, hours=hours)
        delta_low = k_low_aged - k_low_fresh

        k_high_fresh = cb.k_tbc_ysz(T=1573, hours=0)
        k_high_aged = cb.k_tbc_ysz(T=1573, hours=hours)
        delta_high = k_high_aged - k_high_fresh

        # Higher temperature should cause more sintering in same time
        assert delta_high > delta_low

    def test_tbc_bounded_by_dense_limit(self):
        """Verify k never exceeds fully dense YSZ limit."""
        # Extreme aging at very high temperature
        k_extreme = cb.k_tbc_ysz(T=1600, hours=10000)

        # Should not exceed fully dense limit
        assert k_extreme <= 1.85

    def test_tbc_backward_compatibility(self):
        """Verify default behavior matches original APS implementation."""
        # Default should be APS (is_ebpvd=False)
        k_default = cb.k_tbc_ysz(T=1500, hours=100)
        k_aps_explicit = cb.k_tbc_ysz(T=1500, hours=100, is_ebpvd=False)

        assert k_default == pytest.approx(k_aps_explicit, rel=1e-10)


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

    def test_tbc_increases_with_temperature(self):
        """Test that TBC k increases with T (radiation contribution)."""
        k_500 = cb.k_tbc_ysz(500, hours=0)
        k_1000 = cb.k_tbc_ysz(1000, hours=0)
        k_1500 = cb.k_tbc_ysz(1500, hours=0)

        assert k_1000 > k_500
        assert k_1500 > k_1000


class TestMagnitudeRanges:
    """Test that thermal conductivity values are in reasonable ranges."""

    def test_superalloy_magnitude(self):
        """Test superalloy k is in typical range (10-30 W/(m·K))."""
        T = 800.0

        k_inconel = cb.k_inconel718(T)
        k_haynes = cb.k_haynes230(T)

        assert 10 < k_inconel < 30
        assert 10 < k_haynes < 30

    def test_stainless_steel_magnitude(self):
        """Test stainless steel k is in typical range (13-25 W/(m·K))."""
        T = 600.0
        k = cb.k_stainless_steel_316(T)

        assert 13 < k < 25

    def test_aluminum_magnitude(self):
        """Test aluminum k is in typical range (150-180 W/(m·K))."""
        T = 300.0
        k = cb.k_aluminum_6061(T)

        assert 150 < k < 180

    def test_tbc_magnitude(self):
        """Test TBC k is in typical range (0.8-2.0 W/(m·K))."""
        T = 1200.0

        k_fresh = cb.k_tbc_ysz(T, hours=0)
        k_aged = cb.k_tbc_ysz(T, hours=10000)

        assert 0.8 < k_fresh < 2.0
        assert 0.8 < k_aged < 2.0
