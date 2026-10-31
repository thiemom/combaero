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

    def test_haynes230_is_alpha_cp_rho_from_the_brochure(self):
        """Brochure H-3000O: k = alpha cp rho, 8970 kg/m3 (#881)."""
        # 500 C: 38.5e-3 cm2/s, 473 J/kg C.
        assert cb.k_haynes230(500.0 + 273.15) == pytest.approx(16.33, rel=1e-12)
        assert pytest.approx(16.33, abs=0.005) == 38.5e-7 * 473.0 * 8970.0
        # RT (25 C): 24.2e-3 cm2/s, 397 J/kg C.
        assert cb.k_haynes230(25.0 + 273.15) == pytest.approx(8.62, rel=1e-12)
        # 800 C row skipped (cp misprint): linear between 700 and 900 C.
        assert cb.k_haynes230(800.0 + 273.15) == pytest.approx(0.5 * (22.14 + 24.25), rel=1e-12)

    def test_haynes230_sits_inside_the_gen3csp_band(self):
        """Gen3 CSP measured k (95% band) at 500 C: 16.02 +- 1.01; the brochure's
        printed conductivity column (18.4) is outside it, the product is inside."""
        assert abs(cb.k_haynes230(500.0 + 273.15) - 16.02) < 1.01
        assert abs(18.4 - 16.02) > 1.01

    def test_haynes230_temperature_range(self):
        """Valid 25-1000 C, where the brochure table ends."""
        assert cb.k_haynes230(1000.0 + 273.15) == pytest.approx(26.68, rel=1e-12)

        with pytest.raises(RuntimeError):
            cb.k_haynes230(250)

        with pytest.raises(RuntimeError):
            cb.k_haynes230(1274.15)


class TestStructuralAlloys:
    """Test structural alloy thermal conductivity functions."""

    def test_ss316_is_the_gen3csp_table(self):
        """Georgia Tech Gen3 CSP measured k, 260-820 C (#881)."""
        assert cb.k_stainless_steel_316(260.0 + 273.15) == pytest.approx(17.83, rel=1e-12)
        assert cb.k_stainless_steel_316(500.0 + 273.15) == pytest.approx(20.78, rel=1e-12)
        assert cb.k_stainless_steel_316(820.0 + 273.15) == pytest.approx(23.12, rel=1e-12)
        assert cb.k_stainless_steel_316(520.0 + 273.15) == pytest.approx(
            0.5 * (20.78 + 21.16), rel=1e-12
        )

    def test_ss316_refuses_outside_the_measured_range(self):
        with pytest.raises(RuntimeError):
            cb.k_stainless_steel_316(300.0)
        with pytest.raises(RuntimeError):
            cb.k_stainless_steel_316(1200.0)
        # The database lookup holds the end values.
        assert cb.get_material_conductivity("ss316", 300.0) == pytest.approx(17.83, rel=1e-12)

    def test_ss310_is_the_gen3csp_table(self):
        """Georgia Tech Gen3 CSP measured k, 260-800 C (#881)."""
        assert cb.k_stainless_steel_310(260.0 + 273.15) == pytest.approx(17.13, rel=1e-12)
        assert cb.k_stainless_steel_310(500.0 + 273.15) == pytest.approx(20.55, rel=1e-12)
        assert cb.k_stainless_steel_310(800.0 + 273.15) == pytest.approx(25.12, rel=1e-12)
        assert cb.get_material_conductivity("ss310", 600.0 + 273.15) == pytest.approx(
            0.5 * (22.08 + 23.08), rel=1e-12
        )
        with pytest.raises(RuntimeError):
            cb.k_stainless_steel_310(300.0)
        with pytest.raises(RuntimeError):
            cb.k_stainless_steel_310(800.0 + 273.15 + 1.0)

    def test_aluminum_6061_is_deprecated_but_still_answers(self):
        """No cited source (#761): warns, removed in 0.9.0, unchanged until then."""
        with pytest.warns(DeprecationWarning, match="0.9.0"):
            assert cb.k_aluminum_6061(300.0) == pytest.approx(167.0 - 0.012 * 26.85, rel=1e-12)
        with pytest.warns(DeprecationWarning, match="al6061"):
            assert cb.get_material_conductivity("al6061", 300.0) > 0
        with pytest.warns(DeprecationWarning), pytest.raises(RuntimeError):
            cb.k_aluminum_6061(700.0)


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
    assert "ss310" in materials
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

    def test_tbc_magnitude(self):
        """Test TBC k is in typical range (0.8-2.0 W/(m*K))."""
        T = 1200.0

        k_fresh = cb.k_tbc_ysz(T, hours=0)
        k_aged = cb.k_tbc_ysz(T, hours=10000)

        assert 0.8 < k_fresh < 2.0
        assert 0.8 < k_aged < 2.0


class TestLinerAlloysAndBondCoat:
    """Manufacturer brochure tables (Haynes International) and the Zhu &
    Miller 1999 bond-coat line: exact at every printed point."""

    @pytest.mark.parametrize(
        ("func", "table"),
        [
            (
                "k_hastelloy_x",
                [
                    (25, 9.2),
                    (100, 11.2),
                    (200, 14.1),
                    (400, 16.9),
                    (600, 20.9),
                    (650, 21.9),
                    (700, 22.8),
                    (750, 23.8),
                    (800, 24.7),
                    (850, 25.7),
                    (900, 26.7),
                    (1000, 28.7),
                ],
            ),
            (
                "k_haynes188",
                [
                    (100, 12.2),
                    (200, 14.3),
                    (300, 15.9),
                    (400, 17.5),
                    (500, 19.3),
                    (600, 21.1),
                    (700, 23.0),
                    (800, 24.8),
                    (900, 25.5),
                    (1000, 27.6),
                ],
            ),
            (
                "k_haynes263",
                [
                    (100, 13.0),
                    (200, 14.7),
                    (300, 16.3),
                    (400, 18.0),
                    (500, 19.7),
                    (600, 21.4),
                    (700, 23.0),
                    (800, 24.7),
                    (900, 26.8),
                    (1000, 28.5),
                ],
            ),
        ],
    )
    def test_the_brochure_table_is_reproduced(self, func: str, table: list) -> None:
        f = getattr(cb, func)
        for T_C, k in table:
            if T_C + 273.15 >= 300.0:
                assert f(T_C + 273.15) == pytest.approx(k, rel=1e-12)

    def test_between_points_is_linear(self) -> None:
        # Hastelloy X between the data sheet's 400 C and the brochure's 600 C.
        assert cb.k_hastelloy_x(500.0 + 273.15) == pytest.approx(0.5 * (16.9 + 20.9), rel=1e-12)

    def test_the_ranges_are_the_tables(self) -> None:
        for f, T_max in (
            (cb.k_hastelloy_x, 1273.15),
            (cb.k_haynes188, 1273.15),
            (cb.k_haynes263, 1273.15),
        ):
            assert f(T_max) > 0
            with pytest.raises(RuntimeError):
                f(T_max + 1.0)

    def test_the_bond_coat_is_the_zhu_miller_line_over_its_drawn_range(self) -> None:
        for T_C in (575.0, 800.0, 1100.0):
            assert cb.k_bondcoat_nicocraly(T_C + 273.15) == pytest.approx(
                6.729572 + 0.02344 * T_C, rel=1e-12
            )
        with pytest.raises(RuntimeError):
            cb.k_bondcoat_nicocraly(500.0 + 273.15)
        with pytest.raises(RuntimeError):
            cb.k_bondcoat_nicocraly(1200.0 + 273.15)

    def test_they_are_in_the_material_database(self) -> None:
        names = cb.list_materials()
        for n in ("hastelloy_x", "haynes188", "haynes263", "nicocraly"):
            assert n in names
        assert cb.get_material_conductivity("hastelloy_x", 873.15) == pytest.approx(20.9, rel=1e-12)


class TestFurtherAlloysCopperAlumina:
    """Tables from Haynes, Special Metals, Virgamet, Edelstahlwerke
    Suedwestfalen and NSRDS-NBS 8 (#761): exact at every printed point, linear
    between them, refused outside the table."""

    TABLES = [
        (
            "k_haynes25",
            "C",
            298.15,
            1273.15,
            [
                (25, 10.5),
                (100, 12.0),
                (200, 14.0),
                (300, 15.9),
                (400, 17.7),
                (500, 19.5),
                (600, 21.2),
                (700, 22.9),
                (800, 24.5),
                (900, 26.0),
                (1000, 27.5),
            ],
        ),
        (
            "k_inconel625",
            "C",
            116.15,
            1255.15,
            [
                (-157, 7.2),
                (-129, 7.5),
                (-73, 8.4),
                (-18, 9.2),
                (21, 9.8),
                (38, 10.1),
                (93, 10.8),
                (204, 12.5),
                (316, 14.1),
                (427, 15.7),
                (538, 17.5),
                (649, 19.0),
                (760, 20.8),
                (871, 22.8),
                (982, 25.2),
            ],
        ),
        (
            "k_nimonic75",
            "C",
            293.15,
            1273.15,
            [
                (20, 11.7),
                (400, 18.6),
                (600, 22.7),
                (700, 24.7),
                (800, 26.5),
                (900, 28.4),
                (1000, 30.1),
            ],
        ),
        (
            "k_copper",
            "K",
            100.0,
            1200.0,
            [
                (100, 483),
                (150, 428),
                (200, 413),
                (250, 404),
                (273, 401),
                (300, 398),
                (350, 394),
                (400, 392),
                (500, 388),
                (600, 383),
                (700, 377),
                (800, 371),
                (900, 364),
                (1000, 357),
                (1100, 350),
                (1200, 342),
            ],
        ),
        (
            "k_alumina",
            "K",
            100.0,
            2000.0,
            [
                (100, 133.0),
                (150, 77.0),
                (200, 55.0),
                (250, 43.4),
                (273, 39.7),
                (300, 36.0),
                (350, 30.7),
                (400, 26.4),
                (500, 20.2),
                (600, 15.8),
                (700, 12.6),
                (800, 10.4),
                (900, 8.9),
                (1000, 7.85),
                (1100, 7.1),
                (1200, 6.55),
                (1300, 6.13),
                (1400, 5.85),
                (1500, 5.66),
                (1600, 5.56),
                (1700, 5.54),
                (1800, 5.59),
                (1900, 5.74),
                (2000, 6.0),
            ],
        ),
        ("k_x22crmov12_1", "C", 293.15, 923.15, [(20, 24.0), (650, 29.0)]),
    ]

    @pytest.mark.parametrize(("func", "unit", "T_min", "T_max", "table"), TABLES)
    def test_the_table_is_reproduced_and_bounded(
        self, func: str, unit: str, T_min: float, T_max: float, table: list
    ) -> None:
        f = getattr(cb, func)
        to_K = 273.15 if unit == "C" else 0.0
        for x, k in table:
            assert f(x + to_K) == pytest.approx(k, rel=1e-12)
        (x0, k0), (x1, k1) = table[0], table[1]
        assert f(0.5 * (x0 + x1) + to_K) == pytest.approx(0.5 * (k0 + k1), rel=1e-12)
        assert f(T_min) == pytest.approx(table[0][1], rel=1e-12)
        assert f(T_max) == pytest.approx(table[-1][1], rel=1e-12)
        with pytest.raises(RuntimeError):
            f(T_min - 1.0)
        with pytest.raises(RuntimeError):
            f(T_max + 1.0)

    def test_they_are_in_the_material_database(self) -> None:
        names = cb.list_materials()
        for n in ("haynes25", "inconel625", "nimonic75", "copper", "alumina", "x22crmov12_1"):
            assert n in names
            assert cb.get_material_conductivity(n, 700.0) == pytest.approx(
                getattr(cb, f"k_{n}")(700.0), rel=1e-12
            )

    @pytest.mark.parametrize(("func", "unit", "T_min", "T_max", "table"), TABLES)
    def test_the_database_clamps_to_the_table_ends(
        self, func: str, unit: str, T_min: float, T_max: float, table: list
    ) -> None:
        name = func.removeprefix("k_")
        assert cb.get_material_conductivity(name, 1.0) == pytest.approx(table[0][1], rel=1e-12)
        assert cb.get_material_conductivity(name, 5000.0) == pytest.approx(table[-1][1], rel=1e-12)


class TestOptInExtrapolation:
    """Outside the valid range a k_* function raises and the database holds
    the end value; extrapolate=True continues the end slope linearly."""

    def test_default_refuses_or_holds_the_end_value(self) -> None:
        with pytest.raises(RuntimeError):
            cb.k_copper(1300.0)
        assert cb.get_material_conductivity("copper", 1300.0) == pytest.approx(342.0, rel=1e-12)
        assert cb.get_material_conductivity("alumina", 50.0) == pytest.approx(133.0, rel=1e-12)

    def test_a_table_continues_its_end_segment(self) -> None:
        # Copper 1100 -> 1200 K: 350 -> 342, slope -0.08 per K.
        assert cb.k_copper(1300.0, extrapolate=True) == pytest.approx(334.0, rel=1e-9)
        # Alumina 100 -> 150 K: 133 -> 77, slope -1.12 per K, continued below 100 K.
        assert cb.k_alumina(50.0, extrapolate=True) == pytest.approx(189.0, rel=1e-9)
        # Hastelloy X 900 -> 1000 C: 26.7 -> 28.7.
        assert cb.k_hastelloy_x(1100.0 + 273.15, extrapolate=True) == pytest.approx(30.7, rel=1e-9)
        assert cb.get_material_conductivity("copper", 1300.0, extrapolate=True) == pytest.approx(
            334.0, rel=1e-9
        )

    def test_inside_the_range_the_flag_changes_nothing(self) -> None:
        for T in (300.0, 777.7, 1273.15):
            assert cb.k_hastelloy_x(T, extrapolate=True) == cb.k_hastelloy_x(T)

    def test_a_fitted_line_continues_its_tangent(self) -> None:
        # Inconel 718: 11.4 + 0.0146 T_C up to 1200 K.
        k_end = 11.4 + 0.0146 * (1200.0 - 273.15)
        assert cb.k_inconel718(1300.0, extrapolate=True) == pytest.approx(
            k_end + 100.0 * 0.0146, rel=1e-6
        )

    def test_haynes230_continues_its_end_segment_past_1000_c(self) -> None:
        # alpha cp rho, 900 -> 1000 C: 24.25 -> 26.68; the flag continues it.
        assert cb.k_haynes230(1100.0 + 273.15, extrapolate=True) == pytest.approx(
            26.68 + (26.68 - 24.25), rel=1e-9
        )
        assert cb.get_material_conductivity("haynes230", 1400.0) == pytest.approx(26.68, rel=1e-12)

    def test_it_refuses_where_the_continuation_is_unphysical(self) -> None:
        with pytest.raises(RuntimeError, match="<= 0"):
            cb.k_copper(6000.0, extrapolate=True)
        with pytest.raises(RuntimeError, match="> 0 K"):
            cb.k_copper(0.0, extrapolate=True)
        for bad in (float("nan"), float("inf")):
            with pytest.raises(RuntimeError, match="finite"):
                cb.k_copper(bad, extrapolate=True)
            with pytest.raises(RuntimeError, match="finite"):
                cb.get_material_conductivity("copper", bad, extrapolate=True)

    def test_a_wall_layer_opts_in(self) -> None:
        from combaero.network.components import WallLayer

        held = WallLayer(thickness=1e-3, conductivity=1.0, material="copper")
        held.update_conductivity(1300.0)
        assert held.conductivity == pytest.approx(342.0, rel=1e-12)
        continued = WallLayer(thickness=1e-3, conductivity=1.0, material="copper", extrapolate=True)
        continued.update_conductivity(1300.0)
        assert continued.conductivity == pytest.approx(334.0, rel=1e-9)
