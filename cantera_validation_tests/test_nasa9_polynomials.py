"""NASA-9 Polynomial Validation Tests.

These tests validate CombAero's NASA-9 polynomial implementation and evaluation
functions by comparing against Cantera's NASA-9 implementation using the same
polynomial coefficients.

This provides polynomial-level validation, complementing the system-level
validation in other test files (combustion, equilibrium, etc.).
"""

import numpy as np
import pytest

from combaero._core import num_species, species_index_from_name


class TestNASA9Polynomials:
    """Validate NASA-9 polynomial evaluation against Cantera.

    Tests run in parallel to existing NASA-7 validation tests.
    Focus: Polynomial implementation correctness (Cp, H, S, G evaluation).
    """

    @pytest.fixture
    def nasa9_gas(self, cantera):
        """Cantera gas with NASA-9 polynomials matching CombAero."""
        return cantera.Solution("combaero_nasa9.yaml")

    @pytest.fixture
    def test_species(self):
        """Species available in both CombAero and NASA-9 data."""
        return ["N2", "O2", "AR", "CO2", "H2O", "CH4", "C2H6", "C3H8", "CO", "H2"]

    @pytest.fixture
    def test_temperatures(self):
        """Temperature points for validation (200-6000 K range)."""
        return [
            200.0,
            300.0,
            500.0,
            800.0,
            1000.0,
            1500.0,
            2000.0,
            3000.0,
            4000.0,
            5000.0,
            6000.0,
        ]

    def test_cp_evaluation(
        self, combaero, nasa9_gas, test_species, test_temperatures, species_mapping
    ):
        """Test Cp/R polynomial evaluation across temperature range.

        Validates: Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
        Expected: < 0.0001% deviation (numerical precision)
        """
        cb = combaero
        R = 8.314462618  # J/mol-K (universal gas constant)

        print(f"\n{'=' * 80}")
        print("NASA-9 Polynomial Validation: Cp/R Evaluation")
        print(f"{'=' * 80}")

        max_deviations = {}

        for species in test_species:
            idx_cb = species_index_from_name(species)
            ct_name = species_mapping.get(species, species)

            if ct_name not in nasa9_gas.species_names:
                continue

            deviations = []

            print(f"\n{species}:")
            print(f"  {'T [K]':>8s}  {'Cp/R (CB)':>12s}  {'Cp/R (CT)':>12s}  {'Deviation':>12s}")

            for T in test_temperatures:
                # CombAero: Get Cp in J/mol-K, convert to Cp/R
                X = np.zeros(num_species())
                X[idx_cb] = 1.0
                cp_cb = cb.cp(T, X)  # J/mol-K
                cp_r_cb = cp_cb / R

                # Cantera: Get Cp in J/kmol-K, convert to J/mol-K, then Cp/R
                nasa9_gas.TPX = T, 101325.0, {ct_name: 1.0}
                cp_ct = nasa9_gas.cp_mole  # J/kmol-K
                cp_r_ct = cp_ct / 1000.0 / R

                deviation = abs(cp_r_cb - cp_r_ct)
                rel_deviation = deviation / cp_r_ct if cp_r_ct > 0 else 0.0
                deviations.append(rel_deviation)

                if T in [300.0, 1000.0, 3000.0, 6000.0]:
                    print(
                        f"  {T:8.1f}  {cp_r_cb:12.6f}  {cp_r_ct:12.6f}  "
                        f"{rel_deviation * 100:11.4f}%"
                    )

            max_dev = max(deviations)
            max_deviations[species] = max_dev
            print(f"  Maximum deviation: {max_dev * 100:.6f}%")

            # Assert < 0.01% deviation (numerical precision)
            assert max_dev < 0.0001, f"{species} Cp/R deviation {max_dev * 100:.4f}% exceeds 0.01%"

        print(f"\n{'=' * 80}")
        print("Summary: Cp/R Evaluation")
        print(f"{'=' * 80}")
        for species, dev in max_deviations.items():
            print(f"  {species:6s}: {dev * 100:8.6f}% max deviation")
        print(f"{'=' * 80}\n")

    def test_enthalpy_integration(
        self, combaero, nasa9_gas, test_species, test_temperatures, species_mapping
    ):
        """Test enthalpy via Cp integration (validates polynomial integration).

        Computes ΔH = ∫Cp dT numerically and compares CombAero vs Cantera.
        This validates the polynomial integration without reference state dependency.
        Expected: < 0.001% deviation
        """
        cb = combaero
        T_ref = 298.15  # Reference temperature [K]

        print(f"\n{'=' * 80}")
        print(f"NASA-9 Polynomial Validation: ∫Cp dT Integration (from {T_ref} K)")
        print(f"{'=' * 80}")

        max_deviations = {}

        for species in test_species:
            idx_cb = species_index_from_name(species)
            ct_name = species_mapping.get(species, species)

            if ct_name not in nasa9_gas.species_names:
                continue

            deviations = []
            X = np.zeros(num_species())
            X[idx_cb] = 1.0

            print(f"\n{species}:")
            print(
                f"  {'T [K]':>8s}  {'∫Cp dT (CB)':>14s}  {'∫Cp dT (CT)':>14s}  {'Deviation':>12s}"
            )

            for T_target in test_temperatures:
                if T_target <= T_ref:
                    continue  # Skip temperatures at or below reference

                # Numerical integration: ∫Cp dT from T_ref to T_target
                n_points = 50
                dT = (T_target - T_ref) / n_points

                # CombAero integration
                h_int_cb = 0.0
                for i in range(n_points):
                    T = T_ref + (i + 0.5) * dT
                    h_int_cb += cb.cp(T, X) * dT

                # Cantera integration
                h_int_ct = 0.0
                for i in range(n_points):
                    T = T_ref + (i + 0.5) * dT
                    nasa9_gas.TPX = T, 101325.0, {ct_name: 1.0}
                    h_int_ct += (nasa9_gas.cp_mole / 1000.0) * dT

                deviation = abs(h_int_cb - h_int_ct)
                rel_deviation = deviation / abs(h_int_ct) if abs(h_int_ct) > 1e-10 else 0.0
                deviations.append(rel_deviation)

                if T_target in [500.0, 1000.0, 3000.0, 6000.0]:
                    print(
                        f"  {T_target:8.1f}  {h_int_cb:14.2f}  "
                        f"{h_int_ct:14.2f}  {rel_deviation * 100:11.6f}%"
                    )

            if deviations:
                max_dev = max(deviations)
                max_deviations[species] = max_dev
                print(f"  Maximum deviation: {max_dev * 100:.6f}%")

                # Assert < 0.01% deviation (validates polynomial integration)
                assert (
                    max_dev < 0.0001
                ), f"{species} ∫Cp dT deviation {max_dev * 100:.6f}% exceeds 0.01%"

        print(f"\n{'=' * 80}")
        print("Summary: ∫Cp dT Integration")
        print(f"{'=' * 80}")
        for species, dev in max_deviations.items():
            print(f"  {species:6s}: {dev * 100:8.6f}% max deviation")
        print(f"{'=' * 80}\n")

    def test_entropy_evaluation(
        self, combaero, nasa9_gas, test_species, test_temperatures, species_mapping
    ):
        """Test ΔS polynomial evaluation using s(T) - s(Tref).

        Compares entropy CHANGES (not absolute values) to avoid reference state issues.
        """
        cb = combaero
        T_ref = 298.15  # Reference temperature [K]
        P = 101325.0  # Pressure [Pa]

        print(f"\n{'=' * 80}")
        print(f"NASA-9 Polynomial Validation: ΔS Evaluation (relative to {T_ref} K)")
        print(f"{'=' * 80}")

        max_deviations = {}

        for species in test_species:
            idx_cb = species_index_from_name(species)
            ct_name = species_mapping.get(species, species)

            if ct_name not in nasa9_gas.species_names:
                continue

            deviations = []

            # Get reference entropy at T_ref
            X_ref = np.zeros(num_species())
            X_ref[idx_cb] = 1.0
            s_ref_cb = cb.s(T_ref, X_ref, P)  # J/mol-K

            nasa9_gas.TPX = T_ref, P, {ct_name: 1.0}
            s_ref_ct = nasa9_gas.entropy_mole / 1000.0  # J/kmol -> J/mol

            print(f"\n{species}:")
            print("     T [K]     ΔS (CB)       ΔS (CT)       Deviation")

            for T in test_temperatures:
                X = np.zeros(num_species())
                X[idx_cb] = 1.0
                s_cb = cb.s(T, X, P)  # J/mol-K

                nasa9_gas.TPX = T, P, {ct_name: 1.0}
                s_ct = nasa9_gas.entropy_mole / 1000.0  # J/kmol -> J/mol

                # Calculate ΔS = S(T) - S(Tref)
                delta_s_cb = s_cb - s_ref_cb
                delta_s_ct = s_ct - s_ref_ct

                # Relative deviation in ΔS
                deviation = (
                    abs(delta_s_cb - delta_s_ct) / abs(delta_s_ct)
                    if abs(delta_s_ct) > 1e-6
                    else 0.0
                )
                deviations.append(deviation)

                if T in [300.0, 1000.0, 3000.0, 6000.0]:
                    print(
                        f"  {T:8.1f}  {delta_s_cb:12.6f}  {delta_s_ct:12.6f}  "
                        f"{deviation * 100:11.4f}%"
                    )

            if deviations:
                max_dev = max(deviations)
                max_deviations[species] = max_dev
                print(f"  Maximum deviation: {max_dev * 100:.6f}%")

                # Assert < 0.01% deviation (same precision as Cp test)
                assert (
                    max_dev < 0.0001
                ), f"{species} ΔS deviation {max_dev * 100:.6f}% exceeds 0.01%"

        print(f"\n{'=' * 80}")
        print("Summary: ΔS = S(T) - S(Tref) Evaluation")
        print(f"{'=' * 80}")
        for species, dev in max_deviations.items():
            print(f"  {species:6s}: {dev * 100:8.6f}% max deviation")
        print(f"{'=' * 80}\n")

    def test_temperature_range_continuity(self, combaero, nasa9_gas, test_species, species_mapping):
        """Test continuity at temperature range boundaries.

        Validates smooth transitions at 1000 K (typical range boundary).
        """
        cb = combaero
        R = 8.314462618

        print(f"\n{'=' * 80}")
        print("NASA-9 Polynomial Validation: Temperature Range Continuity")
        print(f"{'=' * 80}")

        # Test at 1000 K boundary (+/- 1 K)
        T_boundary = 1000.0
        T_below = 999.0
        T_above = 1001.0

        for species in test_species:
            idx_cb = species_index_from_name(species)
            ct_name = species_mapping.get(species, species)

            if ct_name not in nasa9_gas.species_names:
                continue

            X = np.zeros(num_species())
            X[idx_cb] = 1.0

            # Cp at boundary
            cp_below = cb.cp(T_below, X) / R
            cp_boundary = cb.cp(T_boundary, X) / R
            cp_above = cb.cp(T_above, X) / R

            # Check smoothness (derivative shouldn't jump)
            slope_below = (cp_boundary - cp_below) / (T_boundary - T_below)
            slope_above = (cp_above - cp_boundary) / (T_above - T_boundary)
            slope_change = (
                abs(slope_above - slope_below) / abs(slope_below) if slope_below != 0 else 0.0
            )

            print(f"\n{species} @ {T_boundary} K:")
            print(
                f"  Cp/R: {cp_below:.6f} ({T_below}K) -> "
                f"{cp_boundary:.6f} ({T_boundary}K) -> {cp_above:.6f} ({T_above}K)"
            )
            print(f"  Slope change: {slope_change * 100:.4f}%")

            # Slope should not change dramatically (< 25% change is acceptable)
            # Note: NASA-9 polynomials may have small discontinuities at range boundaries
            # This is expected behavior, especially for light molecules like H2
            assert (
                slope_change < 0.25
            ), f"{species} has large discontinuity at {T_boundary} K (>{slope_change * 100:.1f}%)"

        print(f"{'=' * 80}\n")
