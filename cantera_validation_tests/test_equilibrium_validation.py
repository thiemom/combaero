import numpy as np
import pytest

from combaero._core import num_species, species_index_from_name


class TestWgsEquilibrium:
    """Validate WGS (Water-Gas Shift) equilibrium against Cantera.

    Reaction: CO + H2O ⇌ CO2 + H2

    CombAero uses partial equilibrium (WGS only), Cantera uses restricted
    species to match the same equilibrium model.
    """

    @pytest.fixture
    def wgs_gas(self, cantera):
        """Create Cantera gas with WGS species only."""
        species = {S.name: S for S in cantera.Species.list_from_file("gri30.yaml")}
        wgs_species = [species[S] for S in ("CO", "H2O", "CO2", "H2", "N2", "AR")]
        return cantera.Solution(thermo="ideal-gas", species=wgs_species)

    def set_cantera_composition(self, gas, cb_X, species_mapping):
        """Set Cantera composition from CombAero mole fractions."""
        ct_species = {}
        for i, name in enumerate(
            [
                "N2",
                "O2",
                "AR",
                "CO2",
                "H2O",
                "CH4",
                "C2H6",
                "C3H8",
                "IC4H10",
                "NC5H12",
                "NC6H14",
                "NC7H16",
                "CO",
                "H2",
            ]
        ):
            if cb_X[i] > 1e-10:
                ct_name = species_mapping.get(name, name)
                if ct_name in gas.species_names:
                    ct_species[ct_name] = cb_X[i]
        gas.X = ct_species

    def test_wgs_isothermal_high_temp(
        self, combaero, cantera, wgs_gas, species_mapping, tolerance_config
    ):
        """Test WGS equilibrium at high temperature (1200 K) - favors CO + H2."""
        cb = combaero

        # Initial composition: CO + H2O rich
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        T = 1200.0
        P = 101325.0

        # CombAero WGS equilibrium
        result_cb = cb.wgs_equilibrium(T, X, P)

        # Cantera equilibrium with WGS species only
        self.set_cantera_composition(wgs_gas, X, species_mapping)
        wgs_gas.TP = T, P
        wgs_gas.equilibrate("TP")  # Isothermal equilibrium

        # Compare temperature [K] - should be unchanged (isothermal)
        assert (
            abs(result_cb.T - T) < 0.1
        ), f"Temperature should be unchanged: {result_cb.T:.1f} K vs {T:.1f} K"

        # Compare equilibrium composition [mol/mol]
        # At high T, WGS favors reactants (CO + H2O), so expect shift toward products
        print(f"\nWGS isothermal @ {T} K:")
        for species in ["CO", "H2O", "CO2", "H2"]:
            idx_cb = species_index_from_name(species)
            X_cb = result_cb.X[idx_cb]
            X_ct = wgs_gas[species].X[0] if species in wgs_gas.species_names else 0.0
            diff = abs(X_cb - X_ct)
            print(f"  {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

            # Major species tolerance
            if X_cb > 0.01 or X_ct > 0.01:
                assert (
                    diff < 0.005
                ), f"{species} [mol/mol]: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}"

    def test_wgs_isothermal_low_temp(
        self, combaero, cantera, wgs_gas, species_mapping, tolerance_config
    ):
        """Test WGS equilibrium at low temperature (800 K) - favors CO2 + H2."""
        cb = combaero

        # Initial composition: CO + H2O rich
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        T = 800.0
        P = 101325.0

        # CombAero WGS equilibrium
        result_cb = cb.wgs_equilibrium(T, X, P)

        # Cantera equilibrium
        self.set_cantera_composition(wgs_gas, X, species_mapping)
        wgs_gas.TP = T, P
        wgs_gas.equilibrate("TP")

        # Compare composition
        print(f"\nWGS isothermal @ {T} K:")
        for species in ["CO", "H2O", "CO2", "H2"]:
            idx_cb = species_index_from_name(species)
            X_cb = result_cb.X[idx_cb]
            X_ct = wgs_gas[species].X[0] if species in wgs_gas.species_names else 0.0
            diff = abs(X_cb - X_ct)
            print(f"  {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

            if X_cb > 0.01 or X_ct > 0.01:
                assert diff < 0.005, f"{species} [mol/mol]: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}"

    def test_wgs_isothermal_temperature_sweep(self, combaero, cantera, wgs_gas, species_mapping):
        """Test WGS equilibrium across temperature range (800-1500 K)."""
        cb = combaero

        # Initial composition
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        P = 101325.0

        print("\nWGS isothermal temperature sweep:")
        max_deviations = {"CO": 0.0, "H2O": 0.0, "CO2": 0.0, "H2": 0.0}

        for T in [800.0, 1000.0, 1200.0, 1500.0]:
            result_cb = cb.wgs_equilibrium(T, X, P)

            self.set_cantera_composition(wgs_gas, X, species_mapping)
            wgs_gas.TP = T, P
            wgs_gas.equilibrate("TP")

            print(f"\n  T = {T} K:")
            for species in ["CO", "H2O", "CO2", "H2"]:
                idx_cb = species_index_from_name(species)
                X_cb = result_cb.X[idx_cb]
                X_ct = wgs_gas[species].X[0]
                diff = abs(X_cb - X_ct)
                max_deviations[species] = max(max_deviations[species], diff)
                print(f"    {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

        print("\n  Maximum deviations:")
        for species, max_diff in max_deviations.items():
            print(f"    {species}: {max_diff:.6f}")

    def test_wgs_adiabatic(self, combaero, cantera, wgs_gas, species_mapping, tolerance_config):
        """Test adiabatic WGS equilibrium - temperature changes with reaction."""
        cb = combaero

        # Initial composition at high temperature
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        T_in = 1500.0
        P = 101325.0

        # CombAero adiabatic WGS
        result_cb = cb.wgs_equilibrium_adiabatic(T_in, X, P)

        # Cantera adiabatic equilibrium
        self.set_cantera_composition(wgs_gas, X, species_mapping)
        wgs_gas.TP = T_in, P
        wgs_gas.equilibrate("HP")  # Adiabatic (constant enthalpy, pressure)

        # Compare final temperature [K]
        # WGS is exothermic in forward direction (CO + H2O → CO2 + H2)
        # At high T, equilibrium favors reactants, so forward reaction releases heat
        T_diff = abs(result_cb.T - wgs_gas.T)
        print("\nWGS adiabatic:")
        print(f"  T_in: {T_in:.1f} K")
        print(
            f"  T_final: CombAero={result_cb.T:.1f} K, "
            f"Cantera={wgs_gas.T:.1f} K, diff={T_diff:.1f} K"
        )

        assert T_diff < 10.0, (
            f"Temperature [K]: CombAero={result_cb.T:.1f}, "
            f"Cantera={wgs_gas.T:.1f}, diff={T_diff:.1f}"
        )

        # Compare equilibrium composition
        print("  Composition:")
        for species in ["CO", "H2O", "CO2", "H2"]:
            idx_cb = species_index_from_name(species)
            X_cb = result_cb.X[idx_cb]
            X_ct = wgs_gas[species].X[0]
            diff = abs(X_cb - X_ct)
            print(f"    {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

            if X_cb > 0.01 or X_ct > 0.01:
                assert diff < 0.005, f"{species} [mol/mol]: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}"

    def test_wgs_adiabatic_temperature_sweep(self, combaero, cantera, wgs_gas, species_mapping):
        """Test adiabatic WGS across inlet temperatures."""
        cb = combaero

        # Initial composition
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        P = 101325.0

        print("\nWGS adiabatic temperature sweep:")
        max_T_diff = 0.0

        for T_in in [1000.0, 1200.0, 1500.0, 1800.0]:
            result_cb = cb.wgs_equilibrium_adiabatic(T_in, X, P)

            self.set_cantera_composition(wgs_gas, X, species_mapping)
            wgs_gas.TP = T_in, P
            wgs_gas.equilibrate("HP")

            T_diff = abs(result_cb.T - wgs_gas.T)
            max_T_diff = max(max_T_diff, T_diff)

            print(f"\n  T_in = {T_in} K:")
            print(
                f"    T_final: CombAero={result_cb.T:.1f} K, "
                f"Cantera={wgs_gas.T:.1f} K, diff={T_diff:.1f} K"
            )

            for species in ["CO", "H2O", "CO2", "H2"]:
                idx_cb = species_index_from_name(species)
                X_cb = result_cb.X[idx_cb]
                X_ct = wgs_gas[species].X[0]
                diff = abs(X_cb - X_ct)
                print(f"    {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

        print(f"\n  Maximum temperature deviation: {max_T_diff:.1f} K")

    def test_wgs_equilibrium_constant(self, combaero, cantera, wgs_gas, species_mapping):
        """Test that equilibrium constant Kp matches between CombAero and Cantera."""
        cb = combaero

        # Initial composition
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        P = 101325.0

        print("\nWGS equilibrium constant Kp:")

        for T in [800.0, 1000.0, 1200.0, 1500.0]:
            # CombAero equilibrium
            result_cb = cb.wgs_equilibrium(T, X, P)

            # Calculate Kp from composition: Kp = (X_CO2 * X_H2) / (X_CO * X_H2O)
            X_CO2_cb = result_cb.X[species_index_from_name("CO2")]
            X_H2_cb = result_cb.X[species_index_from_name("H2")]
            X_CO_cb = result_cb.X[species_index_from_name("CO")]
            X_H2O_cb = result_cb.X[species_index_from_name("H2O")]
            Kp_cb = (
                (X_CO2_cb * X_H2_cb) / (X_CO_cb * X_H2O_cb) if (X_CO_cb * X_H2O_cb) > 1e-10 else 0.0
            )

            # Cantera equilibrium
            self.set_cantera_composition(wgs_gas, X, species_mapping)
            wgs_gas.TP = T, P
            wgs_gas.equilibrate("TP")

            X_CO2_ct = wgs_gas["CO2"].X[0]
            X_H2_ct = wgs_gas["H2"].X[0]
            X_CO_ct = wgs_gas["CO"].X[0]
            X_H2O_ct = wgs_gas["H2O"].X[0]
            Kp_ct = (
                (X_CO2_ct * X_H2_ct) / (X_CO_ct * X_H2O_ct) if (X_CO_ct * X_H2O_ct) > 1e-10 else 0.0
            )

            rel_diff = abs(Kp_cb - Kp_ct) / Kp_ct if Kp_ct > 1e-10 else 0.0

            print(
                f"  T = {T} K: CombAero Kp={Kp_cb:.4f}, "
                f"Cantera Kp={Kp_ct:.4f}, diff={rel_diff*100:.2f}%"
            )

            # Kp should match within 5% (due to Gibbs free energy polynomial differences)
            assert rel_diff < 0.05, (
                f"Kp [dimensionless] @ {T} K: CombAero={Kp_cb:.4f}, "
                f"Cantera={Kp_ct:.4f}, diff={rel_diff*100:.1f}%"
            )
