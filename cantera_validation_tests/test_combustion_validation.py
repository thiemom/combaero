import numpy as np
import pytest

from combaero._core import species_index_from_name, num_species


class TestCompleteCombustion:
    """Validate complete combustion calculations against Cantera."""

    @pytest.fixture
    def gri30_gas(self, cantera):
        """Create Cantera gas object with GRI-Mech 3.0 (full mechanism for transport tests)."""
        return cantera.Solution("gri30.yaml")

    def set_cantera_composition(self, gas, cb_X, species_mapping):
        """Set Cantera composition from CombAero mole fractions."""
        ct_species = {}
        for i, name in enumerate(["N2", "O2", "AR", "CO2", "H2O", "CH4", "C2H6", "C3H8", 
                                   "IC4H10", "NC5H12", "NC6H14", "NC7H16", "CO", "H2"]):
            if cb_X[i] > 1e-10:
                ct_name = species_mapping.get(name, name)
                ct_species[ct_name] = cb_X[i]
        gas.X = ct_species
    
    def find_adiabatic_T_enthalpy_balance(self, gas, X_reactants, X_products, T_reactants, P, T_guess=2000.0, tol=0.1):
        """Find adiabatic flame temperature by enthalpy balance.
        
        Finds T_products such that H(T_reactants, X_reactants) = H(T_products, X_products)
        at constant pressure (adiabatic, constant pressure combustion).
        """
        # Calculate reactant enthalpy
        self.set_cantera_composition(gas, X_reactants, {})
        gas.TP = T_reactants, P
        H_reactants = gas.enthalpy_mass
        
        # Binary search for product temperature
        T_low, T_high = 300.0, 4000.0
        
        for _ in range(50):  # Max iterations
            T_mid = (T_low + T_high) / 2.0
            
            self.set_cantera_composition(gas, X_products, {})
            gas.TP = T_mid, P
            H_products = gas.enthalpy_mass
            
            if abs(H_products - H_reactants) < abs(H_reactants) * 1e-6:
                return T_mid
            
            if H_products < H_reactants:
                T_low = T_mid
            else:
                T_high = T_mid
        
        return T_mid

    def test_methane_air_stoichiometric(self, combaero, cantera, complete_combustion_gas, species_mapping, tolerance_config):
        """Test stoichiometric CH4 + air combustion using enthalpy balance."""
        cb = combaero
        
        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0
        
        phi = 1.0
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_air)
        
        T_in = 300.0
        P_in = 101325.0
        
        # CombAero complete combustion
        burned_cb = cb.complete_combustion(T=T_in, X=X_mix, P=P_in)
        
        # Cantera: Find Tad by enthalpy balance (not equilibrium)
        T_cantera = self.find_adiabatic_T_enthalpy_balance(
            complete_combustion_gas, X_mix, burned_cb.X, T_in, P_in
        )
        
        # Compare adiabatic flame temperatures [K]
        # CombAero: NASA-9 polynomials, complete combustion model
        # Cantera: NASA-7 polynomials, enthalpy balance (same products)
        T_diff = abs(burned_cb.T - T_cantera)
        print(f"\nCH4 stoich: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K, diff={T_diff:.1f} K")
        
        assert T_diff < tolerance_config["temperature"], \
            f"Temperature mismatch [K]: CombAero={burned_cb.T:.1f}, Cantera={T_cantera:.1f}, diff={T_diff:.1f}"
        
        # Verify product composition matches stoichiometry
        for species in ["CO2", "H2O", "N2"]:
            idx_cb = species_index_from_name(species)
            X_cb = burned_cb.X[idx_cb]
            assert X_cb > 0.0, f"{species} should be present in products"

    def test_methane_air_lean(self, combaero, cantera, complete_combustion_gas, species_mapping, tolerance_config):
        """Test lean CH4 + air combustion (phi=0.8) using enthalpy balance."""
        cb = combaero
        
        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0
        
        phi = 0.8
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_air)
        
        T_in = 300.0
        P_in = 101325.0
        
        # CombAero complete combustion
        burned_cb = cb.complete_combustion(T=T_in, X=X_mix, P=P_in)
        
        # Cantera: Find Tad by enthalpy balance
        T_cantera = self.find_adiabatic_T_enthalpy_balance(
            complete_combustion_gas, X_mix, burned_cb.X, T_in, P_in
        )
        
        T_diff = abs(burned_cb.T - T_cantera)
        print(f"\nCH4 lean (phi=0.8): CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K, diff={T_diff:.1f} K")
        
        assert T_diff < tolerance_config["temperature"], \
            f"Temperature mismatch: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K"
        
        # Verify excess O2 present
        idx_O2 = species_index_from_name("O2")
        X_O2_cb = burned_cb.X[idx_O2]
        assert X_O2_cb > 0.0, "Lean combustion should have excess O2"

    def test_propane_air_stoichiometric(self, combaero, cantera, complete_combustion_gas, species_mapping, tolerance_config):
        """Test stoichiometric C3H8 + air combustion using enthalpy balance."""
        cb = combaero
        
        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("C3H8")] = 1.0
        
        phi = 1.0
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_air)
        
        T_in = 300.0
        P_in = 101325.0
        
        # CombAero complete combustion
        burned_cb = cb.complete_combustion(T=T_in, X=X_mix, P=P_in)
        
        # Cantera: Find Tad by enthalpy balance
        T_cantera = self.find_adiabatic_T_enthalpy_balance(
            complete_combustion_gas, X_mix, burned_cb.X, T_in, P_in
        )
        
        T_diff = abs(burned_cb.T - T_cantera)
        print(f"\nC3H8 stoich: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K, diff={T_diff:.1f} K")
        
        assert T_diff < tolerance_config["temperature"], \
            f"Temperature mismatch: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K"

    def test_hydrogen_air_stoichiometric(self, combaero, cantera, complete_combustion_gas, species_mapping, tolerance_config):
        """Test stoichiometric H2 + air combustion using enthalpy balance."""
        cb = combaero
        
        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("H2")] = 1.0
        
        phi = 1.0
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_air)
        
        T_in = 300.0
        P_in = 101325.0
        
        # CombAero complete combustion
        burned_cb = cb.complete_combustion(T=T_in, X=X_mix, P=P_in)
        
        # Cantera: Find Tad by enthalpy balance
        T_cantera = self.find_adiabatic_T_enthalpy_balance(
            complete_combustion_gas, X_mix, burned_cb.X, T_in, P_in
        )
        
        T_diff = abs(burned_cb.T - T_cantera)
        print(f"\nH2 stoich: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K, diff={T_diff:.1f} K")
        
        assert T_diff < tolerance_config["temperature"], \
            f"Temperature mismatch: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K"
        
        # Verify H2O production
        idx_H2O = species_index_from_name("H2O")
        X_H2O_cb = burned_cb.X[idx_H2O]
        assert X_H2O_cb > 0.0, "H2O should be present in products"

    def test_temperature_variation(self, combaero, cantera, complete_combustion_gas, species_mapping, tolerance_config):
        """Test combustion at different inlet temperatures using enthalpy balance."""
        cb = combaero
        
        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0
        
        phi = 1.0
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_air)
        P_in = 101325.0
        
        max_diff = 0.0
        for T_in in [300.0, 400.0, 500.0, 600.0]:
            burned_cb = cb.complete_combustion(T=T_in, X=X_mix, P=P_in)
            
            T_cantera = self.find_adiabatic_T_enthalpy_balance(
                complete_combustion_gas, X_mix, burned_cb.X, T_in, P_in
            )
            
            T_diff = abs(burned_cb.T - T_cantera)
            max_diff = max(max_diff, T_diff)
            print(f"\nT_in={T_in} K: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K, diff={T_diff:.1f} K")
            
            assert T_diff < tolerance_config["temperature"], \
                f"T_in={T_in} K: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K"
        
        print(f"\nMax temperature deviation: {max_diff:.1f} K")

    def test_pressure_variation(self, combaero, cantera, complete_combustion_gas, species_mapping, tolerance_config):
        """Test combustion at different pressures using enthalpy balance."""
        cb = combaero
        
        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0
        
        phi = 1.0
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_air)
        T_in = 300.0
        
        max_diff = 0.0
        for P_in in [101325.0, 200000.0, 500000.0, 1000000.0]:
            burned_cb = cb.complete_combustion(T=T_in, X=X_mix, P=P_in)
            
            T_cantera = self.find_adiabatic_T_enthalpy_balance(
                complete_combustion_gas, X_mix, burned_cb.X, T_in, P_in
            )
            
            T_diff = abs(burned_cb.T - T_cantera)
            max_diff = max(max_diff, T_diff)
            print(f"\nP_in={P_in/1e5:.1f} bar: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K, diff={T_diff:.1f} K")
            
            assert T_diff < tolerance_config["temperature"], \
                f"P_in={P_in/1e5:.1f} bar: CombAero={burned_cb.T:.1f} K, Cantera={T_cantera:.1f} K"
        
        print(f"\nMax temperature deviation: {max_diff:.1f} K")


class TestOxygenRequirement:
    """Validate oxygen requirement calculations."""

    def test_methane_oxygen_requirement(self, combaero):
        """Test CH4 oxygen requirement: CH4 + 2 O2 -> CO2 + 2 H2O."""
        cb = combaero
        idx_CH4 = species_index_from_name("CH4")
        
        # Compare stoichiometric O2 requirement [mol O2 / mol fuel]
        O2_per_mol = cb.oxygen_required_per_mol_fuel(idx_CH4)
        assert abs(O2_per_mol - 2.0) < 1e-10, \
            f"O2 requirement [mol O2/mol CH4]: Expected 2.0, got {O2_per_mol}"

    def test_propane_oxygen_requirement(self, combaero):
        """Test C3H8 oxygen requirement: C3H8 + 5 O2 -> 3 CO2 + 4 H2O."""
        cb = combaero
        idx_C3H8 = species_index_from_name("C3H8")
        
        O2_per_mol = cb.oxygen_required_per_mol_fuel(idx_C3H8)
        assert abs(O2_per_mol - 5.0) < 1e-10, f"Expected 5.0, got {O2_per_mol}"

    def test_hydrogen_oxygen_requirement(self, combaero):
        """Test H2 oxygen requirement: 2 H2 + O2 -> 2 H2O."""
        cb = combaero
        idx_H2 = species_index_from_name("H2")
        
        O2_per_mol = cb.oxygen_required_per_mol_fuel(idx_H2)
        assert abs(O2_per_mol - 0.5) < 1e-10, f"Expected 0.5, got {O2_per_mol}"


class TestEquivalenceRatio:
    """Validate equivalence ratio calculations."""

    def test_stoichiometric_mixture(self, combaero):
        """Test that phi=1.0 produces stoichiometric mixture."""
        cb = combaero
        
        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0
        
        phi_target = 1.0
        X_mix = cb.set_equivalence_ratio_mole(phi_target, X_fuel, X_air)
        
        # Compare equivalence ratio [dimensionless]
        # Round-trip test: set φ=1.0, calculate φ from mixture
        phi_calc = cb.equivalence_ratio_mole(X_mix, X_fuel, X_air)
        assert abs(phi_calc - phi_target) < 1e-6, \
            f"Equivalence ratio [dimensionless]: Expected {phi_target}, got {phi_calc}"

    def test_lean_mixture(self, combaero):
        """Test lean mixture (phi < 1)."""
        cb = combaero
        
        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0
        
        phi_target = 0.8
        X_mix = cb.set_equivalence_ratio_mole(phi_target, X_fuel, X_air)
        
        phi_calc = cb.equivalence_ratio_mole(X_mix, X_fuel, X_air)
        assert abs(phi_calc - phi_target) < 1e-6, \
            f"Expected phi={phi_target}, got {phi_calc}"

    def test_rich_mixture(self, combaero):
        """Test rich mixture (phi > 1)."""
        cb = combaero
        
        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0
        
        phi_target = 1.2
        X_mix = cb.set_equivalence_ratio_mole(phi_target, X_fuel, X_air)
        
        phi_calc = cb.equivalence_ratio_mole(X_mix, X_fuel, X_air)
        assert abs(phi_calc - phi_target) < 1e-6, \
            f"Expected phi={phi_target}, got {phi_calc}"
