import numpy as np
import pytest

from combaero._core import num_species, species_index_from_name


class TestStreamMixing:
    """Validate stream mixing calculations against Cantera."""

    @pytest.fixture
    def gri30_gas(self, cantera):
        """Create Cantera gas object with GRI-Mech 3.0."""
        return cantera.Solution("gri30.yaml")

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
                ct_species[ct_name] = cb_X[i]
        gas.X = ct_species

    def test_two_stream_mixing_equal_mass(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test mixing two streams with equal mass flow rates."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0

        stream1 = cb.Stream()
        stream1.T = 300.0
        stream1.P = 101325.0
        stream1.X = X_air
        stream1.mdot = 1.0

        stream2 = cb.Stream()
        stream2.T = 400.0
        stream2.P = 101325.0
        stream2.X = X_fuel
        stream2.mdot = 1.0

        mixed_cb = cb.mix([stream1, stream2])

        self.set_cantera_composition(gri30_gas, X_air, species_mapping)
        gri30_gas.TP = stream1.T, stream1.P
        h1 = gri30_gas.enthalpy_mass

        self.set_cantera_composition(gri30_gas, X_fuel, species_mapping)
        gri30_gas.TP = stream2.T, stream2.P
        h2 = gri30_gas.enthalpy_mass

        h_mixed_target = 0.5 * (h1 + h2)

        X_mixed_target = 0.5 * X_air + 0.5 * X_fuel
        self.set_cantera_composition(gri30_gas, X_mixed_target, species_mapping)
        gri30_gas.HP = h_mixed_target, mixed_cb.P

        # Compare mixed stream temperature [K]
        # Both use enthalpy balance: H_mixed = (m1*H1 + m2*H2) / (m1 + m2)
        # Then find T where H(T, X_mixed) = H_mixed
        T_diff = abs(mixed_cb.T - gri30_gas.T)
        assert T_diff < tolerance_config["temperature"], (
            f"Temperature [K]: CombAero={mixed_cb.T:.1f}, "
            f"Cantera={gri30_gas.T:.1f}, diff={T_diff:.1f}"
        )

        # Compare total mass flow rate [kg/s]
        # Should be exact: mdot_total = mdot1 + mdot2
        assert abs(mixed_cb.mdot - 2.0) < 1e-10, (
            f"Mass flow rate [kg/s]: Expected 2.0, got {mixed_cb.mdot}"
        )

    def test_two_stream_mixing_unequal_mass(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test mixing two streams with different mass flow rates."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0

        mdot1 = 10.0
        mdot2 = 0.5

        stream1 = cb.Stream()
        stream1.T = 400.0
        stream1.P = 101325.0
        stream1.X = X_air
        stream1.mdot = mdot1

        stream2 = cb.Stream()
        stream2.T = 300.0
        stream2.P = 101325.0
        stream2.X = X_fuel
        stream2.mdot = mdot2

        mixed_cb = cb.mix([stream1, stream2])

        self.set_cantera_composition(gri30_gas, X_air, species_mapping)
        gri30_gas.TP = stream1.T, stream1.P
        h1 = gri30_gas.enthalpy_mass

        self.set_cantera_composition(gri30_gas, X_fuel, species_mapping)
        gri30_gas.TP = stream2.T, stream2.P
        h2 = gri30_gas.enthalpy_mass

        h_mixed_target = (mdot1 * h1 + mdot2 * h2) / (mdot1 + mdot2)

        Y1 = cb.mole_to_mass(X_air)
        Y2 = cb.mole_to_mass(X_fuel)
        Y_mixed = (mdot1 * Y1 + mdot2 * Y2) / (mdot1 + mdot2)
        X_mixed_target = cb.mass_to_mole(Y_mixed)

        self.set_cantera_composition(gri30_gas, X_mixed_target, species_mapping)
        gri30_gas.HP = h_mixed_target, mixed_cb.P

        T_diff = abs(mixed_cb.T - gri30_gas.T)
        assert T_diff < tolerance_config["temperature"], (
            f"Temperature mismatch: CombAero={mixed_cb.T:.1f} K, Cantera={gri30_gas.T:.1f} K"
        )

        assert abs(mixed_cb.mdot - (mdot1 + mdot2)) < 1e-10, (
            f"Mass flow should be {mdot1 + mdot2} kg/s, got {mixed_cb.mdot}"
        )

    def test_three_stream_mixing(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test mixing three streams."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0
        X_steam = np.zeros(len(X_air))
        X_steam[species_index_from_name("H2O")] = 1.0

        stream1 = cb.Stream()
        stream1.T = 300.0
        stream1.P = 101325.0
        stream1.X = X_air
        stream1.mdot = 10.0

        stream2 = cb.Stream()
        stream2.T = 300.0
        stream2.P = 101325.0
        stream2.X = X_fuel
        stream2.mdot = 0.5

        stream3 = cb.Stream()
        stream3.T = 500.0
        stream3.P = 101325.0
        stream3.X = X_steam
        stream3.mdot = 1.0

        mixed_cb = cb.mix([stream1, stream2, stream3])

        self.set_cantera_composition(gri30_gas, X_air, species_mapping)
        gri30_gas.TP = stream1.T, stream1.P
        h1 = gri30_gas.enthalpy_mass

        self.set_cantera_composition(gri30_gas, X_fuel, species_mapping)
        gri30_gas.TP = stream2.T, stream2.P
        h2 = gri30_gas.enthalpy_mass

        self.set_cantera_composition(gri30_gas, X_steam, species_mapping)
        gri30_gas.TP = stream3.T, stream3.P
        h3 = gri30_gas.enthalpy_mass

        mdot_total = stream1.mdot + stream2.mdot + stream3.mdot
        h_mixed_target = (stream1.mdot * h1 + stream2.mdot * h2 + stream3.mdot * h3) / mdot_total

        Y1 = cb.mole_to_mass(X_air)
        Y2 = cb.mole_to_mass(X_fuel)
        Y3 = cb.mole_to_mass(X_steam)
        Y_mixed = (stream1.mdot * Y1 + stream2.mdot * Y2 + stream3.mdot * Y3) / mdot_total
        X_mixed_target = cb.mass_to_mole(Y_mixed)

        self.set_cantera_composition(gri30_gas, X_mixed_target, species_mapping)
        gri30_gas.HP = h_mixed_target, mixed_cb.P

        T_diff = abs(mixed_cb.T - gri30_gas.T)
        assert T_diff < tolerance_config["temperature"], (
            f"Temperature mismatch: CombAero={mixed_cb.T:.1f} K, Cantera={gri30_gas.T:.1f} K"
        )

    def test_enthalpy_conservation(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test that mixing conserves enthalpy."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0

        stream1 = cb.Stream()
        stream1.T = 300.0
        stream1.P = 101325.0
        stream1.X = X_air
        stream1.mdot = 10.0

        stream2 = cb.Stream()
        stream2.T = 500.0
        stream2.P = 101325.0
        stream2.X = X_fuel
        stream2.mdot = 0.5

        h1_cb = cb.h(stream1.T, stream1.X)
        h2_cb = cb.h(stream2.T, stream2.X)

        mw1 = cb.mwmix(stream1.X)
        mw2 = cb.mwmix(stream2.X)

        n1 = stream1.mdot / (mw1 / 1000.0)
        n2 = stream2.mdot / (mw2 / 1000.0)

        # Compare molar enthalpy conservation [J/mol]
        # Expected: H_mixed = (n1*H1 + n2*H2) / (n1 + n2) where n = mdot/MW
        # Actual: H(T_mixed, X_mixed) from CombAero
        # Should match within polynomial precision
        h_mixed_expected = (n1 * h1_cb + n2 * h2_cb) / (n1 + n2)

        mixed_cb = cb.mix([stream1, stream2])
        h_mixed_actual = cb.h(mixed_cb.T, mixed_cb.X)

        rel_diff = abs(h_mixed_actual - h_mixed_expected) / abs(h_mixed_expected)
        assert rel_diff < tolerance_config["enthalpy"], (
            f"Enthalpy [J/mol]: Expected {h_mixed_expected:.1f}, "
            f"got {h_mixed_actual:.1f}, diff={rel_diff * 100:.2f}%"
        )


class TestDensityCalculation:
    """Validate density calculations against Cantera."""

    @pytest.fixture
    def gri30_gas(self, cantera):
        """Create Cantera gas object with GRI-Mech 3.0."""
        return cantera.Solution("gri30.yaml")

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
                ct_species[ct_name] = cb_X[i]
        gas.X = ct_species

    def test_air_density(self, combaero, cantera, gri30_gas, species_mapping, tolerance_config):
        """Test air density at standard conditions."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        T = 300.0
        P = 101325.0

        # Compare density [kg/m^3]
        # Ideal gas law: rho = P*MW / (R*T)
        # CombAero: Uses mixture molecular weight from mole fractions
        # Cantera: Same calculation with same ideal gas assumption
        # Expected differences: < 0.1% (only from MW calculation precision)
        rho_cb = cb.density(T, P, X_air)

        self.set_cantera_composition(gri30_gas, X_air, species_mapping)
        gri30_gas.TP = T, P
        rho_ct = gri30_gas.density

        rel_diff = abs(rho_cb - rho_ct) / rho_ct
        assert rel_diff < tolerance_config["density"], (
            f"Density [kg/m³]: CombAero={rho_cb:.4f}, Cantera={rho_ct:.4f}, diff={rel_diff * 100:.2f}%"
        )

    def test_methane_density(self, combaero, cantera, gri30_gas, species_mapping, tolerance_config):
        """Test pure methane density."""
        cb = combaero

        X_fuel = np.zeros(num_species())
        X_fuel[species_index_from_name("CH4")] = 1.0
        T = 300.0
        P = 101325.0

        rho_cb = cb.density(T, P, X_fuel)

        self.set_cantera_composition(gri30_gas, X_fuel, species_mapping)
        gri30_gas.TP = T, P
        rho_ct = gri30_gas.density

        rel_diff = abs(rho_cb - rho_ct) / rho_ct
        assert rel_diff < tolerance_config["density"], (
            f"Density mismatch: CombAero={rho_cb:.4f} kg/m³, Cantera={rho_ct:.4f} kg/m³"
        )

    def test_density_temperature_variation(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test density at various temperatures."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        P = 101325.0

        for T in [300.0, 500.0, 1000.0, 1500.0, 2000.0]:
            rho_cb = cb.density(T, P, X_air)

            self.set_cantera_composition(gri30_gas, X_air, species_mapping)
            gri30_gas.TP = T, P
            rho_ct = gri30_gas.density

            rel_diff = abs(rho_cb - rho_ct) / rho_ct
            assert rel_diff < tolerance_config["density"], (
                f"T={T} K: CombAero={rho_cb:.4f} kg/m³, Cantera={rho_ct:.4f} kg/m³"
            )

    def test_density_pressure_variation(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test density at various pressures."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        T = 300.0

        for P in [101325.0, 200000.0, 500000.0, 1000000.0]:
            rho_cb = cb.density(T, P, X_air)

            self.set_cantera_composition(gri30_gas, X_air, species_mapping)
            gri30_gas.TP = T, P
            rho_ct = gri30_gas.density

            rel_diff = abs(rho_cb - rho_ct) / rho_ct
            assert rel_diff < tolerance_config["density"], (
                f"P={P / 1e5:.1f} bar: CombAero={rho_cb:.4f} kg/m³, Cantera={rho_ct:.4f} kg/m³"
            )
