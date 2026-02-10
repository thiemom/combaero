import numpy as np
import pytest

from combaero._core import num_species, species_index_from_name


class TestTransportProperties:
    """Validate transport property calculations against Cantera."""

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

    def test_air_viscosity(self, combaero, cantera, gri30_gas, species_mapping, tolerance_config):
        """Test air viscosity at standard conditions."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        T = 300.0
        P = 101325.0

        # Compare dynamic viscosity [Pa*s]
        # CombAero: Sutherland/polynomial correlations with specific L-J parameters
        # Cantera: GRI-Mech 3.0 transport data with different L-J parameters
        # Expected differences: 10-20% due to different correlations
        mu_cb = cb.viscosity(T, P, X_air)

        self.set_cantera_composition(gri30_gas, X_air, species_mapping)
        gri30_gas.TP = T, P
        mu_ct = gri30_gas.viscosity

        rel_diff = abs(mu_cb - mu_ct) / mu_ct
        assert rel_diff < tolerance_config["transport"], (
            f"Viscosity [Pa·s]: CombAero={mu_cb:.2e}, Cantera={mu_ct:.2e}, "  # noqa: E501
            f"diff={rel_diff * 100:.1f}%"
        )

    def test_viscosity_temperature_variation(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test viscosity at various temperatures."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        P = 101325.0

        for T in [300.0, 500.0, 1000.0, 1500.0, 2000.0]:
            mu_cb = cb.viscosity(T, P, X_air)

            self.set_cantera_composition(gri30_gas, X_air, species_mapping)
            gri30_gas.TP = T, P
            mu_ct = gri30_gas.viscosity

            rel_diff = abs(mu_cb - mu_ct) / mu_ct
            assert rel_diff < tolerance_config["transport"], (
                f"T={T} K: CombAero={mu_cb:.2e} Pa*s, "
                f"Cantera={mu_ct:.2e} Pa*s, diff={rel_diff * 100:.1f}%"
            )

    def test_methane_viscosity(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test pure methane viscosity."""
        cb = combaero

        X_fuel = np.zeros(num_species())
        X_fuel[species_index_from_name("CH4")] = 1.0
        T = 300.0
        P = 101325.0

        mu_cb = cb.viscosity(T, P, X_fuel)

        self.set_cantera_composition(gri30_gas, X_fuel, species_mapping)
        gri30_gas.TP = T, P
        mu_ct = gri30_gas.viscosity

        rel_diff = abs(mu_cb - mu_ct) / mu_ct
        assert rel_diff < tolerance_config["transport"], (
            f"Viscosity mismatch: CombAero={mu_cb:.2e} Pa·s, Cantera={mu_ct:.2e} Pa·s"
        )

    def test_air_thermal_conductivity(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test air thermal conductivity at standard conditions."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        T = 300.0
        P = 101325.0

        # Compare thermal conductivity [W/(m*K)]
        # CombAero: Eucken or modified Eucken formula with specific transport data
        # Cantera: GRI-Mech 3.0 transport data with different mixing rules
        # Expected differences: 15-20% due to different correlations and mixing rules
        k_cb = cb.thermal_conductivity(T, P, X_air)

        self.set_cantera_composition(gri30_gas, X_air, species_mapping)
        gri30_gas.TP = T, P
        k_ct = gri30_gas.thermal_conductivity

        rel_diff = abs(k_cb - k_ct) / k_ct
        assert rel_diff < tolerance_config["transport"], (
            f"Thermal conductivity [W/(m*K)]: CombAero={k_cb:.4f}, "
            f"Cantera={k_ct:.4f}, diff={rel_diff * 100:.1f}%"
        )

    def test_thermal_conductivity_temperature_variation(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test thermal conductivity at various temperatures."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        P = 101325.0

        for T in [300.0, 500.0, 1000.0, 1500.0, 2000.0]:
            k_cb = cb.thermal_conductivity(T, P, X_air)

            self.set_cantera_composition(gri30_gas, X_air, species_mapping)
            gri30_gas.TP = T, P
            k_ct = gri30_gas.thermal_conductivity

            rel_diff = abs(k_cb - k_ct) / k_ct
            assert rel_diff < tolerance_config["transport"], (
                f"T={T} K: CombAero={k_cb:.4f} W/(m*K), "
                f"Cantera={k_ct:.4f} W/(m*K), diff={rel_diff * 100:.1f}%"
            )

    def test_methane_thermal_conductivity(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test pure methane thermal conductivity."""
        cb = combaero

        X_fuel = np.zeros(num_species())
        X_fuel[species_index_from_name("CH4")] = 1.0
        T = 300.0
        P = 101325.0

        k_cb = cb.thermal_conductivity(T, P, X_fuel)

        self.set_cantera_composition(gri30_gas, X_fuel, species_mapping)
        gri30_gas.TP = T, P
        k_ct = gri30_gas.thermal_conductivity

        rel_diff = abs(k_cb - k_ct) / k_ct
        assert rel_diff < tolerance_config["transport"], (
            f"Thermal conductivity mismatch: CombAero={k_cb:.4f} W/(m·K), Cantera={k_ct:.4f} W/(m·K)"
        )

    def test_prandtl_number(self, combaero, cantera, gri30_gas, species_mapping, tolerance_config):
        """Test Prandtl number calculation."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        T = 300.0
        P = 101325.0

        # Compare Prandtl number [dimensionless]
        # Pr = (Cp * μ) / k
        # CombAero: Uses CombAero's Cp, μ, k
        # Cantera: Uses Cantera's Cp, μ, k
        # Differences propagate from transport property differences
        Pr_cb = cb.prandtl(T, P, X_air)

        self.set_cantera_composition(gri30_gas, X_air, species_mapping)
        gri30_gas.TP = T, P

        cp_ct = gri30_gas.cp_mass  # [J/(kg*K)]
        mu_ct = gri30_gas.viscosity  # [Pa*s]
        k_ct = gri30_gas.thermal_conductivity  # [W/(m*K)]
        Pr_ct = cp_ct * mu_ct / k_ct  # [dimensionless]

        rel_diff = abs(Pr_cb - Pr_ct) / Pr_ct
        assert rel_diff < tolerance_config["transport"], (
            f"Prandtl number [dimensionless]: CombAero={Pr_cb:.4f}, "
            f"Cantera={Pr_ct:.4f}, diff={rel_diff * 100:.1f}%"
        )

    def test_mixture_transport_properties(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test transport properties for fuel-air mixture."""
        cb = combaero
        from combaero._core import species_index_from_name

        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0

        phi = 0.80
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_air)

        T = 300.0
        P = 101325.0

        mu_cb = cb.viscosity(T, P, X_mix)
        k_cb = cb.thermal_conductivity(T, P, X_mix)

        self.set_cantera_composition(gri30_gas, X_mix, species_mapping)
        gri30_gas.TP = T, P
        mu_ct = gri30_gas.viscosity
        k_ct = gri30_gas.thermal_conductivity

        mu_rel_diff = abs(mu_cb - mu_ct) / mu_ct
        k_rel_diff = abs(k_cb - k_ct) / k_ct

        assert mu_rel_diff < tolerance_config["transport"], (
            f"Viscosity mismatch: CombAero={mu_cb:.2e} Pa·s, Cantera={mu_ct:.2e} Pa·s"
        )
        assert k_rel_diff < tolerance_config["transport"], (
            f"Thermal conductivity mismatch: CombAero={k_cb:.4f} W/(m·K), Cantera={k_ct:.4f} W/(m·K)"  # noqa: E501
        )

    def test_high_temperature_transport(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test transport properties at high temperatures (post-combustion)."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        X_fuel = np.zeros(len(X_air))
        X_fuel[species_index_from_name("CH4")] = 1.0

        phi = 1.0
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_air)

        burned = cb.complete_combustion(T=300.0, X=X_mix, P=101325.0)

        mu_cb = cb.viscosity(burned.T, burned.P, burned.X)
        k_cb = cb.thermal_conductivity(burned.T, burned.P, burned.X)

        self.set_cantera_composition(gri30_gas, burned.X, species_mapping)
        gri30_gas.TP = burned.T, burned.P
        mu_ct = gri30_gas.viscosity
        k_ct = gri30_gas.thermal_conductivity

        mu_rel_diff = abs(mu_cb - mu_ct) / mu_ct
        k_rel_diff = abs(k_cb - k_ct) / k_ct

        assert mu_rel_diff < tolerance_config["transport"], (
            f"High-T viscosity mismatch at {burned.T:.0f} K: "
            f"CombAero={mu_cb:.2e} Pa*s, Cantera={mu_ct:.2e} Pa*s"
        )
        assert k_rel_diff < tolerance_config["transport"], (
            f"High-T thermal conductivity mismatch at {burned.T:.0f} K: "
            f"CombAero={k_cb:.4f} W/(m*K), Cantera={k_ct:.4f} W/(m*K)"
        )


class TestThermodynamicProperties:
    """Validate thermodynamic property calculations against Cantera."""

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

    def test_air_cp(self, combaero, cantera, gri30_gas, species_mapping, tolerance_config):
        """Test air specific heat capacity."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        T = 300.0
        P = 101325.0

        # Compare molar heat capacity at constant pressure [J/(mol*K)]
        # CombAero: NASA-9 polynomials (7 coefficients, 200-6000 K)
        # Cantera: NASA-7 polynomials (7 coefficients per range, with break point)
        # Expected differences: < 1% for well-fitted polynomials
        cp_cb = cb.cp(T, X_air)  # [J/(mol*K)]

        self.set_cantera_composition(gri30_gas, X_air, species_mapping)
        gri30_gas.TP = T, P
        cp_ct = gri30_gas.cp_mole / 1000.0  # Convert [J/(kmol*K)] to [J/(mol*K)]

        rel_diff = abs(cp_cb - cp_ct) / cp_ct
        assert rel_diff < tolerance_config["enthalpy"], (
            f"Cp [J/(mol·K)]: CombAero={cp_cb:.2f}, Cantera={cp_ct:.2f}, diff={rel_diff * 100:.2f}%"
        )

    def test_cp_temperature_variation(
        self, combaero, cantera, gri30_gas, species_mapping, tolerance_config
    ):
        """Test Cp at various temperatures."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        P = 101325.0

        for T in [300.0, 500.0, 1000.0, 1500.0, 2000.0]:
            cp_cb = cb.cp(T, X_air)

            self.set_cantera_composition(gri30_gas, X_air, species_mapping)
            gri30_gas.TP = T, P
            cp_ct = gri30_gas.cp_mole / 1000.0  # Convert J/(kmol*K) to J/(mol*K)

            rel_diff = abs(cp_cb - cp_ct) / cp_ct
            assert rel_diff < tolerance_config["enthalpy"], (
                f"T={T} K: CombAero={cp_cb:.2f} J/(mol·K), Cantera={cp_ct:.2f} J/(mol·K)"
            )

    def test_enthalpy(self, combaero, cantera, gri30_gas, species_mapping, tolerance_config):
        """Test enthalpy calculation."""
        cb = combaero

        X_air = cb.standard_dry_air_composition()
        T = 300.0
        P = 101325.0

        # Compare molar enthalpy [J/mol]
        # H(T) = H_ref(298.15K) + integral(Cp(T)dT) from 298.15K to T
        # CombAero: NASA-9 polynomials integrated
        # Cantera: NASA-7 polynomials integrated
        # Expected differences: < 1.5% due to polynomial differences
        h_cb = cb.h(T, X_air)  # [J/mol]

        self.set_cantera_composition(gri30_gas, X_air, species_mapping)
        gri30_gas.TP = T, P
        h_ct = gri30_gas.enthalpy_mole / 1000.0  # Convert [J/kmol] to [J/mol]

        rel_diff = abs(h_cb - h_ct) / abs(h_ct)
        assert rel_diff < tolerance_config["enthalpy"], (
            f"Enthalpy [J/mol]: CombAero={h_cb:.1f}, Cantera={h_ct:.1f}, diff={rel_diff * 100:.2f}%"
        )
