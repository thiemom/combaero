#include <gtest/gtest.h>
#include "../include/thermo.h"
#include "../include/transport.h"
#include "../include/combustion.h"
#include "../include/equilibrium.h"
#include <vector>
#include <cmath>

// Test fixture for thermo transport tests
class ThermoTransportTest : public ::testing::Test {
protected:
    // Set up common test data
    void SetUp() override {
        size_t n_species = species_names.size();
        
        // Create generic test vectors sized for current species list
        air_composition.resize(n_species, 0.0);
        non_normalized.resize(n_species, 0.0);
        humid_air.resize(n_species, 0.0);
        water_vapor.resize(n_species, 0.0);
        all_zeros.resize(n_species, 0.0);
        
        // Set up a simple air-like composition using species that should exist
        // N2: ~78%, O2: ~21%, trace others
        try {
            air_composition[species_index_from_name("N2")] = 0.78;
            air_composition[species_index_from_name("O2")] = 0.21;
            if (species_index.count("AR") || species_index.count("Ar")) {
                const std::size_t ar_idx = species_index.count("AR") ? species_index_from_name("AR") : species_index_from_name("Ar");
                air_composition[ar_idx] = 0.01;
            }
        } catch (...) {
            // If species don't exist, just use zeros
        }
        
        // Non-normalized: arbitrary non-zero values
        if (n_species > 0) non_normalized[0] = 0.1;
        if (n_species > 1) non_normalized[1] = 0.4;
        if (n_species > 2) non_normalized[2] = 0.05;
        
        // Humid air: similar to air but with water
        humid_air = air_composition;
        try {
            const std::size_t h2o_idx = species_index_from_name("H2O");
            // Reduce others by 4% and add 4% water
            for (auto& val : humid_air) val *= 0.96;
            humid_air[h2o_idx] = 0.04;
        } catch (...) {
            // If H2O doesn't exist, just use air composition
        }
        
        // Pure water vapor
        try {
            std::size_t h2o_idx = species_index_from_name("H2O");
            water_vapor[h2o_idx] = 1.0;
        } catch (...) {
            // If H2O doesn't exist, leave as zeros
        }
    }
    
    std::vector<double> air_composition;
    std::vector<double> non_normalized;
    std::vector<double> humid_air;
    std::vector<double> water_vapor;
    std::vector<double> all_zeros;
    
    // Helper function to check if two vectors are approximately equal
    bool vectors_approx_equal(const std::vector<double>& a, const std::vector<double>& b, double tolerance = 1e-6) {
        if (a.size() != b.size()) return false;
        
        for (size_t i = 0; i < a.size(); ++i) {
            if (std::abs(a[i] - b[i]) > tolerance) return false;
        }
        
        return true;
    }
    
    // Helper function to check if sum of vector elements is approximately equal to a value
    bool sum_approx_equal(const std::vector<double>& vec, double value, double tolerance = 1e-6) {
        double sum = 0.0;
        for (double x : vec) sum += x;
        return std::abs(sum - value) < tolerance;
    }
};

// Test normalize_fractions function with normal input
TEST_F(ThermoTransportTest, NormalizeNormalInput) {
    auto result = normalize_fractions(non_normalized);
    
    // Check that sum is 1.0
    EXPECT_TRUE(sum_approx_equal(result, 1.0));
    
    // Check that relative proportions are preserved
    for (size_t i = 0; i < non_normalized.size(); ++i) {
        if (std::abs(non_normalized[i]) > 1e-10) {
            double ratio1 = non_normalized[i] / non_normalized[1]; // Compare to O2
            double ratio2 = result[i] / result[1];                // Compare to O2
            EXPECT_NEAR(ratio1, ratio2, 1e-6);
        }
    }
}

// Test normalize_fractions function with already normalized input
TEST_F(ThermoTransportTest, NormalizeNormalizedInput) {
    auto result = normalize_fractions(air_composition);
    
    // Should be approximately the same as input
    EXPECT_TRUE(vectors_approx_equal(result, air_composition));
    
    // Sum should still be 1.0
    EXPECT_TRUE(sum_approx_equal(result, 1.0));
}

// Test normalize_fractions function with all zeros input
TEST_F(ThermoTransportTest, NormalizeAllZeros) {
    // Redirect cerr to capture warning
    testing::internal::CaptureStderr();
    
    auto result = normalize_fractions(all_zeros);
    
    // Check that warning was issued
    std::string output = testing::internal::GetCapturedStderr();
    EXPECT_TRUE(output.find("Warning") != std::string::npos);
    
    // Result should be all zeros
    EXPECT_TRUE(vectors_approx_equal(result, all_zeros));
}

// Test convert_to_dry_fractions function with humid air
TEST_F(ThermoTransportTest, ConvertToDryFractions) {
    auto result = convert_to_dry_fractions(humid_air);
    
    // Check that water vapor is zero (if H2O exists)
    try {
        std::size_t h2o_idx = species_index_from_name("H2O");
        EXPECT_DOUBLE_EQ(result[h2o_idx], 0.0);
    } catch (...) {
        // H2O doesn't exist in species list, skip this check
    }
    
    // Check that sum is 1.0
    EXPECT_TRUE(sum_approx_equal(result, 1.0));
    
    // Check that sum is 1.0
    EXPECT_TRUE(sum_approx_equal(result, 1.0));
}

// Test convert_to_dry_fractions function with pure water vapor
TEST_F(ThermoTransportTest, ConvertPureWaterVaporToDry) {
    // Redirect cerr to capture warning
    testing::internal::CaptureStderr();
    
    auto result = convert_to_dry_fractions(water_vapor);
    
    // Check that warning was issued
    std::string output = testing::internal::GetCapturedStderr();
    EXPECT_TRUE(output.find("Warning") != std::string::npos);
    
    // Result should be all zeros
    EXPECT_TRUE(vectors_approx_equal(result, all_zeros));
}

// Test basic thermodynamic properties
TEST_F(ThermoTransportTest, BasicThermodynamicProperties) {
    // Test at standard conditions (300 K to avoid boundary warnings)
    double T = 300.0;
    double P = 101325.0;
    
    // Test specific heat capacity - should be positive and reasonable
    double cp_value = cp(T, air_composition);
    EXPECT_GT(cp_value, 20.0);  // Reasonable lower bound
    EXPECT_LT(cp_value, 50.0);  // Reasonable upper bound
    
    // Test enthalpy - should be finite
    double h_value = h(T, air_composition);
    EXPECT_TRUE(std::isfinite(h_value));
    
    // Test entropy - should be positive and finite
    double s_value = s(T, air_composition, P);
    EXPECT_GT(s_value, 100.0);
    EXPECT_TRUE(std::isfinite(s_value));
    
    // Test density - should be positive and reasonable for a gas
    double rho = density(T, P, air_composition);
    EXPECT_GT(rho, 0.5);  // Reasonable for any gas mixture
    EXPECT_LT(rho, 5.0);  // Reasonable upper bound
    
    // Test internal energy: u = h - R*T for ideal gas
    double u_value = u(T, air_composition);
    double expected_u = h_value - 8.31446261815324 * T;
    EXPECT_NEAR(u_value, expected_u, 1.0e-6);
}

// Test molecular weight calculation
TEST_F(ThermoTransportTest, MolecularWeight) {
    // Molecular weight should be positive and reasonable
    double mw = mwmix(air_composition);
    EXPECT_GT(mw, 10.0);  // Lighter than any common gas mixture
    EXPECT_LT(mw, 100.0); // Heavier than most common mixtures
}

// Test per-species dimensional properties
TEST_F(ThermoTransportTest, PerSpeciesProperties) {
    double T = 300.0;
    std::size_t i_N2 = species_index_from_name("N2");
    
    // cp_species = cp_R * R
    double cp_dim = cp_species(i_N2, T);
    double cp_nondim = cp_R(i_N2, T);
    EXPECT_NEAR(cp_dim, cp_nondim * 8.31446261815324, 1.0e-10);
    
    // h_species = h_RT * R * T
    double h_dim = h_species(i_N2, T);
    double h_nondim = h_RT(i_N2, T);
    EXPECT_NEAR(h_dim, h_nondim * 8.31446261815324 * T, 1.0e-6);
    
    // s_species = s_R * R
    double s_dim = s_species(i_N2, T);
    double s_nondim = s_R(i_N2, T);
    EXPECT_NEAR(s_dim, s_nondim * 8.31446261815324, 1.0e-10);
}

// Test transport properties
TEST_F(ThermoTransportTest, TransportProperties) {
    // Test at 300 K to avoid boundary warnings
    double T = 300.0;
    double P = 101325.0;
    
    // Test viscosity - should be positive and reasonable for gases
    double mu = viscosity(T, P, air_composition);
    EXPECT_GT(mu, 1.0e-6);  // Lower bound for gas viscosity
    EXPECT_LT(mu, 1.0e-4);  // Upper bound for gas viscosity
    
    // Test thermal diffusivity: α = k / (ρ * cp)
    double alpha = thermal_diffusivity(T, P, air_composition);
    double k_val = thermal_conductivity(T, P, air_composition);
    double rho_val = density(T, P, air_composition);
    double cp_val = cp(T, air_composition);
    double MW = mwmix(air_composition) / 1000.0;  // kg/mol
    double cp_mass = cp_val / MW;
    double expected_alpha = k_val / (rho_val * cp_mass);
    EXPECT_NEAR(alpha, expected_alpha, 1.0e-12);
    
    // Test thermal conductivity - should be positive
    double k = thermal_conductivity(T, P, air_composition);
    EXPECT_GT(k, 0.001);
    EXPECT_LT(k, 10.0);  // Wide range for different species mixtures
    
    // Test Prandtl number - should be positive
    double pr = prandtl(T, P, air_composition);
    EXPECT_GT(pr, 1e-6);  // Just check it's positive and finite
    EXPECT_LT(pr, 100.0);
}

// Helper to build a CH4/O2 mixture (other species zero)
static std::vector<double> make_CH4_O2_mixture(double n_CH4, double n_O2) {
    std::vector<double> X(species_names.size(), 0.0);

    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2  = species_index_from_name("O2");

    double n_tot = n_CH4 + n_O2;
    X[idx_CH4] = n_CH4 / n_tot;
    X[idx_O2]  = n_O2  / n_tot;

    return X;
}

// Equivalence-ratio consistency test: multi-species fuel + oxidizer streams
TEST_F(ThermoTransportTest, EquivalenceRatioConsistency) {
    const std::size_t n = species_names.size();

    // Fuel stream: pure CH4
    std::vector<double> X_fuel(n, 0.0);
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    X_fuel[idx_CH4] = 1.0;

    // Oxidizer stream: simple N2/O2 "air" (79% N2, 21% O2)
    std::vector<double> X_ox(n, 0.0);
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    X_ox[idx_N2] = 0.79;
    X_ox[idx_O2] = 0.21;

    // A few representative equivalence ratios
    const double phis[] = {0.5, 1.0, 2.0};

    for (double phi_target : phis) {
        auto X_mix = set_equivalence_ratio_mole(phi_target, X_fuel, X_ox);

        // Mixture should be normalized
        double sum = 0.0;
        for (double x : X_mix) sum += x;
        EXPECT_NEAR(sum, 1.0, 1e-12);

        double phi_back = equivalence_ratio_mole(X_mix, X_fuel, X_ox);
        EXPECT_NEAR(phi_back, phi_target, 1e-10);
    }
}

// Mass-basis equivalence-ratio consistency test: mirrors mole-based test
TEST_F(ThermoTransportTest, EquivalenceRatioMassConsistency) {
    const std::size_t n = species_names.size();

    // Fuel stream (mole basis): pure CH4
    std::vector<double> X_fuel(n, 0.0);
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    X_fuel[idx_CH4] = 1.0;

    // Oxidizer stream (mole basis): simple N2/O2 "air" (79% N2, 21% O2)
    std::vector<double> X_ox(n, 0.0);
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    X_ox[idx_N2] = 0.79;
    X_ox[idx_O2] = 0.21;

    // Convert to mass fractions for mass-basis helpers
    std::vector<double> Y_fuel = mole_to_mass(X_fuel);
    std::vector<double> Y_ox   = mole_to_mass(X_ox);

    const double phis[] = {0.5, 1.0, 2.0};

    for (double phi_target : phis) {
        auto Y_mix = set_equivalence_ratio_mass(phi_target, Y_fuel, Y_ox);

        // Mixture should be normalized
        double sum = 0.0;
        for (double y : Y_mix) sum += y;
        EXPECT_NEAR(sum, 1.0, 1e-12);

        double phi_back = equivalence_ratio_mass(Y_mix, Y_fuel, Y_ox);
        EXPECT_NEAR(phi_back, phi_target, 1e-10);
    }
}

// Bilger Z <-> phi (mass basis) consistency test
TEST_F(ThermoTransportTest, BilgerZPhiMassConsistency) {
    const std::size_t n = species_names.size();

    // Fuel stream (mole basis): pure CH4
    std::vector<double> X_fuel(n, 0.0);
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    X_fuel[idx_CH4] = 1.0;

    // Oxidizer stream (mole basis): simple N2/O2 "air" (79% N2, 21% O2)
    std::vector<double> X_ox(n, 0.0);
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    X_ox[idx_N2] = 0.79;
    X_ox[idx_O2] = 0.21;

    // Convert to mass fractions
    std::vector<double> Y_fuel = mole_to_mass(X_fuel);
    std::vector<double> Y_ox   = mole_to_mass(X_ox);

    // Round-trip test for several phi values
    const double phis[] = {0.5, 1.0, 2.0};

    for (double phi_target : phis) {
        const double Z = bilger_Z_from_equivalence_ratio_mass(phi_target, Y_fuel, Y_ox);
        EXPECT_GT(Z, 0.0);
        EXPECT_LT(Z, 1.0);

        const double phi_back = equivalence_ratio_from_bilger_Z_mass(Z, Y_fuel, Y_ox);
        EXPECT_NEAR(phi_back, phi_target, 1e-10);
    }

    // Check that Z_st corresponds to phi ~ 1
    const double Z_st = bilger_stoich_mixture_fraction_mass(Y_fuel, Y_ox);
    EXPECT_GT(Z_st, 0.0);
    EXPECT_LT(Z_st, 1.0);

    const double phi_from_Z_st = equivalence_ratio_from_bilger_Z_mass(Z_st, Y_fuel, Y_ox);
    EXPECT_NEAR(phi_from_Z_st, 1.0, 1e-10);
}

// No O2: mixture should be unchanged, f = 0
TEST_F(ThermoTransportTest, Combustion_NoOxygen) {
    auto X_in = make_CH4_O2_mixture(1.0, 0.0);

    double f = -1.0;
    auto X_out = complete_combustion_to_CO2_H2O(X_in, f);

    EXPECT_NEAR(f, 0.0, 1e-12);
    EXPECT_TRUE(vectors_approx_equal(X_in, X_out));
}

// Exactly stoichiometric O2: CH4 + 2 O2 -> CO2 + 2 H2O, f = 1, no O2 remaining
TEST_F(ThermoTransportTest, Combustion_StoichiometricOxygen) {
    // 1 mol CH4, 2 mol O2
    auto X_in = make_CH4_O2_mixture(1.0, 2.0);

    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_O2  = species_index_from_name("O2");
    const std::size_t idx_CH4 = species_index_from_name("CH4");

    double f = -1.0;
    auto X_out = complete_combustion_to_CO2_H2O(X_in, f);

    EXPECT_NEAR(f, 1.0, 1e-12);

    // Final moles: 1 CO2 + 2 H2O = 3 mol
    // Mole fractions: CO2 = 1/3, H2O = 2/3
    EXPECT_NEAR(X_out[idx_CO2], 1.0 / 3.0, 1e-8);
    EXPECT_NEAR(X_out[idx_H2O], 2.0 / 3.0, 1e-8);
    EXPECT_NEAR(X_out[idx_O2], 0.0, 1e-12);
    EXPECT_NEAR(X_out[idx_CH4], 0.0, 1e-12);
}

// Excess O2 (strongly lean): CH4 + 3 O2, f = 1, O2 remaining
TEST_F(ThermoTransportTest, Combustion_ExcessOxygenStronglyLean) {
    auto X_in = make_CH4_O2_mixture(1.0, 3.0);

    std::size_t idx_CO2 = species_index_from_name("CO2");
    std::size_t idx_H2O = species_index_from_name("H2O");
    std::size_t idx_O2  = species_index_from_name("O2");
    std::size_t idx_CH4 = species_index_from_name("CH4");

    double f = -1.0;
    auto X_out = complete_combustion_to_CO2_H2O(X_in, f);

    EXPECT_NEAR(f, 1.0, 1e-12);

    // After reaction: 1 CO2 + 2 H2O + 1 O2 = 4 mol
    // Mole fractions: O2 = 1/4, CO2 = 1/4, H2O = 1/2
    EXPECT_NEAR(X_out[idx_O2],  0.25, 1e-8);
    EXPECT_NEAR(X_out[idx_CO2], 0.25, 1e-8);
    EXPECT_NEAR(X_out[idx_H2O], 0.50, 1e-8);
    EXPECT_NEAR(X_out[idx_CH4], 0.0,  1e-12);
}

// Intermediate lean case: CH4 + 2.5 O2 (still fuel-limited), f = 1
TEST_F(ThermoTransportTest, Combustion_IntermediateLean) {
    auto X_in = make_CH4_O2_mixture(1.0, 2.5);

    std::size_t idx_CO2 = species_index_from_name("CO2");
    std::size_t idx_H2O = species_index_from_name("H2O");
    std::size_t idx_O2  = species_index_from_name("O2");
    std::size_t idx_CH4 = species_index_from_name("CH4");

    double f = -1.0;
    auto X_out = complete_combustion_to_CO2_H2O(X_in, f);

    EXPECT_NEAR(f, 1.0, 1e-12);

    // After reaction: 1 CO2 + 2 H2O + 0.5 O2 = 3.5 mol
    // Mole fractions: O2 = 0.5/3.5, CO2 = 1/3.5, H2O = 2/3.5
    double denom = 3.5;
    EXPECT_NEAR(X_out[idx_O2],  0.5 / denom, 1e-8);
    EXPECT_NEAR(X_out[idx_CO2], 1.0 / denom, 1e-8);
    EXPECT_NEAR(X_out[idx_H2O], 2.0 / denom, 1e-8);
    EXPECT_NEAR(X_out[idx_CH4], 0.0,         1e-12);
}

// Intermediate rich case: CH4 + 1.5 O2 (O2-limited), 0 < f < 1, CH4 remains
TEST_F(ThermoTransportTest, Combustion_IntermediateRich) {
    auto X_in = make_CH4_O2_mixture(1.0, 1.5);

    std::size_t idx_CO2 = species_index_from_name("CO2");
    std::size_t idx_H2O = species_index_from_name("H2O");
    std::size_t idx_O2  = species_index_from_name("O2");
    std::size_t idx_CH4 = species_index_from_name("CH4");

    double f = -1.0;
    auto X_out = complete_combustion_to_CO2_H2O(X_in, f);

    // Stoichiometric O2 requirement is 2 mol per mol CH4
    // Here we have 1.5 mol O2 -> f = 1.5 / 2 = 0.75
    EXPECT_NEAR(f, 0.75, 1e-12);

    // Reacted CH4: 0.75 mol, remaining CH4: 0.25 mol
    // Products: 0.75 CO2, 1.5 H2O, all O2 consumed
    // Total moles after reaction: 0.25 CH4 + 0.75 CO2 + 1.5 H2O = 2.5
    double denom = 2.5;
    EXPECT_NEAR(X_out[idx_CH4], 0.25 / denom, 1e-8);
    EXPECT_NEAR(X_out[idx_CO2], 0.75 / denom, 1e-8);
    EXPECT_NEAR(X_out[idx_H2O], 1.50 / denom, 1e-8);
    EXPECT_NEAR(X_out[idx_O2],  0.0,          1e-12);
}

// Test State-based thermo functions
TEST_F(ThermoTransportTest, StateBasedThermo) {
    State s;
    s.T = 300.0;
    s.P = 101325.0;
    s.X = air_composition;
    
    // State-based functions should match vector-based functions
    EXPECT_NEAR(cp(s), cp(s.T, s.X), 1e-12);
    EXPECT_NEAR(h(s), h(s.T, s.X), 1e-12);
    EXPECT_NEAR(::s(s), ::s(s.T, s.X, s.P), 1e-12);
    EXPECT_NEAR(cv(s), cv(s.T, s.X), 1e-12);
    EXPECT_NEAR(u(s), u(s.T, s.X), 1e-12);
    EXPECT_NEAR(density(s), density(s.T, s.P, s.X), 1e-12);
    EXPECT_NEAR(mwmix(s), mwmix(s.X), 1e-12);
}

// Test State-based combustion functions
TEST_F(ThermoTransportTest, StateBasedCombustion) {
    const std::size_t n = species_names.size();
    
    // Create stoichiometric CH4 + air mixture
    std::vector<double> X_mix(n, 0.0);
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // CH4 + 2 O2 -> CO2 + 2 H2O, stoichiometric with air (21% O2)
    // For phi=1: need 2 mol O2 per mol CH4
    // Air is 21% O2, so need 2/0.21 = 9.52 mol air per mol CH4
    double n_CH4 = 1.0;
    double n_air = 9.52;
    double n_total = n_CH4 + n_air;
    X_mix[idx_CH4] = n_CH4 / n_total;
    X_mix[idx_O2] = 0.21 * n_air / n_total;
    X_mix[idx_N2] = 0.79 * n_air / n_total;
    
    State in;
    in.T = 300.0;
    in.P = 101325.0;
    in.X = X_mix;
    
    // Isothermal combustion should preserve T
    State out_iso = complete_combustion_isothermal(in);
    EXPECT_NEAR(out_iso.T, in.T, 1e-10);
    EXPECT_NEAR(out_iso.P, in.P, 1e-10);
    
    // Adiabatic combustion should increase T significantly
    State out_ad = complete_combustion(in);
    EXPECT_GT(out_ad.T, in.T + 1000.0);  // Flame temp should be >1300 K
    EXPECT_NEAR(out_ad.P, in.P, 1e-10);
    
    // Enthalpy should be conserved in adiabatic case
    double H_in = h(in);
    double H_out = h(out_ad);
    EXPECT_NEAR(H_in, H_out, 1.0);  // Within 1 J/mol
}

// Test State-based WGS equilibrium functions
TEST_F(ThermoTransportTest, StateBasedWgsEquilibrium) {
    const std::size_t n = species_names.size();
    
    // Create a mixture with CO, H2O, CO2, H2
    std::vector<double> X_mix(n, 0.0);
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2 = species_index_from_name("H2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    X_mix[idx_CO] = 0.1;
    X_mix[idx_H2O] = 0.2;
    X_mix[idx_CO2] = 0.05;
    X_mix[idx_H2] = 0.05;
    X_mix[idx_N2] = 0.6;  // Diluent
    
    State in;
    in.T = 1000.0;  // High T for WGS
    in.P = 101325.0;
    in.X = X_mix;
    
    // Isothermal WGS should preserve T
    State out_iso = wgs_equilibrium(in);
    EXPECT_NEAR(out_iso.T, in.T, 1e-10);
    EXPECT_NEAR(out_iso.P, in.P, 1e-10);
    
    // Mole fractions should still sum to 1
    double sum = 0.0;
    for (double x : out_iso.X) sum += x;
    EXPECT_NEAR(sum, 1.0, 1e-10);
    
    // Adiabatic WGS
    State out_ad = wgs_equilibrium_adiabatic(in);
    EXPECT_NEAR(out_ad.P, in.P, 1e-10);
    
    // Enthalpy should be conserved
    double H_in = h(in);
    double H_out = h(out_ad);
    EXPECT_NEAR(H_in, H_out, 1.0);  // Within 1 J/mol
}

// Test dry air requirements for combustion
TEST_F(ThermoTransportTest, DryAirRequirements) {
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    
    // CH4 + 2 O2 -> CO2 + 2 H2O
    // O2 required per mol CH4 = 2
    double O2_per_mol = oxygen_required_per_mol_fuel(idx_CH4);
    EXPECT_NEAR(O2_per_mol, 2.0, 1e-10);
    
    // Dry air is ~20.95% O2, so air required = O2_required / 0.2095
    double air_per_mol = dryair_required_per_mol_fuel(idx_CH4);
    double expected_air = O2_per_mol / 0.2095;
    EXPECT_NEAR(air_per_mol, expected_air, 1e-10);
    
    // Test mass basis consistency
    double O2_per_kg = oxygen_required_per_kg_fuel(idx_CH4);
    double air_per_kg = dryair_required_per_kg_fuel(idx_CH4);
    
    // Both should be positive
    EXPECT_GT(O2_per_kg, 0.0);
    EXPECT_GT(air_per_kg, 0.0);
    
    // Air requirement should be larger than O2 requirement (air is ~21% O2)
    EXPECT_GT(air_per_kg, O2_per_kg);
    
    // Test mixture functions with pure CH4
    const std::size_t n = species_names.size();
    std::vector<double> X_fuel(n, 0.0);
    X_fuel[idx_CH4] = 1.0;
    
    double air_per_mol_mix = dryair_required_per_mol_mixture(X_fuel);
    EXPECT_NEAR(air_per_mol_mix, air_per_mol, 1e-10);
    
    double air_per_kg_mix = dryair_required_per_kg_mixture(X_fuel);
    EXPECT_NEAR(air_per_kg_mix, air_per_kg, 1e-10);
}
