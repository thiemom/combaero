#include <gtest/gtest.h>
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"
#include "../include/transport.h"
#include "../include/combustion.h"
#include "../include/equilibrium.h"
#include "../include/humidair.h"
#include "../include/state.h"
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

// Test that temperature below valid range issues a warning but still returns a result
TEST_F(ThermoTransportTest, TemperatureBelowValidRangeWarning) {
    // Use 298.15 K which is below the 300 K minimum for most species
    double T = 298.15;
    
    // Redirect cerr to capture warning
    testing::internal::CaptureStderr();
    
    // Calculate cp - should still work but issue warning
    double cp_value = cp(T, air_composition);
    
    // Check that warning was issued
    std::string output = testing::internal::GetCapturedStderr();
    EXPECT_TRUE(output.find("Warning") != std::string::npos) 
        << "Expected temperature range warning for T=298.15 K";
    EXPECT_TRUE(output.find("below valid range") != std::string::npos)
        << "Warning should mention 'below valid range'";
    
    // Result should still be reasonable (extrapolated)
    EXPECT_GT(cp_value, 20.0);
    EXPECT_LT(cp_value, 50.0);
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

// Test State property getters
TEST_F(ThermoTransportTest, StatePropertyGetters) {
    State s;
    s.T = 300.0;
    s.P = 101325.0;
    s.X = air_composition;
    
    // Test property getters match free functions
    EXPECT_NEAR(s.mw(), mwmix(s.X), 1e-12);
    EXPECT_NEAR(s.cp(), cp(s.T, s.X), 1e-12);
    EXPECT_NEAR(s.h(), h(s.T, s.X), 1e-12);
    EXPECT_NEAR(s.s(), ::s(s.T, s.X, s.P), 1e-12);
    EXPECT_NEAR(s.cv(), cv(s.T, s.X), 1e-12);
    EXPECT_NEAR(s.rho(), density(s.T, s.P, s.X), 1e-12);
    EXPECT_NEAR(s.gamma(), isentropic_expansion_coefficient(s.T, s.X), 1e-12);
    EXPECT_NEAR(s.a(), speed_of_sound(s.T, s.X), 1e-12);
    
    // Test transport properties
    EXPECT_NEAR(s.mu(), viscosity(s.T, s.P, s.X), 1e-12);
    EXPECT_NEAR(s.k(), thermal_conductivity(s.T, s.P, s.X), 1e-12);
    EXPECT_NEAR(s.nu(), kinematic_viscosity(s.T, s.P, s.X), 1e-12);
    EXPECT_NEAR(s.Pr(), prandtl(s.T, s.P, s.X), 1e-12);
    EXPECT_NEAR(s.alpha(), thermal_diffusivity(s.T, s.P, s.X), 1e-12);
    
    // Test setters with chaining
    State s2;
    s2.set_T(400.0).set_P(200000.0).set_X(air_composition);
    EXPECT_NEAR(s2.T, 400.0, 1e-12);
    EXPECT_NEAR(s2.P, 200000.0, 1e-12);
}

// Test Stream mixing
TEST_F(ThermoTransportTest, StreamMixing) {
    const std::size_t n = species_names.size();
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    
    // Stream 1: Hot air at 500 K
    std::vector<double> X_air(n, 0.0);
    X_air[idx_N2] = 0.79;
    X_air[idx_O2] = 0.21;
    
    Stream s1;
    s1.state.T = 500.0;
    s1.state.P = 101325.0;
    s1.state.X = X_air;
    s1.mdot = 1.0;  // 1 kg/s
    
    // Stream 2: Cold CO2 at 300 K
    std::vector<double> X_co2(n, 0.0);
    X_co2[idx_CO2] = 1.0;
    
    Stream s2;
    s2.state.T = 300.0;
    s2.state.P = 90000.0;  // Lower pressure
    s2.state.X = X_co2;
    s2.mdot = 0.5;  // 0.5 kg/s
    
    // Mix streams
    Stream mixed = mix({s1, s2});
    
    // Check mass balance
    EXPECT_NEAR(mixed.mdot, 1.5, 1e-12);
    
    // Check pressure is minimum (90000 Pa)
    EXPECT_NEAR(mixed.P(), 90000.0, 1e-12);
    
    // Check temperature is between inputs (energy balance)
    EXPECT_GT(mixed.T(), 300.0);
    EXPECT_LT(mixed.T(), 500.0);
    
    // Check composition has both species
    EXPECT_GT(mixed.X()[idx_N2], 0.0);
    EXPECT_GT(mixed.X()[idx_O2], 0.0);
    EXPECT_GT(mixed.X()[idx_CO2], 0.0);
    
    // Mole fractions should sum to 1
    double sum = 0.0;
    for (double x : mixed.X()) sum += x;
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

// Test Stream mixing with pressure override
TEST_F(ThermoTransportTest, StreamMixingPressureOverride) {
    const std::size_t n = species_names.size();
    
    Stream s1;
    s1.state.T = 400.0;
    s1.state.P = 100000.0;
    s1.state.X = air_composition;
    s1.mdot = 1.0;
    
    Stream s2;
    s2.state.T = 300.0;
    s2.state.P = 80000.0;
    s2.state.X = air_composition;
    s2.mdot = 1.0;
    
    // Mix with explicit pressure
    double P_override = 150000.0;
    Stream mixed = mix({s1, s2}, P_override);
    
    EXPECT_NEAR(mixed.P(), P_override, 1e-12);
}

// Test Stream mixing enthalpy conservation
TEST_F(ThermoTransportTest, StreamMixingEnthalpyConservation) {
    const std::size_t n = species_names.size();
    
    Stream s1;
    s1.state.T = 600.0;
    s1.state.P = 101325.0;
    s1.state.X = air_composition;
    s1.mdot = 2.0;
    
    Stream s2;
    s2.state.T = 300.0;
    s2.state.P = 101325.0;
    s2.state.X = air_composition;
    s2.mdot = 1.0;
    
    Stream mixed = mix({s1, s2});
    
    // Calculate enthalpy flows (J/s)
    auto H_flow = [](const Stream& st) {
        double h_molar = st.h();  // J/mol
        double MW = st.mw();      // g/mol
        double h_mass = h_molar / (MW * 1e-3);  // J/kg
        return h_mass * st.mdot;  // W
    };
    
    double H_in = H_flow(s1) + H_flow(s2);
    double H_out = H_flow(mixed);
    
    // Enthalpy should be conserved
    EXPECT_NEAR(H_in, H_out, std::abs(H_in) * 1e-6);
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

// Test set_fuel_stream_for_phi
TEST_F(ThermoTransportTest, SetFuelStreamForPhi) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    
    // Fuel stream: pure CH4
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = standard_dry_air_composition();
    air.mdot = 10.0;
    
    // Test stoichiometric case (phi = 1.0)
    Stream fuel_stoich = set_fuel_stream_for_phi(1.0, fuel, air);
    EXPECT_GT(fuel_stoich.mdot, 0.0);
    
    // Mix and verify equivalence ratio
    Stream mixed = mix({fuel_stoich, air});
    std::vector<double> Y_fuel = mole_to_mass(fuel.X());
    std::vector<double> Y_air = mole_to_mass(air.X());
    std::vector<double> Y_mix = mole_to_mass(mixed.X());
    
    double phi_check = equivalence_ratio_mass(Y_mix, Y_fuel, Y_air);
    EXPECT_NEAR(phi_check, 1.0, 0.01);
    
    // Test lean case (phi = 0.5)
    Stream fuel_lean = set_fuel_stream_for_phi(0.5, fuel, air);
    EXPECT_LT(fuel_lean.mdot, fuel_stoich.mdot);
    
    Stream mixed_lean = mix({fuel_lean, air});
    std::vector<double> Y_mix_lean = mole_to_mass(mixed_lean.X());
    double phi_lean = equivalence_ratio_mass(Y_mix_lean, Y_fuel, Y_air);
    EXPECT_NEAR(phi_lean, 0.5, 0.01);
}

// =============================================================================
// Hydrocarbon Mixture Tests (Combustion, Mixing, Thermo)
// =============================================================================

// Test thermo properties for natural gas mixture
TEST_F(ThermoTransportTest, ThermoPropertiesNaturalGasMixture) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    
    // Natural gas composition: ~90% CH4, 5% C2H6, 2% C3H8, 2% N2, 1% CO2
    State ng;
    ng.T = 300.0;
    ng.P = 101325.0;
    ng.X = std::vector<double>(n, 0.0);
    ng.X[idx_CH4] = 0.90;
    ng.X[idx_C2H6] = 0.05;
    ng.X[idx_C3H8] = 0.02;
    ng.X[idx_N2] = 0.02;
    ng.X[idx_CO2] = 0.01;
    
    // Verify thermo properties are reasonable
    EXPECT_GT(ng.cp(), 30.0);  // J/(mol·K)
    EXPECT_LT(ng.cp(), 50.0);
    EXPECT_GT(ng.rho(), 0.5);  // kg/m³
    EXPECT_LT(ng.rho(), 1.5);
    EXPECT_GT(ng.mw(), 16.0);  // g/mol (CH4 is 16)
    EXPECT_LT(ng.mw(), 20.0);
    
    // Enthalpy should be negative (formation enthalpy of hydrocarbons)
    EXPECT_LT(h(ng.T, ng.X), 0.0);
    
    // Transport properties
    EXPECT_GT(ng.mu(), 1e-6);  // Pa·s
    EXPECT_LT(ng.mu(), 2e-5);
    EXPECT_GT(ng.Pr(), 0.6);
    EXPECT_LT(ng.Pr(), 1.0);
}

// Test mixing of natural gas with humid air
TEST_F(ThermoTransportTest, MixingNaturalGasWithHumidAir) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_AR = species_index_from_name("AR");
    
    // Natural gas fuel
    State ng;
    ng.T = 300.0;
    ng.P = 101325.0;
    ng.X = std::vector<double>(n, 0.0);
    ng.X[idx_CH4] = 0.90;
    ng.X[idx_C2H6] = 0.05;
    ng.X[idx_C3H8] = 0.02;
    ng.X[idx_N2] = 0.02;
    ng.X[idx_CO2] = 0.01;
    
    Stream fuel;
    fuel.state = ng;
    fuel.mdot = 1.0;  // 1 kg/s
    
    // Humid air
    std::vector<double> X_air = humid_air_composition(300.0, 101325.0, 0.5);
    State air_state;
    air_state.T = 300.0;
    air_state.P = 101325.0;
    air_state.X = X_air;
    Stream air;
    air.state = air_state;
    air.mdot = 20.0;  // 20 kg/s
    
    // Mix streams
    Stream mixed = mix({fuel, air});
    
    // Verify mixing conserves mass
    EXPECT_NEAR(mixed.mdot, fuel.mdot + air.mdot, 1e-10);
    
    // Verify all species are present
    EXPECT_GT(mixed.X()[idx_CH4], 0.0);
    EXPECT_GT(mixed.X()[idx_C2H6], 0.0);
    EXPECT_GT(mixed.X()[idx_C3H8], 0.0);
    EXPECT_GT(mixed.X()[idx_O2], 0.0);
    EXPECT_GT(mixed.X()[idx_N2], 0.0);
    EXPECT_GT(mixed.X()[idx_H2O], 0.0);
    EXPECT_GT(mixed.X()[idx_AR], 0.0);
    
    // Mole fractions should sum to 1
    double sum = 0.0;
    for (double x : mixed.X()) sum += x;
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

// Test complete combustion of natural gas mixture
TEST_F(ThermoTransportTest, CompleteCombustionNaturalGas) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_AR = species_index_from_name("AR");
    
    // Natural gas fuel
    State ng;
    ng.T = 300.0;
    ng.P = 101325.0;
    ng.X = std::vector<double>(n, 0.0);
    ng.X[idx_CH4] = 0.90;
    ng.X[idx_C2H6] = 0.05;
    ng.X[idx_C3H8] = 0.02;
    ng.X[idx_N2] = 0.02;
    ng.X[idx_CO2] = 0.01;
    
    Stream fuel;
    fuel.state = ng;
    fuel.mdot = 1.0;
    
    // Humid air
    std::vector<double> X_air = humid_air_composition(300.0, 101325.0, 0.5);
    State air_state;
    air_state.T = 300.0;
    air_state.P = 101325.0;
    air_state.X = X_air;
    Stream air;
    air.state = air_state;
    air.mdot = 20.0;  // Excess air for complete combustion
    
    // Mix and combust
    Stream mixed = mix({fuel, air});
    State burned = complete_combustion(mixed.state);
    
    // All hydrocarbons should be consumed
    EXPECT_LT(burned.X[idx_CH4], 1e-10);
    EXPECT_LT(burned.X[idx_C2H6], 1e-10);
    EXPECT_LT(burned.X[idx_C3H8], 1e-10);
    
    // CO2 and H2O should be produced
    EXPECT_GT(burned.X[idx_CO2], 0.05);
    EXPECT_GT(burned.X[idx_H2O], 0.10);
    
    // Temperature should increase significantly
    EXPECT_GT(burned.T, 1500.0);
    
    // Inerts should still be present
    EXPECT_GT(burned.X[idx_N2], 0.5);
    EXPECT_GT(burned.X[idx_AR], 0.001);
}

// Test equivalence ratio calculation with natural gas
TEST_F(ThermoTransportTest, EquivalenceRatioNaturalGas) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    
    // Natural gas fuel
    State ng;
    ng.T = 300.0;
    ng.P = 101325.0;
    ng.X = std::vector<double>(n, 0.0);
    ng.X[idx_CH4] = 0.90;
    ng.X[idx_C2H6] = 0.05;
    ng.X[idx_C3H8] = 0.02;
    ng.X[idx_N2] = 0.02;
    ng.X[idx_CO2] = 0.01;
    
    Stream fuel;
    fuel.state = ng;
    fuel.mdot = 1.0;
    
    // Humid air
    std::vector<double> X_air = humid_air_composition(300.0, 101325.0, 0.5);
    State air_state;
    air_state.T = 300.0;
    air_state.P = 101325.0;
    air_state.X = X_air;
    Stream air;
    air.state = air_state;
    air.mdot = 10.0;
    
    // Test set_fuel_stream_for_phi with natural gas
    Stream fuel_stoich = set_fuel_stream_for_phi(1.0, fuel, air);
    EXPECT_GT(fuel_stoich.mdot, 0.0);
    
    // Mix and verify equivalence ratio
    Stream mixed = mix({fuel_stoich, air});
    std::vector<double> Y_fuel = mole_to_mass(fuel.X());
    std::vector<double> Y_air = mole_to_mass(air.X());
    std::vector<double> Y_mix = mole_to_mass(mixed.X());
    
    double phi_check = equivalence_ratio_mass(Y_mix, Y_fuel, Y_air);
    EXPECT_NEAR(phi_check, 1.0, 0.02);
    
    // Test rich case
    Stream fuel_rich = set_fuel_stream_for_phi(1.2, fuel, air);
    Stream mixed_rich = mix({fuel_rich, air});
    std::vector<double> Y_mix_rich = mole_to_mass(mixed_rich.X());
    double phi_rich = equivalence_ratio_mass(Y_mix_rich, Y_fuel, Y_air);
    EXPECT_NEAR(phi_rich, 1.2, 0.02);
}

// Test LPG (propane/butane) combustion
TEST_F(ThermoTransportTest, CompleteCombustionLPG) {
    const std::size_t n = species_names.size();
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_IC4H10 = species_index_from_name("IC4H10");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_AR = species_index_from_name("AR");
    
    // LPG: ~60% propane, 40% butane
    State lpg;
    lpg.T = 300.0;
    lpg.P = 101325.0;
    lpg.X = std::vector<double>(n, 0.0);
    lpg.X[idx_C3H8] = 0.60;
    lpg.X[idx_IC4H10] = 0.40;
    
    Stream fuel;
    fuel.state = lpg;
    fuel.mdot = 1.0;
    
    // Dry air
    std::vector<double> X_air = humid_air_composition(300.0, 101325.0, 0.0);
    State air_state;
    air_state.T = 300.0;
    air_state.P = 101325.0;
    air_state.X = X_air;
    Stream air;
    air.state = air_state;
    air.mdot = 25.0;  // Excess air
    
    // Mix and combust
    Stream mixed = mix({fuel, air});
    State burned = complete_combustion(mixed.state);
    
    // All hydrocarbons consumed
    EXPECT_LT(burned.X[idx_C3H8], 1e-10);
    EXPECT_LT(burned.X[idx_IC4H10], 1e-10);
    
    // Products formed
    EXPECT_GT(burned.X[idx_CO2], 0.05);
    EXPECT_GT(burned.X[idx_H2O], 0.08);
    
    // Temperature increase
    EXPECT_GT(burned.T, 1500.0);
    
    // Inerts present
    EXPECT_GT(burned.X[idx_N2], 0.5);
    EXPECT_GT(burned.X[idx_AR], 0.001);
}

// Test rich combustion with natural gas produces unburned hydrocarbons
TEST_F(ThermoTransportTest, RichCombustionNaturalGasUnburnedHC) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_H2 = species_index_from_name("H2");
    
    // Natural gas fuel
    State ng;
    ng.T = 300.0;
    ng.P = 101325.0;
    ng.X = std::vector<double>(n, 0.0);
    ng.X[idx_CH4] = 0.90;
    ng.X[idx_C2H6] = 0.05;
    ng.X[idx_C3H8] = 0.02;
    ng.X[idx_N2] = 0.02;
    ng.X[idx_CO2] = 0.01;
    
    Stream fuel;
    fuel.state = ng;
    fuel.mdot = 1.0;
    
    // Humid air
    std::vector<double> X_air = humid_air_composition(300.0, 101325.0, 0.5);
    State air_state;
    air_state.T = 300.0;
    air_state.P = 101325.0;
    air_state.X = X_air;
    Stream air;
    air.state = air_state;
    air.mdot = 10.0;
    
    // Rich mixture (phi = 1.3)
    Stream fuel_rich = set_fuel_stream_for_phi(1.3, fuel, air);
    Stream mixed = mix({fuel_rich, air});
    State burned = complete_combustion(mixed.state);
    
    // Some unburned hydrocarbons should remain
    double total_hc = burned.X[idx_CH4] + burned.X[idx_C2H6] + burned.X[idx_C3H8];
    EXPECT_GT(total_hc, 0.005);
    
    // Apply reforming equilibrium
    State eq = reforming_equilibrium_adiabatic(burned);
    
    // Hydrocarbons should be reformed
    double total_hc_eq = eq.X[idx_CH4] + eq.X[idx_C2H6] + eq.X[idx_C3H8];
    EXPECT_LT(total_hc_eq, total_hc * 0.5);
    
    // CO and H2 should be produced
    EXPECT_GT(eq.X[idx_CO], 0.01);
    EXPECT_GT(eq.X[idx_H2], 0.01);
    
    // Temperature should drop (endothermic reforming)
    EXPECT_LT(eq.T, burned.T);
}

// =============================================================================
// SMR+WGS Equilibrium Tests
// =============================================================================

// Test SMR+WGS equilibrium at high temperature with CH4
TEST_F(ThermoTransportTest, SmrWgsEquilibriumIsothermal) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2 = species_index_from_name("H2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Create a mixture with CH4, H2O, CO2, and N2 (typical rich combustion products)
    State in;
    in.T = 2000.0;  // High temperature favors SMR
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.02;   // 2% CH4
    in.X[idx_H2O] = 0.20;   // 20% H2O
    in.X[idx_CO2] = 0.10;   // 10% CO2
    in.X[idx_N2] = 0.68;    // 68% N2
    
    State out = smr_wgs_equilibrium(in);
    
    // Temperature should be unchanged (isothermal)
    EXPECT_NEAR(out.T, in.T, 1e-6);
    
    // At high T, SMR should convert CH4 to CO + H2
    // CH4 should decrease
    EXPECT_LT(out.X[idx_CH4], in.X[idx_CH4]);
    
    // CO and H2 should increase
    EXPECT_GT(out.X[idx_CO], in.X[idx_CO]);
    EXPECT_GT(out.X[idx_H2], in.X[idx_H2]);
    
    // Mole fractions should sum to 1
    double sum = 0.0;
    for (double x : out.X) sum += x;
    EXPECT_NEAR(sum, 1.0, 1e-10);
    
    // Element balance: use N2 as reference (inert)
    // SMR changes total moles, so ratios to N2 should be conserved
    double n2_in = in.X[idx_N2];
    double n2_out = out.X[idx_N2];
    
    double C_per_N2_in = (in.X[idx_CH4] + in.X[idx_CO] + in.X[idx_CO2]) / n2_in;
    double C_per_N2_out = (out.X[idx_CH4] + out.X[idx_CO] + out.X[idx_CO2]) / n2_out;
    EXPECT_NEAR(C_per_N2_in, C_per_N2_out, 1e-6);
    
    double H_per_N2_in = (4.0 * in.X[idx_CH4] + 2.0 * in.X[idx_H2O] + 2.0 * in.X[idx_H2]) / n2_in;
    double H_per_N2_out = (4.0 * out.X[idx_CH4] + 2.0 * out.X[idx_H2O] + 2.0 * out.X[idx_H2]) / n2_out;
    EXPECT_NEAR(H_per_N2_in, H_per_N2_out, 1e-6);
}

// Test SMR+WGS adiabatic equilibrium
TEST_F(ThermoTransportTest, SmrWgsEquilibriumAdiabatic) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Create a mixture with CH4 (typical rich combustion products)
    State in;
    in.T = 2200.0;  // Start at high temperature
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.02;
    in.X[idx_H2O] = 0.20;
    in.X[idx_CO2] = 0.10;
    in.X[idx_N2] = 0.68;
    
    State out = smr_wgs_equilibrium_adiabatic(in);
    
    // SMR is endothermic, so temperature should decrease
    EXPECT_LT(out.T, in.T);
    
    // Temperature should still be reasonable (not too low)
    EXPECT_GT(out.T, 1500.0);
    
    // Enthalpy per mole of N2 should be conserved (N2 is inert)
    // H_total / n_N2 = h(T,X) / X_N2 should be constant
    // Note: tolerance is ~5% due to nested Newton solver convergence
    double H_per_N2_in = h(in.T, in.X) / in.X[idx_N2];
    double H_per_N2_out = h(out.T, out.X) / out.X[idx_N2];
    EXPECT_NEAR(H_per_N2_in, H_per_N2_out, std::abs(H_per_N2_in) * 0.05);
}

// Test SMR+WGS falls back to WGS when no CH4
TEST_F(ThermoTransportTest, SmrWgsFallbackToWgs) {
    const std::size_t n = species_names.size();
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Create a mixture without CH4 (stoichiometric combustion products)
    State in;
    in.T = 2200.0;
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_H2O] = 0.20;
    in.X[idx_CO2] = 0.10;
    in.X[idx_N2] = 0.70;
    
    State out_smr_wgs = smr_wgs_equilibrium_adiabatic(in);
    State out_wgs = wgs_equilibrium_adiabatic(in);
    
    // Results should be identical when no CH4
    EXPECT_NEAR(out_smr_wgs.T, out_wgs.T, 1e-6);
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(out_smr_wgs.X[i], out_wgs.X[i], 1e-10);
    }
}

// Test SMR+WGS with complete combustion products at rich conditions
TEST_F(ThermoTransportTest, SmrWgsRichCombustion) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_H2 = species_index_from_name("H2");
    
    // Create fuel and air streams
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = standard_dry_air_composition();
    air.mdot = 10.0;
    
    // Rich mixture (phi = 1.2)
    Stream fuel_rich = set_fuel_stream_for_phi(1.2, fuel, air);
    Stream mixed = mix({fuel_rich, air});
    
    // Complete combustion
    State burned = complete_combustion(mixed.state);
    
    // Verify there's unburned CH4
    EXPECT_GT(burned.X[idx_CH4], 0.01);
    
    // Apply SMR+WGS equilibrium
    State eq = smr_wgs_equilibrium_adiabatic(burned);
    
    // CH4 should be mostly reformed
    EXPECT_LT(eq.X[idx_CH4], burned.X[idx_CH4] * 0.5);
    
    // CO and H2 should be present
    EXPECT_GT(eq.X[idx_CO], 0.01);
    EXPECT_GT(eq.X[idx_H2], 0.01);
    
    // Temperature should drop (SMR is endothermic)
    EXPECT_LT(eq.T, burned.T);
    
    // Enthalpy per mole of N2 should be conserved
    // Note: tolerance is ~10% due to nested Newton solver convergence
    const std::size_t idx_N2 = species_index_from_name("N2");
    double H_per_N2_burned = h(burned.T, burned.X) / burned.X[idx_N2];
    double H_per_N2_eq = h(eq.T, eq.X) / eq.X[idx_N2];
    EXPECT_NEAR(H_per_N2_burned, H_per_N2_eq, std::abs(H_per_N2_burned) * 0.10);
}

// Test general reforming equilibrium with multiple hydrocarbons
TEST_F(ThermoTransportTest, ReformingEquilibriumMultipleHydrocarbons) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2 = species_index_from_name("H2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Create a mixture with multiple hydrocarbons (simulating natural gas combustion)
    State in;
    in.T = 2000.0;  // High temperature
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.015;   // 1.5% CH4
    in.X[idx_C2H6] = 0.003;  // 0.3% C2H6
    in.X[idx_C3H8] = 0.002;  // 0.2% C3H8
    in.X[idx_H2O] = 0.20;    // 20% H2O
    in.X[idx_CO2] = 0.08;    // 8% CO2
    in.X[idx_N2] = 0.70;     // 70% N2
    
    State out = reforming_equilibrium(in);
    
    // All hydrocarbons should be mostly reformed at high T
    EXPECT_LT(out.X[idx_CH4], in.X[idx_CH4] * 0.5);
    EXPECT_LT(out.X[idx_C2H6], in.X[idx_C2H6] * 0.5);
    EXPECT_LT(out.X[idx_C3H8], in.X[idx_C3H8] * 0.5);
    
    // CO and H2 should be produced
    EXPECT_GT(out.X[idx_CO], 0.01);
    EXPECT_GT(out.X[idx_H2], 0.01);
    
    // Element balance using N2 as reference
    double n2_in = in.X[idx_N2];
    double n2_out = out.X[idx_N2];
    
    // C atoms: CH4 + 2*C2H6 + 3*C3H8 + CO + CO2
    double C_per_N2_in = (in.X[idx_CH4] + 2*in.X[idx_C2H6] + 3*in.X[idx_C3H8] 
                        + in.X[idx_CO] + in.X[idx_CO2]) / n2_in;
    double C_per_N2_out = (out.X[idx_CH4] + 2*out.X[idx_C2H6] + 3*out.X[idx_C3H8]
                         + out.X[idx_CO] + out.X[idx_CO2]) / n2_out;
    EXPECT_NEAR(C_per_N2_in, C_per_N2_out, 1e-6);
}

// Test general reforming adiabatic equilibrium
TEST_F(ThermoTransportTest, ReformingEquilibriumAdiabatic) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Create a mixture with multiple hydrocarbons
    State in;
    in.T = 2200.0;
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.015;
    in.X[idx_C2H6] = 0.005;
    in.X[idx_H2O] = 0.20;
    in.X[idx_CO2] = 0.08;
    in.X[idx_N2] = 0.70;
    
    State out = reforming_equilibrium_adiabatic(in);
    
    // Reforming is endothermic, so temperature should decrease
    EXPECT_LT(out.T, in.T);
    
    // Temperature should be reasonable
    EXPECT_GT(out.T, 1500.0);
}

// Test reforming with natural gas composition (CH4 + C2H6 + C3H8 + inerts)
TEST_F(ThermoTransportTest, ReformingNaturalGasCombustion) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2 = species_index_from_name("H2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_AR = species_index_from_name("AR");
    
    // Simulate rich combustion products from natural gas
    // Natural gas: ~90% CH4, ~5% C2H6, ~2% C3H8, ~2% N2, ~1% CO2
    // After rich combustion, some unburned hydrocarbons remain
    State in;
    in.T = 2100.0;
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.012;   // Unburned CH4
    in.X[idx_C2H6] = 0.002;  // Unburned C2H6
    in.X[idx_C3H8] = 0.001;  // Unburned C3H8
    in.X[idx_H2O] = 0.18;    // Combustion product
    in.X[idx_CO2] = 0.09;    // Combustion product + fuel inert
    in.X[idx_N2] = 0.705;    // Air N2 + fuel N2
    in.X[idx_AR] = 0.01;     // Air Ar
    
    State out = reforming_equilibrium_adiabatic(in);
    
    // All hydrocarbons should be significantly reformed
    EXPECT_LT(out.X[idx_CH4], in.X[idx_CH4] * 0.3);
    EXPECT_LT(out.X[idx_C2H6], in.X[idx_C2H6] * 0.3);
    EXPECT_LT(out.X[idx_C3H8], in.X[idx_C3H8] * 0.3);
    
    // CO and H2 should be produced
    EXPECT_GT(out.X[idx_CO], 0.01);
    EXPECT_GT(out.X[idx_H2], 0.01);
    
    // Inerts should be unchanged (relative to total moles via N2 reference)
    double ar_per_n2_in = in.X[idx_AR] / in.X[idx_N2];
    double ar_per_n2_out = out.X[idx_AR] / out.X[idx_N2];
    EXPECT_NEAR(ar_per_n2_in, ar_per_n2_out, 1e-10);
    
    // Temperature should drop (endothermic reforming)
    EXPECT_LT(out.T, in.T);
}

// Test reforming with heavier hydrocarbons (C4+)
TEST_F(ThermoTransportTest, ReformingHeavyHydrocarbons) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_IC4H10 = species_index_from_name("IC4H10");
    const std::size_t idx_NC5H12 = species_index_from_name("NC5H12");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2 = species_index_from_name("H2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_AR = species_index_from_name("AR");
    
    // Mixture with heavier hydrocarbons (simulating LPG or gasoline-like fuel)
    State in;
    in.T = 2000.0;
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.005;    // Small amount of CH4
    in.X[idx_C3H8] = 0.008;   // Propane
    in.X[idx_IC4H10] = 0.004; // Isobutane
    in.X[idx_NC5H12] = 0.002; // n-Pentane
    in.X[idx_H2O] = 0.20;
    in.X[idx_CO2] = 0.08;
    in.X[idx_N2] = 0.69;
    in.X[idx_AR] = 0.011;
    
    State out = reforming_equilibrium(in);
    
    // All hydrocarbons should be reformed at high T
    EXPECT_LT(out.X[idx_CH4], in.X[idx_CH4] * 0.5);
    EXPECT_LT(out.X[idx_C3H8], in.X[idx_C3H8] * 0.5);
    EXPECT_LT(out.X[idx_IC4H10], in.X[idx_IC4H10] * 0.5);
    EXPECT_LT(out.X[idx_NC5H12], in.X[idx_NC5H12] * 0.5);
    
    // CO and H2 should be produced
    EXPECT_GT(out.X[idx_CO], 0.01);
    EXPECT_GT(out.X[idx_H2], 0.01);
    
    // Element balance: C atoms per N2 should be conserved
    double n2_in = in.X[idx_N2];
    double n2_out = out.X[idx_N2];
    
    // C atoms: 1*CH4 + 3*C3H8 + 4*iC4H10 + 5*nC5H12 + CO + CO2
    double C_per_N2_in = (in.X[idx_CH4] + 3*in.X[idx_C3H8] + 4*in.X[idx_IC4H10] 
                        + 5*in.X[idx_NC5H12] + in.X[idx_CO] + in.X[idx_CO2]) / n2_in;
    double C_per_N2_out = (out.X[idx_CH4] + 3*out.X[idx_C3H8] + 4*out.X[idx_IC4H10]
                         + 5*out.X[idx_NC5H12] + out.X[idx_CO] + out.X[idx_CO2]) / n2_out;
    EXPECT_NEAR(C_per_N2_in, C_per_N2_out, 1e-6);
    
    // H atoms per N2 should be conserved
    // H atoms: 4*CH4 + 8*C3H8 + 10*iC4H10 + 12*nC5H12 + 2*H2O + 2*H2
    double H_per_N2_in = (4*in.X[idx_CH4] + 8*in.X[idx_C3H8] + 10*in.X[idx_IC4H10]
                        + 12*in.X[idx_NC5H12] + 2*in.X[idx_H2O] + 2*in.X[idx_H2]) / n2_in;
    double H_per_N2_out = (4*out.X[idx_CH4] + 8*out.X[idx_C3H8] + 10*out.X[idx_IC4H10]
                         + 12*out.X[idx_NC5H12] + 2*out.X[idx_H2O] + 2*out.X[idx_H2]) / n2_out;
    EXPECT_NEAR(H_per_N2_in, H_per_N2_out, 1e-6);
    
    // Ar should be unchanged relative to N2
    double ar_per_n2_in = in.X[idx_AR] / n2_in;
    double ar_per_n2_out = out.X[idx_AR] / n2_out;
    EXPECT_NEAR(ar_per_n2_in, ar_per_n2_out, 1e-10);
}

// Test reforming with only C2+ hydrocarbons (no CH4)
TEST_F(ThermoTransportTest, ReformingC2PlusOnly) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_C3H8 = species_index_from_name("C3H8");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2 = species_index_from_name("H2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Mixture with only C2+ hydrocarbons (no CH4)
    State in;
    in.T = 2000.0;
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.0;      // No CH4!
    in.X[idx_C2H6] = 0.010;   // Ethane only
    in.X[idx_C3H8] = 0.005;   // Propane
    in.X[idx_H2O] = 0.20;
    in.X[idx_CO2] = 0.085;
    in.X[idx_N2] = 0.70;
    
    State out = reforming_equilibrium(in);
    
    // C2H6 and C3H8 should be reformed even without CH4
    EXPECT_LT(out.X[idx_C2H6], in.X[idx_C2H6] * 0.5);
    EXPECT_LT(out.X[idx_C3H8], in.X[idx_C3H8] * 0.5);
    
    // CO and H2 should be produced
    EXPECT_GT(out.X[idx_CO], 0.01);
    EXPECT_GT(out.X[idx_H2], 0.01);
    
    // Element balance: C atoms per N2
    double n2_in = in.X[idx_N2];
    double n2_out = out.X[idx_N2];
    double C_per_N2_in = (in.X[idx_CH4] + 2*in.X[idx_C2H6] + 3*in.X[idx_C3H8]
                        + in.X[idx_CO] + in.X[idx_CO2]) / n2_in;
    double C_per_N2_out = (out.X[idx_CH4] + 2*out.X[idx_C2H6] + 3*out.X[idx_C3H8]
                         + out.X[idx_CO] + out.X[idx_CO2]) / n2_out;
    EXPECT_NEAR(C_per_N2_in, C_per_N2_out, 1e-6);
}

// Test reforming with pre-existing CO and H2
TEST_F(ThermoTransportTest, ReformingWithExistingCOH2) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2 = species_index_from_name("H2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_AR = species_index_from_name("AR");
    
    // Mixture with pre-existing CO and H2 (partial reforming already occurred)
    State in;
    in.T = 2000.0;
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.010;
    in.X[idx_C2H6] = 0.003;
    in.X[idx_H2O] = 0.18;
    in.X[idx_CO] = 0.02;     // Pre-existing CO
    in.X[idx_CO2] = 0.07;
    in.X[idx_H2] = 0.01;     // Pre-existing H2
    in.X[idx_N2] = 0.70;
    in.X[idx_AR] = 0.007;
    
    State out = reforming_equilibrium(in);
    
    // Hydrocarbons should still be reformed
    EXPECT_LT(out.X[idx_CH4], in.X[idx_CH4] * 0.5);
    EXPECT_LT(out.X[idx_C2H6], in.X[idx_C2H6] * 0.5);
    
    // CO and H2 should increase further
    EXPECT_GT(out.X[idx_CO], in.X[idx_CO]);
    EXPECT_GT(out.X[idx_H2], in.X[idx_H2]);
    
    // Element balance
    double n2_in = in.X[idx_N2];
    double n2_out = out.X[idx_N2];
    
    double C_per_N2_in = (in.X[idx_CH4] + 2*in.X[idx_C2H6] + in.X[idx_CO] + in.X[idx_CO2]) / n2_in;
    double C_per_N2_out = (out.X[idx_CH4] + 2*out.X[idx_C2H6] + out.X[idx_CO] + out.X[idx_CO2]) / n2_out;
    EXPECT_NEAR(C_per_N2_in, C_per_N2_out, 1e-6);
    
    // Ar unchanged
    EXPECT_NEAR(in.X[idx_AR] / n2_in, out.X[idx_AR] / n2_out, 1e-10);
}

// Test element conservation in SMR+WGS
// Note: SMR changes total moles (Δn=2), so we need to account for this
// when checking element balance. We use the ratio of elements to N2
// (which is inert) to verify conservation.
TEST_F(ThermoTransportTest, SmrWgsElementConservation) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2 = species_index_from_name("H2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Test at multiple temperatures
    for (double T : {1500.0, 2000.0, 2500.0}) {
        State in;
        in.T = T;
        in.P = 101325.0;
        in.X = std::vector<double>(n, 0.0);
        in.X[idx_CH4] = 0.03;
        in.X[idx_H2O] = 0.18;
        in.X[idx_CO2] = 0.08;
        in.X[idx_CO] = 0.01;
        in.X[idx_H2] = 0.02;
        in.X[idx_N2] = 0.68;
        
        State out = smr_wgs_equilibrium(in);
        
        // Use N2 as reference (inert, doesn't participate in reactions)
        // Element ratios relative to N2 should be conserved
        double n2_in = in.X[idx_N2];
        double n2_out = out.X[idx_N2];
        
        // C atoms per N2: (CH4 + CO + CO2) / N2
        double C_per_N2_in = (in.X[idx_CH4] + in.X[idx_CO] + in.X[idx_CO2]) / n2_in;
        double C_per_N2_out = (out.X[idx_CH4] + out.X[idx_CO] + out.X[idx_CO2]) / n2_out;
        EXPECT_NEAR(C_per_N2_in, C_per_N2_out, 1e-6) << "C/N2 not conserved at T=" << T;
        
        // H atoms per N2: (4*CH4 + 2*H2O + 2*H2) / N2
        double H_per_N2_in = (4.0 * in.X[idx_CH4] + 2.0 * in.X[idx_H2O] + 2.0 * in.X[idx_H2]) / n2_in;
        double H_per_N2_out = (4.0 * out.X[idx_CH4] + 2.0 * out.X[idx_H2O] + 2.0 * out.X[idx_H2]) / n2_out;
        EXPECT_NEAR(H_per_N2_in, H_per_N2_out, 1e-6) << "H/N2 not conserved at T=" << T;
        
        // O atoms per N2: (H2O + CO + 2*CO2) / N2
        double O_per_N2_in = (in.X[idx_H2O] + in.X[idx_CO] + 2.0 * in.X[idx_CO2]) / n2_in;
        double O_per_N2_out = (out.X[idx_H2O] + out.X[idx_CO] + 2.0 * out.X[idx_CO2]) / n2_out;
        EXPECT_NEAR(O_per_N2_in, O_per_N2_out, 1e-6) << "O/N2 not conserved at T=" << T;
    }
}

// =============================================================================
// Combustion + Equilibrium Convenience Function Tests
// =============================================================================

// Test combustion_equilibrium with unburned fuel+air mixture
TEST_F(ThermoTransportTest, CombustionEquilibriumFromUnburned) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_H2 = species_index_from_name("H2");
    
    // Rich CH4 + air mixture (unburned)
    State in;
    in.T = 300.0;
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.11;   // Rich
    in.X[idx_O2] = 0.19;
    in.X[idx_N2] = 0.70;
    
    // One-step: combustion_equilibrium
    State result = combustion_equilibrium(in);
    
    // Two-step equivalent: complete_combustion + reforming_equilibrium_adiabatic
    State burned = complete_combustion(in);
    State eq = reforming_equilibrium_adiabatic(burned);
    
    // Results should be identical
    EXPECT_NEAR(result.T, eq.T, 0.1);
    EXPECT_NEAR(result.X[idx_CO], eq.X[idx_CO], 1e-6);
    EXPECT_NEAR(result.X[idx_H2], eq.X[idx_H2], 1e-6);
    EXPECT_NEAR(result.X[idx_CO2], eq.X[idx_CO2], 1e-6);
    EXPECT_NEAR(result.X[idx_H2O], eq.X[idx_H2O], 1e-6);
    
    // Temperature should be high (combustion occurred)
    EXPECT_GT(result.T, 2000.0);
    
    // CH4 should be fully consumed
    EXPECT_LT(result.X[idx_CH4], 1e-6);
    
    // CO and H2 should be present (rich combustion + reforming)
    EXPECT_GT(result.X[idx_CO], 0.01);
    EXPECT_GT(result.X[idx_H2], 0.01);
}

// Test combustion_equilibrium with stoichiometric mixture
TEST_F(ThermoTransportTest, CombustionEquilibriumStoichiometric) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    const std::size_t idx_H2O = species_index_from_name("H2O");
    const std::size_t idx_CO = species_index_from_name("CO");
    const std::size_t idx_H2 = species_index_from_name("H2");
    
    // Stoichiometric CH4 + air mixture
    State in;
    in.T = 300.0;
    in.P = 101325.0;
    in.X = std::vector<double>(n, 0.0);
    in.X[idx_CH4] = 0.095;  // Stoichiometric
    in.X[idx_O2] = 0.19;
    in.X[idx_N2] = 0.715;
    
    State result = combustion_equilibrium(in);
    
    // Temperature should be high
    EXPECT_GT(result.T, 2200.0);
    
    // CH4 and O2 should be fully consumed
    EXPECT_LT(result.X[idx_CH4], 1e-6);
    EXPECT_LT(result.X[idx_O2], 1e-6);
    
    // Products should be mainly CO2 and H2O
    EXPECT_GT(result.X[idx_CO2], 0.05);
    EXPECT_GT(result.X[idx_H2O], 0.10);
    
    // Little CO/H2 at stoichiometric (WGS equilibrium)
    EXPECT_LT(result.X[idx_CO], 0.01);
    EXPECT_LT(result.X[idx_H2], 0.01);
}

// -------------------------------------------------------------
// Inverse solvers for fuel stream tests
// -------------------------------------------------------------

// Test set_fuel_stream_for_Tad
TEST_F(ThermoTransportTest, SetFuelStreamForTad) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Fuel stream: pure CH4 at 300 K
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target Tad = 1500 K (lean combustion)
    double T_ad_target = 1500.0;
    Stream fuel_result = set_fuel_stream_for_Tad(T_ad_target, fuel, air);
    
    // Verify fuel mdot is positive and reasonable
    EXPECT_GT(fuel_result.mdot, 0.0);
    EXPECT_LT(fuel_result.mdot, 1.0);  // Should be lean
    
    // Verify the result by mixing and combusting
    Stream mixed = mix({fuel_result, air});
    State burned = complete_combustion(mixed.state);
    
    // Check that achieved Tad is close to target
    EXPECT_NEAR(burned.T, T_ad_target, 2.0);  // Within 2 K
    
    // Verify O2 is present in products (lean combustion)
    EXPECT_GT(burned.X[idx_O2], 0.0);
}

// Test set_fuel_stream_for_O2
TEST_F(ThermoTransportTest, SetFuelStreamForO2) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Fuel stream: pure CH4 at 300 K
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target O2 = 10% in burned products (lean combustion)
    double X_O2_target = 0.10;
    Stream fuel_result = set_fuel_stream_for_O2(X_O2_target, fuel, air);
    
    // Verify fuel mdot is positive and reasonable
    EXPECT_GT(fuel_result.mdot, 0.0);
    EXPECT_LT(fuel_result.mdot, 1.0);  // Should be lean
    
    // Verify the result by mixing and combusting
    Stream mixed = mix({fuel_result, air});
    std::vector<double> X_burned = complete_combustion_to_CO2_H2O(mixed.state.X);
    
    // Check that achieved O2 is close to target
    EXPECT_NEAR(X_burned[idx_O2], X_O2_target, 1e-5);
}

// Test set_fuel_stream_for_CO2
TEST_F(ThermoTransportTest, SetFuelStreamForCO2) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    
    // Fuel stream: pure CH4 at 300 K
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target CO2 = 5% in burned products (lean combustion)
    double X_CO2_target = 0.05;
    Stream fuel_result = set_fuel_stream_for_CO2(X_CO2_target, fuel, air);
    
    // Verify fuel mdot is positive and reasonable
    EXPECT_GT(fuel_result.mdot, 0.0);
    EXPECT_LT(fuel_result.mdot, 1.0);  // Should be lean
    
    // Verify the result by mixing and combusting
    Stream mixed = mix({fuel_result, air});
    std::vector<double> X_burned = complete_combustion_to_CO2_H2O(mixed.state.X);
    
    // Check that achieved CO2 is close to target
    EXPECT_NEAR(X_burned[idx_CO2], X_CO2_target, 1e-5);
}

// Test set_fuel_stream_for_Tad with different fuel temperature
TEST_F(ThermoTransportTest, SetFuelStreamForTadHotFuel) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Fuel stream: pure CH4 at 500 K (preheated)
    Stream fuel;
    fuel.state.T = 500.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target Tad = 1500 K
    double T_ad_target = 1500.0;
    Stream fuel_result = set_fuel_stream_for_Tad(T_ad_target, fuel, air);
    
    // Verify the result
    Stream mixed = mix({fuel_result, air});
    State burned = complete_combustion(mixed.state);
    
    EXPECT_NEAR(burned.T, T_ad_target, 2.0);
}

// Test set_fuel_stream_for_O2 with multi-component fuel
TEST_F(ThermoTransportTest, SetFuelStreamForO2MultiComponentFuel) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_C2H6 = species_index_from_name("C2H6");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Fuel stream: 90% CH4 + 10% C2H6 at 300 K
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 0.90;
    fuel.state.X[idx_C2H6] = 0.10;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target O2 = 8% in burned products
    double X_O2_target = 0.08;
    Stream fuel_result = set_fuel_stream_for_O2(X_O2_target, fuel, air);
    
    // Verify the result
    Stream mixed = mix({fuel_result, air});
    std::vector<double> X_burned = complete_combustion_to_CO2_H2O(mixed.state.X);
    
    EXPECT_NEAR(X_burned[idx_O2], X_O2_target, 1e-5);
}

// Test set_oxidizer_stream_for_Tad
TEST_F(ThermoTransportTest, SetOxidizerStreamForTad) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Fuel stream: pure CH4 at 300 K, 0.5 kg/s
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    fuel.mdot = 0.5;
    
    // Oxidizer stream: dry air at 300 K (mdot to be solved)
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    
    // Target Tad = 1500 K (lean combustion)
    double T_ad_target = 1500.0;
    Stream air_result = set_oxidizer_stream_for_Tad(T_ad_target, fuel, air);
    
    // Verify air mdot is positive and reasonable
    EXPECT_GT(air_result.mdot, 0.0);
    
    // Verify the result by mixing and combusting
    Stream mixed = mix({fuel, air_result});
    State burned = complete_combustion(mixed.state);
    
    // Check that achieved Tad is close to target
    EXPECT_NEAR(burned.T, T_ad_target, 2.0);
    
    // Verify O2 is present in products (lean combustion)
    EXPECT_GT(burned.X[idx_O2], 0.0);
}

// Test set_fuel_stream_for_O2_dry
TEST_F(ThermoTransportTest, SetFuelStreamForO2Dry) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Fuel stream: pure CH4 at 300 K
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target dry O2 = 12% in burned products
    double X_O2_dry_target = 0.12;
    Stream fuel_result = set_fuel_stream_for_O2_dry(X_O2_dry_target, fuel, air);
    
    // Verify fuel mdot is positive
    EXPECT_GT(fuel_result.mdot, 0.0);
    
    // Verify the result
    Stream mixed = mix({fuel_result, air});
    std::vector<double> X_burned = complete_combustion_to_CO2_H2O(mixed.state.X);
    std::vector<double> X_dry = convert_to_dry_fractions(X_burned);
    
    EXPECT_NEAR(X_dry[idx_O2], X_O2_dry_target, 1e-5);
}

// Test set_fuel_stream_for_CO2_dry
TEST_F(ThermoTransportTest, SetFuelStreamForCO2Dry) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    
    // Fuel stream: pure CH4 at 300 K
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target dry CO2 = 6% in burned products
    double X_CO2_dry_target = 0.06;
    Stream fuel_result = set_fuel_stream_for_CO2_dry(X_CO2_dry_target, fuel, air);
    
    // Verify fuel mdot is positive
    EXPECT_GT(fuel_result.mdot, 0.0);
    
    // Verify the result
    Stream mixed = mix({fuel_result, air});
    std::vector<double> X_burned = complete_combustion_to_CO2_H2O(mixed.state.X);
    std::vector<double> X_dry = convert_to_dry_fractions(X_burned);
    
    EXPECT_NEAR(X_dry[idx_CO2], X_CO2_dry_target, 1e-5);
}

// -------------------------------------------------------------
// Round-trip tests: Create lean mixture at known phi, combust, 
// extract properties, then verify inverse solvers recover the same mdot
// -------------------------------------------------------------

TEST_F(ThermoTransportTest, InverseSolversRoundTripFromPhi) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    
    // Fuel stream: pure CH4 at 350 K
    Stream fuel;
    fuel.state.T = 350.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Test at phi = 0.7 (lean)
    double phi = 0.7;
    Stream fuel_phi = set_fuel_stream_for_phi(phi, fuel, air);
    double mdot_fuel_original = fuel_phi.mdot;
    
    // Mix and combust
    Stream mixed = mix({fuel_phi, air});
    State burned = complete_combustion(mixed.state);
    std::vector<double> X_burned = burned.X;
    std::vector<double> X_dry = convert_to_dry_fractions(X_burned);
    
    // Extract properties from burned products
    double Tad = burned.T;
    double X_O2 = X_burned[idx_O2];
    double X_CO2 = X_burned[idx_CO2];
    double X_O2_dry = X_dry[idx_O2];
    double X_CO2_dry = X_dry[idx_CO2];
    
    // Verify we have lean combustion products
    EXPECT_GT(X_O2, 0.0);
    EXPECT_GT(X_CO2, 0.0);
    EXPECT_GT(Tad, air.T());
    
    // Now test each inverse solver recovers the original fuel mdot
    const double tol_T = 1.0;
    const double tol_X = 1e-6;
    const double tol_mdot_rel = 0.01;  // 1% relative tolerance on mdot
    
    // Test set_fuel_stream_for_Tad
    Stream fuel_from_Tad = set_fuel_stream_for_Tad(Tad, fuel, air, tol_T);
    EXPECT_NEAR(fuel_from_Tad.mdot, mdot_fuel_original, mdot_fuel_original * tol_mdot_rel);
    
    // Test set_fuel_stream_for_O2
    Stream fuel_from_O2 = set_fuel_stream_for_O2(X_O2, fuel, air, tol_X);
    EXPECT_NEAR(fuel_from_O2.mdot, mdot_fuel_original, mdot_fuel_original * tol_mdot_rel);
    
    // Test set_fuel_stream_for_CO2
    Stream fuel_from_CO2 = set_fuel_stream_for_CO2(X_CO2, fuel, air, tol_X);
    EXPECT_NEAR(fuel_from_CO2.mdot, mdot_fuel_original, mdot_fuel_original * tol_mdot_rel);
    
    // Test set_fuel_stream_for_O2_dry
    Stream fuel_from_O2_dry = set_fuel_stream_for_O2_dry(X_O2_dry, fuel, air, tol_X);
    EXPECT_NEAR(fuel_from_O2_dry.mdot, mdot_fuel_original, mdot_fuel_original * tol_mdot_rel);
    
    // Test set_fuel_stream_for_CO2_dry
    Stream fuel_from_CO2_dry = set_fuel_stream_for_CO2_dry(X_CO2_dry, fuel, air, tol_X);
    EXPECT_NEAR(fuel_from_CO2_dry.mdot, mdot_fuel_original, mdot_fuel_original * tol_mdot_rel);
}

// Test round-trip for oxidizer solvers
TEST_F(ThermoTransportTest, InverseSolversRoundTripOxidizer) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    const std::size_t idx_CO2 = species_index_from_name("CO2");
    
    // Fuel stream: pure CH4 at 300 K, 0.3 kg/s (fixed)
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    fuel.mdot = 0.3;
    
    // Oxidizer stream: dry air at 320 K
    Stream air;
    air.state.T = 320.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    
    // Set air mdot for phi = 0.6 (lean)
    // phi = (mdot_fuel / mdot_air) / (mdot_fuel / mdot_air)_stoich
    // For CH4: stoich air/fuel ratio ~ 17.2 kg air / kg fuel
    // At phi=0.6: air/fuel = 17.2 / 0.6 = 28.7
    // mdot_air = mdot_fuel * 28.7 = 0.3 * 28.7 = 8.6 kg/s (approximate)
    // Use set_fuel_stream_for_phi in reverse to get correct air mdot
    Stream air_template = air;
    air_template.mdot = 10.0;  // temporary
    Stream fuel_at_phi = set_fuel_stream_for_phi(0.6, fuel, air_template);
    // Now scale: we want fuel.mdot = 0.3, so air.mdot = 10.0 * (0.3 / fuel_at_phi.mdot)
    double mdot_air_original = 10.0 * (0.3 / fuel_at_phi.mdot);
    air.mdot = mdot_air_original;
    
    // Mix and combust
    Stream mixed = mix({fuel, air});
    State burned = complete_combustion(mixed.state);
    std::vector<double> X_burned = burned.X;
    std::vector<double> X_dry = convert_to_dry_fractions(X_burned);
    
    double Tad = burned.T;
    double X_O2 = X_burned[idx_O2];
    double X_CO2 = X_burned[idx_CO2];
    double X_O2_dry = X_dry[idx_O2];
    double X_CO2_dry = X_dry[idx_CO2];
    
    // Verify lean combustion
    EXPECT_GT(X_O2, 0.0);
    
    const double tol_T = 1.0;
    const double tol_X = 1e-6;
    const double tol_mdot_rel = 0.01;
    
    // Test set_oxidizer_stream_for_Tad
    Stream air_from_Tad = set_oxidizer_stream_for_Tad(Tad, fuel, air, tol_T);
    EXPECT_NEAR(air_from_Tad.mdot, mdot_air_original, mdot_air_original * tol_mdot_rel);
    
    // Test set_oxidizer_stream_for_O2
    Stream air_from_O2 = set_oxidizer_stream_for_O2(X_O2, fuel, air, tol_X);
    EXPECT_NEAR(air_from_O2.mdot, mdot_air_original, mdot_air_original * tol_mdot_rel);
    
    // Test set_oxidizer_stream_for_CO2
    Stream air_from_CO2 = set_oxidizer_stream_for_CO2(X_CO2, fuel, air, tol_X);
    EXPECT_NEAR(air_from_CO2.mdot, mdot_air_original, mdot_air_original * tol_mdot_rel);
    
    // Test set_oxidizer_stream_for_O2_dry
    Stream air_from_O2_dry = set_oxidizer_stream_for_O2_dry(X_O2_dry, fuel, air, tol_X);
    EXPECT_NEAR(air_from_O2_dry.mdot, mdot_air_original, mdot_air_original * tol_mdot_rel);
    
    // Test set_oxidizer_stream_for_CO2_dry
    Stream air_from_CO2_dry = set_oxidizer_stream_for_CO2_dry(X_CO2_dry, fuel, air, tol_X);
    EXPECT_NEAR(air_from_CO2_dry.mdot, mdot_air_original, mdot_air_original * tol_mdot_rel);
}

// -------------------------------------------------------------
// Edge case tests: Verify solvers reject invalid/ambiguous inputs
// -------------------------------------------------------------

// Test that Tad solver rejects target below oxidizer temperature
TEST_F(ThermoTransportTest, SetFuelStreamForTadRejectsBelowOxidizerT) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    Stream air;
    air.state.T = 400.0;  // Hot air
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target Tad = 350 K, which is below air temperature (400 K)
    // This should throw because Tad cannot be below oxidizer T
    EXPECT_THROW(set_fuel_stream_for_Tad(350.0, fuel, air), std::invalid_argument);
}

// Test that oxidizer Tad solver rejects target below fuel temperature
TEST_F(ThermoTransportTest, SetOxidizerStreamForTadRejectsBelowFuelT) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    Stream fuel;
    fuel.state.T = 500.0;  // Hot fuel
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    fuel.mdot = 0.5;
    
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    
    // Target Tad = 400 K, which is below fuel temperature (500 K)
    EXPECT_THROW(set_oxidizer_stream_for_Tad(400.0, fuel, air), std::invalid_argument);
}

// Test that O2 solver rejects target >= oxidizer O2 (impossible)
TEST_F(ThermoTransportTest, SetFuelStreamForO2RejectsAboveOxidizerO2) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target O2 = 0.25, which is above air O2 (0.21) - impossible
    EXPECT_THROW(set_fuel_stream_for_O2(0.25, fuel, air), std::invalid_argument);
    
    // Target O2 = 0.0 - also invalid (would require stoich or rich)
    EXPECT_THROW(set_fuel_stream_for_O2(0.0, fuel, air), std::invalid_argument);
    
    // Target O2 negative - invalid
    EXPECT_THROW(set_fuel_stream_for_O2(-0.05, fuel, air), std::invalid_argument);
}

// Test that CO2 solver rejects zero or negative target
TEST_F(ThermoTransportTest, SetFuelStreamForCO2RejectsInvalidTarget) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target CO2 = 0.0 - invalid
    EXPECT_THROW(set_fuel_stream_for_CO2(0.0, fuel, air), std::invalid_argument);
    
    // Target CO2 negative - invalid
    EXPECT_THROW(set_fuel_stream_for_CO2(-0.05, fuel, air), std::invalid_argument);
}

// Round-trip test for rich combustion (phi > 1)
TEST_F(ThermoTransportTest, InverseSolversRoundTripRich) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Fuel stream: pure CH4 at 300 K
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Test at phi = 1.3 (rich)
    double phi = 1.3;
    Stream fuel_phi = set_fuel_stream_for_phi(phi, fuel, air);
    double mdot_fuel_original = fuel_phi.mdot;
    
    // Mix and combust
    Stream mixed = mix({fuel_phi, air});
    State burned = complete_combustion(mixed.state);
    
    // Extract Tad from burned products
    double Tad = burned.T;
    
    // Verify we have rich combustion (no O2 in products)
    EXPECT_LT(burned.X[idx_O2], 0.001);
    EXPECT_GT(Tad, air.T());
    
    // Test set_fuel_stream_for_Tad with lean=false recovers the original fuel mdot
    const double tol_T = 1.0;
    const double tol_mdot_rel = 0.01;  // 1% relative tolerance on mdot
    
    Stream fuel_from_Tad = set_fuel_stream_for_Tad(Tad, fuel, air, tol_T, 100, false);
    EXPECT_NEAR(fuel_from_Tad.mdot, mdot_fuel_original, mdot_fuel_original * tol_mdot_rel);
    
    // Verify the recovered fuel produces the same Tad
    Stream mixed_check = mix({fuel_from_Tad, air});
    State burned_check = complete_combustion(mixed_check.state);
    EXPECT_NEAR(burned_check.T, Tad, tol_T);
}

// Test set_fuel_stream_for_Tad with lean=false (rich side)
TEST_F(ThermoTransportTest, SetFuelStreamForTadRichSide) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    // Fuel stream: pure CH4 at 300 K
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    
    // Oxidizer stream: dry air at 300 K, 10 kg/s
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 10.0;
    
    // Target Tad = 1800 K (achievable on both lean and rich sides)
    double T_ad_target = 1800.0;
    
    // Get lean solution
    Stream fuel_lean = set_fuel_stream_for_Tad(T_ad_target, fuel, air, 1.0, 100, true);
    
    // Get rich solution
    Stream fuel_rich = set_fuel_stream_for_Tad(T_ad_target, fuel, air, 1.0, 100, false);
    
    // Verify both achieve the target Tad
    Stream mixed_lean = mix({fuel_lean, air});
    State burned_lean = complete_combustion(mixed_lean.state);
    EXPECT_NEAR(burned_lean.T, T_ad_target, 2.0);
    
    Stream mixed_rich = mix({fuel_rich, air});
    State burned_rich = complete_combustion(mixed_rich.state);
    EXPECT_NEAR(burned_rich.T, T_ad_target, 2.0);
    
    // Verify lean has O2 in products, rich has no O2
    EXPECT_GT(burned_lean.X[idx_O2], 0.01);  // Lean: O2 excess
    EXPECT_LT(burned_rich.X[idx_O2], 0.001); // Rich: no O2 (all consumed)
    
    // Verify rich has more fuel than lean
    EXPECT_GT(fuel_rich.mdot, fuel_lean.mdot);
}

// Test that solvers reject zero mdot on the fixed stream
TEST_F(ThermoTransportTest, InverseSolversRejectZeroMdot) {
    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_O2 = species_index_from_name("O2");
    const std::size_t idx_N2 = species_index_from_name("N2");
    
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n, 0.0);
    fuel.state.X[idx_CH4] = 1.0;
    fuel.mdot = 0.0;  // Zero mdot
    
    Stream air;
    air.state.T = 300.0;
    air.state.P = 101325.0;
    air.state.X = std::vector<double>(n, 0.0);
    air.state.X[idx_O2] = 0.21;
    air.state.X[idx_N2] = 0.79;
    air.mdot = 0.0;  // Zero mdot
    
    // Fuel solvers should reject zero oxidizer mdot
    EXPECT_THROW(set_fuel_stream_for_Tad(1500.0, fuel, air), std::invalid_argument);
    EXPECT_THROW(set_fuel_stream_for_O2(0.10, fuel, air), std::invalid_argument);
    EXPECT_THROW(set_fuel_stream_for_CO2(0.05, fuel, air), std::invalid_argument);
    
    // Set air mdot for oxidizer solver tests
    air.mdot = 10.0;
    fuel.mdot = 0.0;
    
    // Oxidizer solvers should reject zero fuel mdot
    EXPECT_THROW(set_oxidizer_stream_for_Tad(1500.0, fuel, air), std::invalid_argument);
    EXPECT_THROW(set_oxidizer_stream_for_O2(0.10, fuel, air), std::invalid_argument);
    EXPECT_THROW(set_oxidizer_stream_for_CO2(0.05, fuel, air), std::invalid_argument);
}
