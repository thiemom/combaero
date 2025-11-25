#include <gtest/gtest.h>
#include "../include/thermo_transport.h"
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
                int ar_idx = species_index.count("AR") ? species_index_from_name("AR") : species_index_from_name("Ar");
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
            int h2o_idx = species_index_from_name("H2O");
            // Reduce others by 4% and add 4% water
            for (auto& val : humid_air) val *= 0.96;
            humid_air[h2o_idx] = 0.04;
        } catch (...) {
            // If H2O doesn't exist, just use air composition
        }
        
        // Pure water vapor
        try {
            int h2o_idx = species_index_from_name("H2O");
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
        int h2o_idx = species_index_from_name("H2O");
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
    double s_value = s(T, P, air_composition);
    EXPECT_GT(s_value, 100.0);
    EXPECT_TRUE(std::isfinite(s_value));
    
    // Test density - should be positive and reasonable for a gas
    double rho = density(T, P, air_composition);
    EXPECT_GT(rho, 0.5);  // Reasonable for any gas mixture
    EXPECT_LT(rho, 5.0);  // Reasonable upper bound
}

// Test molecular weight calculation
TEST_F(ThermoTransportTest, MolecularWeight) {
    // Molecular weight should be positive and reasonable
    double mw = mwmix(air_composition);
    EXPECT_GT(mw, 10.0);  // Lighter than any common gas mixture
    EXPECT_LT(mw, 100.0); // Heavier than most common mixtures
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

    int idx_CH4 = species_index_from_name("CH4");
    int idx_O2  = species_index_from_name("O2");

    double n_tot = n_CH4 + n_O2;
    X[idx_CH4] = n_CH4 / n_tot;
    X[idx_O2]  = n_O2  / n_tot;

    return X;
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

    int idx_CO2 = species_index_from_name("CO2");
    int idx_H2O = species_index_from_name("H2O");
    int idx_O2  = species_index_from_name("O2");
    int idx_CH4 = species_index_from_name("CH4");

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

    int idx_CO2 = species_index_from_name("CO2");
    int idx_H2O = species_index_from_name("H2O");
    int idx_O2  = species_index_from_name("O2");
    int idx_CH4 = species_index_from_name("CH4");

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

    int idx_CO2 = species_index_from_name("CO2");
    int idx_H2O = species_index_from_name("H2O");
    int idx_O2  = species_index_from_name("O2");
    int idx_CH4 = species_index_from_name("CH4");

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

    int idx_CO2 = species_index_from_name("CO2");
    int idx_H2O = species_index_from_name("H2O");
    int idx_O2  = species_index_from_name("O2");
    int idx_CH4 = species_index_from_name("CH4");

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
