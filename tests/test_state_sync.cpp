#include <gtest/gtest.h>
#include "../include/state.h"
#include "../include/thermo.h"
#include "../include/composition.h"
#include <vector>
#include <cmath>

using namespace combaero;

TEST(StateSyncTest, XtoY) {
    State s;
    size_t n = num_species();
    // Standard air-like composition
    std::vector<double> X(n, 0.0);
    if (n > 3) {
        X[1] = 0.78; // N2
        X[3] = 0.21; // O2
        if (n > 4) X[4] = 0.01; // AR
    } else {
        X[0] = 1.0;
    }
    s.set_X(X);

    // Check synchronization (Y should be calculated from X)
    std::vector<double> X_back = s.X;
    std::vector<double> Y = s.Y;

    std::vector<double> Y_calc = mole_to_mass(X_back);

    ASSERT_EQ(Y.size(), Y_calc.size());
    for (size_t i = 0; i < Y.size(); ++i) {
        EXPECT_NEAR(Y[i], Y_calc[i], 1e-10);
    }
}

TEST(StateSyncTest, YtoX) {
    State s;
    size_t n = num_species();
    // Pure oxygen mass fractions (assuming O2 is at index 3)
    std::vector<double> Y(n, 0.0);
    if (n > 3) Y[3] = 1.0;
    else Y[0] = 1.0;

    s.set_Y(Y);

    std::vector<double> X = s.X;
    std::vector<double> X_calc = mass_to_mole(Y);

    ASSERT_EQ(X.size(), X_calc.size());
    for (size_t i = 0; i < X.size(); ++i) {
        EXPECT_NEAR(X[i], X_calc[i], 1e-10);
    }
}

TEST(StateSyncTest, Normalization) {
    State s;
    size_t n = num_species();
    // Non-normalized mole fractions
    std::vector<double> X(n, 0.0);
    if (n > 3) {
        X[1] = 2.0;
        X[3] = 8.0;
    } else {
        X[0] = 5.0;
    }
    s.set_X(X);

    double sum = 0.0;
    for (double x : s.X) sum += x;
    EXPECT_NEAR(sum, 1.0, 1e-10);
}
