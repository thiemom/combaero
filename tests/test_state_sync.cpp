#include <gtest/gtest.h>
#include "../include/state.h"
#include "../include/thermo.h"
#include "../include/composition.h"
#include "../include/compressible.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

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

// ---------------------------------------------------------------------------
// The X/Y invariant across the library boundary (issue #352)
//
// State::X and State::Y are both public, and only set_X()/set_Y() keep them in
// sync. Every State property getter (h, cp, rho, mw, ...) reads X, so a State
// whose Y never got populated computes properties correctly and looks healthy.
// Y is read in exactly one place a caller is likely to reach -- mix() -- which
// is why the defect stayed invisible: it surfaces as a segfault in mixing,
// arbitrarily far from the code that built the state.
// ---------------------------------------------------------------------------

namespace {

std::vector<double> air_like() {
    std::vector<double> X(num_species(), 0.0);
    if (num_species() > 3) {
        X[1] = 0.78;
        X[3] = 0.21;
        if (num_species() > 4) X[4] = 0.01;
    } else {
        X[0] = 1.0;
    }
    return X;
}

}  // namespace

// A state assembled by writing the X member directly is the mistake this guard
// exists for. It must report what to do, not read off the end of an empty Y.
TEST(StateSyncTest, MixRejectsAStateWhoseYWasNeverPopulated) {
    Stream a;
    a.state.T = 500.0;
    a.state.P = 2.0e5;
    a.state.X = air_like();  // deliberate: bypasses set_X, leaves Y empty
    a.mdot = 1.0;

    ASSERT_TRUE(a.state.Y.empty());

    Stream b = a;
    EXPECT_THROW(mix({a, b}), std::invalid_argument);

    try {
        mix({a, b});
        FAIL() << "expected mix() to reject the unsynced state";
    } catch (const std::invalid_argument &e) {
        // The size mismatch alone does not tell a caller what they did wrong.
        EXPECT_NE(std::string(e.what()).find("set_X"), std::string::npos);
    }
}

// The same state built through the setter must still mix.
TEST(StateSyncTest, MixAcceptsAProperlyBuiltState) {
    Stream a;
    a.state.T = 500.0;
    a.state.P = 2.0e5;
    a.state.set_X(air_like());
    a.mdot = 1.0;

    Stream b = a;
    Stream m = mix({a, b});
    EXPECT_DOUBLE_EQ(m.mdot, 2.0);
    EXPECT_EQ(m.state.Y.size(), m.state.X.size());
}

// The library itself must not hand back a state that cannot be mixed. Every
// compressible solver assigned the X member directly and returned states with
// an empty Y, so a caller who did nothing wrong got a segfault out of mix().
TEST(StateSyncTest, CompressibleSolversReturnStatesThatCanBeMixed) {
    const std::vector<double> X = air_like();

    const FannoSolution fanno = fanno_channel(500.0, 2.0e5, 120.0, 1.0, 0.05, 0.02, X);
    EXPECT_EQ(fanno.inlet.Y.size(), fanno.inlet.X.size()) << "fanno inlet";
    EXPECT_EQ(fanno.outlet.Y.size(), fanno.outlet.X.size()) << "fanno outlet";

    const CompressibleFlowSolution nozzle = nozzle_flow(1000.0, 5.0e5, 1.0e5, 0.01, X);
    EXPECT_EQ(nozzle.stagnation.Y.size(), nozzle.stagnation.X.size()) << "nozzle stagnation";
    EXPECT_EQ(nozzle.outlet.Y.size(), nozzle.outlet.X.size()) << "nozzle outlet";

    // The point of the invariant: these mix without reading off the end of Y.
    Stream a;
    a.state = fanno.outlet;
    a.mdot = 1.0;
    Stream b;
    b.state = nozzle.outlet;
    b.mdot = 2.0;
    const Stream m = mix({a, b});
    EXPECT_DOUBLE_EQ(m.mdot, 3.0);
    EXPECT_EQ(m.state.Y.size(), m.state.X.size());
}
