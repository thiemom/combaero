// C++ Mynard junction closure against the Python that ships.
//
// The port's gate is equivalence, so every reference number here is produced
// by combaero.network._mynard2010 (see
// validation/junction/data/generate_mynard_reference.py) rather than typed in.
// The derivative rows are CENTRAL DIFFERENCES of the Python, which makes them
// an independent check on the C++ dual partials: two methods, two languages,
// no shared derivation.
//
// Cases are chosen to reach the branches, not to look tidy -- both flow
// directions, both mask orientations, equal and extreme area ratios, angles
// either side of the second-quadrant flip, a collinear lateral and one that
// just misses it, and each empirical term on and off. See the generator.
//
// Coverage, measured (llvm-cov, -fcoverage-mapping): 100% of regions and
// lines in include/mynard_junction.h, 98% of branches. The two branches not
// taken are defensive and unreachable with valid input -- the `a_max + a_min
// > 0` guard, which needs a non-positive area, and the common-port search's
// loop-exhaustion exit, which needs no common port to exist after the
// validity check has already passed. Recorded so the next reader does not
// hunt for them. Note that a percentage alone would not have been enough:
// four branches were initially unreached, and two of those WERE real gaps in
// the case list (an angle below -pi, and the joining term with one supplier).
//
// Issue #271, step 4.2/4.3.

#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "dual_number.h"
#include "math_constants.h"
#include "mynard_junction.h"
#include "validation/junction/data/mynard_reference_data.h"

using combaero::solver::DualN;
using combaero::solver::kMynardPorts;
using combaero::solver::mynard_junction_loss_coefficient;
using combaero::validation::junction::kMynardCases;
using combaero::validation::junction::MynardCase;

namespace {

using D3 = DualN<3>; // one seed per port velocity

// Seeded on the three velocities, which is what the golden derivatives vary.
combaero::solver::MynardResult<3> Evaluate(const MynardCase& c) {
  std::array<D3, kMynardPorts> u{};
  for (int i = 0; i < kMynardPorts; ++i) u[i] = D3::seed(c.u[i], i);
  return mynard_junction_loss_coefficient<3>(u, c.area, c.theta,
                                             c.joining_etransfer_alpha, c.eta_scale);
}

// Values are compared tightly; derivatives against a central difference get
// the looser tolerance the difference itself deserves.
constexpr double kValueTol = 1e-12;
constexpr double kDerivTol = 1e-6;

} // namespace

TEST(MynardJunction, EveryGoldenCaseIsValid) {
  for (const auto& c : kMynardCases) {
    EXPECT_TRUE(Evaluate(c).valid) << c.id;
  }
}

TEST(MynardJunction, LossCoefficientsMatchThePython) {
  for (const auto& c : kMynardCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    for (int i = 0; i < kMynardPorts; ++i) {
      EXPECT_NEAR(r.C[i].v, c.c[i], kValueTol)
          << c.id << " branch " << i;
    }
  }
}

TEST(MynardJunction, KCoefficientsMatchThePython) {
  for (const auto& c : kMynardCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    for (int j = 0; j < kMynardPorts - 1; ++j) {
      ASSERT_FALSE(std::isnan(c.k[j])) << c.id << ": reference K missing";
      EXPECT_NEAR(r.K[j].v, c.k[j], kValueTol) << c.id << " K entry " << j;
    }
  }
}

TEST(MynardJunction, DualPartialsOfCMatchAFiniteDifferenceOfThePython) {
  for (const auto& c : kMynardCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    for (int i = 0; i < kMynardPorts; ++i) {
      for (int j = 0; j < kMynardPorts; ++j) {
        double reference = c.dc_du[i * kMynardPorts + j];
        EXPECT_NEAR(r.C[i].d[j], reference,
                    kDerivTol * std::max(1.0, std::abs(reference)))
            << c.id << ": dC[" << i << "]/dU[" << j << "]";
      }
    }
  }
}

TEST(MynardJunction, DualPartialsOfKMatchAFiniteDifferenceOfThePython) {
  for (const auto& c : kMynardCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    for (int i = 0; i < kMynardPorts - 1; ++i) {
      for (int j = 0; j < kMynardPorts; ++j) {
        double reference = c.dk_du[i * kMynardPorts + j];
        EXPECT_NEAR(r.K[i].d[j], reference,
                    kDerivTol * std::max(1.0, std::abs(reference)))
            << c.id << ": dK[" << i << "]/dU[" << j << "]";
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Structure the golden table cannot express
// -----------------------------------------------------------------------------

TEST(MynardJunction, SupplierBranchesCarryNoC) {
  // C is defined per collector; the Python leaves supplier entries at zero and
  // the element relies on that when it maps K back onto ports.
  for (const auto& c : kMynardCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    for (int i = 0; i < kMynardPorts; ++i) {
      if (c.u[i] * c.area[i] < 0.0) continue; // collector
      EXPECT_DOUBLE_EQ(r.C[i].v, 0.0) << c.id << " supplier branch " << i;
      for (int j = 0; j < kMynardPorts; ++j) {
        EXPECT_DOUBLE_EQ(r.C[i].d[j], 0.0) << c.id << " supplier branch " << i;
      }
    }
  }
}

TEST(MynardJunction, KIsIndexedByTheNonCommonPorts) {
  // Which ports K describes changes with the flow pattern, and getting it
  // wrong is a silent relabelling rather than a numerical error -- exactly
  // the class of defect that cost this arc two rounds on the Python side.
  for (const auto& c : kMynardCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    const bool collector_is_common = (r.n_collector == 1);
    for (int j = 0; j < kMynardPorts - 1; ++j) {
      int port = r.k_port[j];
      bool port_is_collector = c.u[port] * c.area[port] < 0.0;
      EXPECT_NE(port_is_collector, collector_is_common)
          << c.id << ": K entry " << j << " points at the common port";
    }
    EXPECT_LT(r.k_port[0], r.k_port[1]) << c.id << ": K entries are not in port order";
  }
}

TEST(MynardJunction, AllInflowOrAllOutflowIsRefusedNotGuessed) {
  std::array<double, kMynardPorts> area{0.01, 0.01, 0.01};
  std::array<double, kMynardPorts> theta{0.0, M_PI, M_PI / 2.0};
  for (double sign : {1.0, -1.0}) {
    std::array<D3, kMynardPorts> u{};
    for (int i = 0; i < kMynardPorts; ++i) u[i] = D3::seed(sign * (1.0 + i), i);
    EXPECT_FALSE(mynard_junction_loss_coefficient<3>(u, area, theta).valid)
        << "sign " << sign;
  }
}

TEST(MynardJunction, TheDividingStreamlineTermFiresOnlyWhenCollinear) {
  // Pinned as a pair: the collinear case and the one that misses the
  // tolerance must differ by exactly the recovery term, or the threshold has
  // moved. Falsifying either half alone would not catch a shifted tolerance.
  const MynardCase* collinear = nullptr;
  const MynardCase* near = nullptr;
  for (const auto& c : kMynardCases) {
    if (std::string(c.id) == "div_collinear") collinear = &c;
    if (std::string(c.id) == "div_near_collinear") near = &c;
  }
  ASSERT_NE(collinear, nullptr);
  ASSERT_NE(near, nullptr);

  auto on = Evaluate(*collinear);
  auto off = Evaluate(*near);
  ASSERT_TRUE(on.valid);
  ASSERT_TRUE(off.valid);

  // The straight collector is port 1 in both, and it is the one that continues.
  bool any_lower = false;
  for (int j = 0; j < kMynardPorts - 1; ++j) {
    if (on.K[j].v < off.K[j].v - 0.1) any_lower = true;
  }
  EXPECT_TRUE(any_lower) << "the recovery term did not fire on the collinear case";
}

TEST(MynardJunction, TheJoiningAsymmetryTermIsInertWithOneSupplier) {
  // The correction is defined over SUPPLIER area asymmetry, so a dividing
  // junction must be untouched by it however extreme its areas. Pinned as a
  // pair against the alpha = 0 twin: asserting the alpha-on case alone would
  // pass even if the term were applied and happened to be small.
  const MynardCase* with_alpha = nullptr;
  const MynardCase* without = nullptr;
  for (const auto& c : kMynardCases) {
    if (std::string(c.id) == "div_alpha_on_is_inert") with_alpha = &c;
    if (std::string(c.id) == "div_area_ratio_4") without = &c;
  }
  ASSERT_NE(with_alpha, nullptr);
  ASSERT_NE(without, nullptr);
  ASSERT_GT(with_alpha->joining_etransfer_alpha, 0.0);
  ASSERT_DOUBLE_EQ(without->joining_etransfer_alpha, 0.0);

  auto on = Evaluate(*with_alpha);
  auto off = Evaluate(*without);
  ASSERT_TRUE(on.valid);
  ASSERT_TRUE(off.valid);
  for (int i = 0; i < kMynardPorts; ++i) {
    EXPECT_DOUBLE_EQ(on.C[i].v, off.C[i].v) << "branch " << i;
  }
  for (int j = 0; j < kMynardPorts - 1; ++j) {
    EXPECT_DOUBLE_EQ(on.K[j].v, off.K[j].v) << "K entry " << j;
  }
}

TEST(MynardJunction, AnAngleBelowMinusPiWrapsTheSameAsItsEquivalent) {
  // Nothing forbids a caller declaring -270 degrees, and the initial wrap has
  // a negative-fmod path that no ordinary geometry reaches. -1.5*pi and
  // +0.5*pi are the same direction, so the closure must not tell them apart.
  std::array<double, kMynardPorts> area{0.01, 0.01, 0.01};
  std::array<D3, kMynardPorts> u{D3::seed(10.0, 0), D3::seed(-6.0, 1), D3::seed(-4.0, 2)};

  std::array<double, kMynardPorts> below{0.0, M_PI, -1.5 * M_PI};
  std::array<double, kMynardPorts> equivalent{0.0, M_PI, 0.5 * M_PI};
  auto a = mynard_junction_loss_coefficient<3>(u, area, below);
  auto b = mynard_junction_loss_coefficient<3>(u, area, equivalent);

  ASSERT_TRUE(a.valid);
  ASSERT_TRUE(b.valid);
  for (int i = 0; i < kMynardPorts; ++i) {
    EXPECT_NEAR(a.C[i].v, b.C[i].v, 1e-12) << "branch " << i;
  }
  for (int j = 0; j < kMynardPorts - 1; ++j) {
    EXPECT_NEAR(a.K[j].v, b.K[j].v, 1e-12) << "K entry " << j;
  }
}
