// Whole-element (f, J) for the momentum-CV junction, against the Python.
//
// Residual values are compared with the shipping Python element, so the gate
// is equivalence with what runs today.
//
// The Jacobian is compared against CENTRAL DIFFERENCES of that residual, NOT
// against the Python's own analytic Jacobian. The Python assembles its rows
// from an explicit dKQ/dmdot block plus a hand-derived dR/dP column for the
// common port alone, and misses the other two ports' static-pressure columns
// -- K depends on every port's velocity and every velocity on its own
// density. Measured at a converged state those columns were off by 4.4e-3 and
// 2.3e-3. Whole-element seeding supplies them, so comparing against the
// Python's Jacobian would mark the C++ wrong exactly where it is right.
//
// Coverage, measured (llvm-cov): 100% of regions and lines in
// include/mpce_junction.h, 97% of branches. The one branch not taken is the
// common-port search's loop-exhaustion exit, unreachable once the validity
// check has passed. Measuring was not a formality: it found the 1e-9 dead
// band untested, which is the case the axial-reassignment comment makes a
// claim about.
//
// Issue #271, step 4.4.

#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "math_constants.h"
#include "mpce_junction.h"
#include "validation/junction/data/mpce_reference_data.h"

using combaero::solver::kMpcePorts;
using combaero::solver::kMpceRows;
using combaero::solver::kMpceSeeds;
using combaero::solver::mpce_v2_residuals_and_jacobian;
using combaero::solver::MpceGeometry;
using combaero::validation::junction::kMpceCases;
using combaero::validation::junction::MpceCase;

namespace {

combaero::solver::MpceResidualJacobian Evaluate(const MpceCase& c) {
  MpceGeometry geom;
  geom.area = c.area;
  geom.theta_rad = c.theta_rad;
  geom.port_sign = c.port_sign;
  geom.joining_etransfer_alpha = c.joining_etransfer_alpha;
  geom.eta_scale = c.eta_scale;
  return mpce_v2_residuals_and_jacobian(c.p_static, c.p_total, c.rho, c.drho_dp,
                                        c.outer_mdot, c.pt_jct, geom);
}

// Residuals are pressures of order 1e5 Pa, so a relative tolerance is the
// meaningful one.
constexpr double kResidualRelTol = 1e-9;
// Central differences of a 1e5-scale quantity; 1e-5 relative is what the
// difference itself is good for.
constexpr double kJacobianRelTol = 1e-5;

} // namespace

TEST(MpceJunction, EveryGoldenCaseIsAJunction) {
  for (const auto& c : kMpceCases) {
    EXPECT_TRUE(Evaluate(c).valid) << c.id;
  }
}

TEST(MpceJunction, ResidualsMatchThePython) {
  for (const auto& c : kMpceCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    for (int i = 0; i < kMpceRows; ++i) {
      EXPECT_NEAR(r.residual[i], c.residual[i],
                  kResidualRelTol * std::max(1.0, std::abs(c.residual[i])))
          << c.id << " row " << i;
    }
  }
}

TEST(MpceJunction, EveryJacobianEntryMatchesAFiniteDifference) {
  for (const auto& c : kMpceCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    for (int i = 0; i < kMpceRows; ++i) {
      for (int j = 0; j < kMpceSeeds; ++j) {
        double reference = c.jacobian_fd[i * kMpceSeeds + j];
        EXPECT_NEAR(r.jacobian[i][j], reference,
                    kJacobianRelTol * std::max(1.0, std::abs(reference)))
            << c.id << ": dR[" << i << "]/dx[" << j << "]";
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Structure the golden table cannot express
// -----------------------------------------------------------------------------

TEST(MpceJunction, ThePtColumnsAreTheIdentityAndTheJunctionColumnIsMinusOne) {
  // R_i = Pt_i - Pt_jct + ..., and the loss term carries no Pt dependence.
  // Cheap to state, and it would catch a seed-layout slip that the finite
  // differences would also catch but far less legibly.
  for (const auto& c : kMpceCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    for (int i = 0; i < kMpcePorts; ++i) {
      for (int k = 0; k < kMpcePorts; ++k) {
        double expected = (i == k) ? 1.0 : 0.0;
        EXPECT_DOUBLE_EQ(r.jacobian[i][kMpcePorts + k], expected)
            << c.id << " dR[" << i << "]/dPt[" << k << "]";
      }
      EXPECT_DOUBLE_EQ(r.jacobian[i][3 * kMpcePorts], -1.0) << c.id << " row " << i;
    }
  }
}

TEST(MpceJunction, TheMassRowIsTheSignedPortMap) {
  // Exactly the port signs, and nothing else: no pressure dependence at all.
  for (const auto& c : kMpceCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    const auto& row = r.jacobian[kMpcePorts];
    for (int j = 0; j < kMpcePorts; ++j) {
      EXPECT_DOUBLE_EQ(row[j], 0.0) << c.id << " mass row dP[" << j << "]";
      EXPECT_DOUBLE_EQ(row[kMpcePorts + j], 0.0) << c.id << " mass row dPt[" << j << "]";
      EXPECT_DOUBLE_EQ(row[2 * kMpcePorts + j], c.port_sign[j])
          << c.id << " mass row dmdot[" << j << "]";
    }
    EXPECT_DOUBLE_EQ(row[3 * kMpcePorts], 0.0) << c.id << " mass row dPt_jct";
  }
}

TEST(MpceJunction, TheCommonPortCarriesNoLossTerm) {
  // K is defined per NON-common port; the common row reduces to
  // Pt_common = Pt_jct. If this ever fails the port labelling has slipped,
  // which is a silent relabelling rather than a numerical error.
  for (const auto& c : kMpceCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    ASSERT_GE(r.common_port, 0);
    EXPECT_DOUBLE_EQ(r.k_per_port[r.common_port], 0.0) << c.id;
    bool other_is_nonzero = false;
    for (int i = 0; i < kMpcePorts; ++i) {
      if (i != r.common_port && r.k_per_port[i] != 0.0) other_is_nonzero = true;
    }
    EXPECT_TRUE(other_is_nonzero) << c.id << ": every K is zero, so nothing is being tested";
  }
}

TEST(MpceJunction, TheKSignFollowsTheFlowDirection) {
  // Separating: the common port is the single supplier and holds the higher
  // Pt, so the term signs +1. Joining: the common port is the single
  // collector and holds the lower Pt, so -1.
  for (const auto& c : kMpceCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    bool joining = std::string(c.id).rfind("join", 0) == 0;
    EXPECT_DOUBLE_EQ(r.k_term_sign, joining ? -1.0 : 1.0) << c.id;
  }
}

TEST(MpceJunction, TheNonCommonStaticPressureColumnsAreNotZero) {
  // The whole point of seeding the element rather than hand-assembling it.
  // The Python leaves these empty; if the C++ ever does too, this port has
  // stopped buying anything.
  for (const auto& c : kMpceCases) {
    auto r = Evaluate(c);
    ASSERT_TRUE(r.valid) << c.id;
    bool any = false;
    for (int i = 0; i < kMpcePorts; ++i) {
      for (int j = 0; j < kMpcePorts; ++j) {
        if (j == r.common_port) continue;
        if (std::abs(r.jacobian[i][j]) > 0.0) any = true;
      }
    }
    EXPECT_TRUE(any) << c.id << ": the non-common dR/dP columns are all zero";
  }
}

TEST(MpceJunction, AllInflowOrAllOutflowIsRefusedNotGuessed) {
  MpceGeometry geom;
  geom.area = {0.01, 0.01, 0.01};
  geom.theta_rad = {0.0, 0.0, M_PI / 2.0};
  geom.port_sign = {-1.0, 1.0, 1.0};
  std::array<double, kMpcePorts> p{2.0e5, 2.0e5, 2.0e5};
  std::array<double, kMpcePorts> pt{2.05e5, 2.05e5, 2.05e5};
  std::array<double, kMpcePorts> rho{2.2, 2.2, 2.2};
  std::array<double, kMpcePorts> drho{1.1e-5, 1.1e-5, 1.1e-5};

  // port_sign maps these so every port flows the same way.
  std::array<double, kMpcePorts> all_out{-1.0, 1.0, 1.0};
  std::array<double, kMpcePorts> all_in{1.0, -1.0, -1.0};
  EXPECT_FALSE(
      mpce_v2_residuals_and_jacobian(p, pt, rho, drho, all_out, 2.05e5, geom).valid);
  EXPECT_FALSE(
      mpce_v2_residuals_and_jacobian(p, pt, rho, drho, all_in, 2.05e5, geom).valid);
}

TEST(MpceJunction, APortInsideTheDeadBandLeavesTheAnglesAlone) {
  // The axial-back reassignment counts ports with a 1e-9 dead band, matching
  // the Python exactly -- which is STRICTER than the closure's own `Q < 0`
  // test. A port inside the band belongs to neither count, so neither arm
  // fires and every angle stays as declared.
  //
  // Observed as a discontinuity rather than by reaching inside: moving the
  // third port from inside the band to just outside it flips the count from
  // (1 positive, 1 negative) to (1 positive, 2 negative), which fires the
  // reassignment and re-points a port by a large angle. The residual has to
  // jump by far more than a 1e-12 kg/s change could otherwise explain.
  //
  // This state is one the ELEMENT would snap before calling in -- the kernel
  // deliberately does not snap, since that is a guard and guards stay in
  // Python. It is still a legal input here and its behaviour is pinned.
  MpceGeometry geom;
  geom.area = {0.01, 0.01, 0.01};
  geom.theta_rad = {0.0, 0.0, M_PI / 2.0};
  geom.port_sign = {-1.0, 1.0, 1.0};
  std::array<double, kMpcePorts> p{2.00e5, 1.97e5, 1.95e5};
  std::array<double, kMpcePorts> pt{2.06e5, 2.02e5, 2.00e5};
  std::array<double, kMpcePorts> rho{2.18, 2.15, 2.12};
  std::array<double, kMpcePorts> drho{1.09e-5, 1.09e-5, 1.09e-5};

  std::array<double, kMpcePorts> in_band{0.9, 0.9, 1.0e-12};
  std::array<double, kMpcePorts> outside{0.9, 0.9, 1.0e-3};
  auto banded = mpce_v2_residuals_and_jacobian(p, pt, rho, drho, in_band, 2.05e5, geom);
  auto fired = mpce_v2_residuals_and_jacobian(p, pt, rho, drho, outside, 2.05e5, geom);

  ASSERT_TRUE(banded.valid);
  ASSERT_TRUE(fired.valid);
  for (int i = 0; i < kMpceRows; ++i) {
    EXPECT_TRUE(std::isfinite(banded.residual[i])) << "row " << i;
  }
  // A 1e-3 kg/s change cannot move a 1e5-scale residual by thousands of Pa on
  // its own; the reassignment firing is what does.
  bool jumped = false;
  for (int i = 0; i < kMpcePorts; ++i) {
    if (std::abs(banded.residual[i] - fired.residual[i]) > 1.0e3) jumped = true;
  }
  EXPECT_TRUE(jumped) << "the axial reassignment did not change behaviour across the dead band";
}
