// Ideal momentum junction (#517): the parameter-free closure against its own
// closed forms, and its Jacobian against central differences.
//
// The closed forms are written out here independently of the kernel:
//   - merge: Blaisdell & Manson's theory, the same algebra as
//     validation/junction/joining_variants.py::bassett_cv(wall=False);
//   - split: the vector Borda-Carnot loss rho/2 |v_j - v_s|^2 on every
//     receiver, which for a collinear equal-area receiver is (1 - q)^2.
// The kernel is evaluated incompressibly (drho/dP = 0) for those, so the only
// departure is the role smoothing, O((eps/m)^2) ~ 1e-7.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "ideal_junction.h"
#include "math_constants.h"

using combaero::solver::ideal_junction_residuals_and_jacobian;
using combaero::solver::IdealJunctionGeometry;
using combaero::solver::IdealJunctionResult;

namespace {

constexpr double kRho = 1.2;
// Role-smoothing bias on a row, relative to the dynamic head: about
// (IDEAL_ROLE_EPS_REL)^2 m_total / (4 |m_weakest|) times a velocity ratio.
// Measured <= 5e-6 over the B&M merge grid (worst: a 5% supplier at psi 16).
constexpr double kSmoothTol = 5e-5;

struct Ports {
  std::vector<double> m;  // junction convention, + out
  IdealJunctionGeometry g;
};

// Port static pressures that zero the port rows at p_jn = 0 (incompressible,
// so the rows are P_i - p_jn + g_i(m) and P_i = -g_i).
std::vector<double> solved_statics(const Ports& p) {
  const auto n = p.m.size();
  const std::vector<double> zero(n, 0.0);
  const std::vector<double> rho(n, kRho);
  const IdealJunctionResult r =
      ideal_junction_residuals_and_jacobian(zero, rho, zero, p.m, 0.0, p.g);
  std::vector<double> P(n);
  for (std::size_t i = 0; i < n; ++i) P[i] = -r.residual[i];
  return P;
}

double total_pressure(double P, double m, double A) {
  const double u = m / (kRho * A);
  return P + 0.5 * kRho * u * u;
}

// Merge: collector 0 on +x, straight inlet 1 (axis pi), lateral 2 whose FLOW
// is deflected by d from the outlet axis (its outward axis is d + pi).
Ports merge_tee(double q_lat, double psi, double d_rad) {
  Ports p;
  p.m = {1.0, -(1.0 - q_lat), -q_lat};
  p.g.area = {1.0, 1.0, 1.0 / psi};
  p.g.axis_rad = {0.0, M_PI, d_rad + M_PI};
  p.g.port_sign = {1.0, 1.0, 1.0};
  return p;
}

// Split: supply 0 (axis pi), straight receiver 1 (axis 0), branch 2 at theta.
Ports split_tee(double q_straight, double psi, double theta) {
  Ports p;
  p.m = {-1.0, q_straight, 1.0 - q_straight};
  p.g.area = {1.0, 1.0, 1.0 / psi};
  p.g.axis_rad = {M_PI, 0.0, theta};
  p.g.port_sign = {1.0, 1.0, 1.0};
  return p;
}

}  // namespace

// Merge: K_j = (Pt_j - Pt_c) / (rho u_c^2 / 2) equals B&M's theory,
//   dp = 2 (A_c - M) / A_c,  M = sum_j (U_j/u_c)^2 A_j cos d_j,
//   K_j = dp + (U_j/u_c)^2 - 1.
TEST(IdealJunction, MergeReducesToBlaisdellMansonTheory) {
  for (double psi : {1.0, 2.25, 4.0}) {
    for (double d_deg : {30.0, 45.0, 60.0, 90.0, 120.0, 150.0}) {
      for (double q : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        const double d = d_deg * M_PI / 180.0;
        const Ports p = merge_tee(q, psi, d);
        const std::vector<double> P = solved_statics(p);
        const double uc = 1.0 / kRho;  // m_c = 1, A_c = 1
        const double head = 0.5 * kRho * uc * uc;
        const double Ptc = total_pressure(P[0], p.m[0], p.g.area[0]);

        const double Us = (1.0 - q) / kRho;
        const double Ul = q / (kRho * p.g.area[2]);
        const double Mm = (Us / uc) * (Us / uc) * 1.0 +
                          (Ul / uc) * (Ul / uc) * p.g.area[2] * std::cos(d);
        const double dp = 2.0 * (1.0 - Mm) / 1.0;
        const double K_s = dp + (Us / uc) * (Us / uc) - 1.0;
        const double K_l = dp + (Ul / uc) * (Ul / uc) - 1.0;

        const double Pts = total_pressure(P[1], p.m[1], p.g.area[1]);
        const double Ptl = total_pressure(P[2], p.m[2], p.g.area[2]);
        EXPECT_NEAR((Pts - Ptc) / head, K_s, kSmoothTol)
            << "psi " << psi << " d " << d_deg << " q " << q;
        EXPECT_NEAR((Ptl - Ptc) / head, K_l, kSmoothTol)
            << "psi " << psi << " d " << d_deg << " q " << q;
        // Suppliers share the junction static pressure, up to the smoothing.
        EXPECT_NEAR(P[1], 0.0, kSmoothTol * head);
        EXPECT_NEAR(P[2], 0.0, kSmoothTol * head);
      }
    }
  }
}

// Split: every receiver loses the vector Borda head, never negative; the
// collinear equal-area receiver gives (1 - q)^2 on the straight fraction.
TEST(IdealJunction, SplitIsVectorBordaCarnot) {
  for (double psi : {1.0, 3.0}) {
    for (double th_deg : {30.0, 45.0, 60.0, 90.0, 120.0}) {
      for (double qs : {0.0, 0.2, 0.5, 0.8, 1.0}) {
        const double th = th_deg * M_PI / 180.0;
        const Ports p = split_tee(qs, psi, th);
        const std::vector<double> P = solved_statics(p);
        const double u0 = 1.0 / kRho;
        const double head = 0.5 * kRho * u0 * u0;
        const double Pt0 = total_pressure(P[0], p.m[0], p.g.area[0]);
        const double v0x = u0;  // arriving along +x

        for (int j : {1, 2}) {
          const double uj = p.m[j] / (kRho * p.g.area[j]);
          const double dx = uj * std::cos(p.g.axis_rad[j]) - v0x;
          const double dy = uj * std::sin(p.g.axis_rad[j]);
          const double K_borda = (dx * dx + dy * dy) / (u0 * u0);
          const double Ptj = total_pressure(P[j], p.m[j], p.g.area[j]);
          const double K = (Pt0 - Ptj) / head;
          EXPECT_NEAR(K, K_borda, kSmoothTol)
              << "port " << j << " psi " << psi << " th " << th_deg << " qs " << qs;
          EXPECT_GE(K, -1e-6);
        }
        const double Pt1 = total_pressure(P[1], p.m[1], p.g.area[1]);
        EXPECT_NEAR((Pt0 - Pt1) / head, (1.0 - qs) * (1.0 - qs), kSmoothTol);
      }
    }
  }
}

// A port carrying no flow sees the junction static pressure.
TEST(IdealJunction, DeadPortSitsAtJunctionStatic) {
  const Ports p = split_tee(1.0, 1.0, M_PI / 2.0);  // branch carries nothing
  const std::vector<double> P = solved_statics(p);
  const double u0 = 1.0 / kRho;
  EXPECT_NEAR(P[2], 0.0, 1e-3 * 0.5 * kRho * u0 * u0);
}

// Merge energy: the mechanical power lost, sum_in |m| Pt - sum_out m Pt (all
// over rho), is never negative, although a weak inlet can gain Pt (ejector
// action) -- the verifier's admissibility check must accept every merge.
TEST(IdealJunction, MergeDissipatesMechanicalEnergy) {
  for (double psi : {0.5, 1.0, 2.0, 4.0, 9.0}) {
    for (double d_deg : {0.0, 15.0, 30.0, 45.0, 60.0, 90.0, 120.0, 150.0, 180.0}) {
      for (double q = 0.02; q < 1.0; q += 0.04) {
        const Ports p = merge_tee(q, psi, d_deg * M_PI / 180.0);
        const std::vector<double> P = solved_statics(p);
        double loss = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
          loss -= p.m[i] * total_pressure(P[i], p.m[i], p.g.area[i]);
        }
        EXPECT_GE(loss, -1e-9) << "psi " << psi << " d " << d_deg << " q " << q;
      }
    }
  }
}

namespace {

// Central differences of the full residual, with rho following P linearly as
// the kernel's first-order seeding assumes.
void expect_jacobian_matches_fd(const std::vector<double>& P0,
                                const std::vector<double>& rho0,
                                const std::vector<double>& drho,
                                const std::vector<double>& m0, double pjn0,
                                const IdealJunctionGeometry& g) {
  const auto n = P0.size();
  const IdealJunctionResult r =
      ideal_junction_residuals_and_jacobian(P0, rho0, drho, m0, pjn0, g);
  auto eval = [&](std::size_t col, double h) {
    std::vector<double> P = P0, m = m0, rho = rho0;
    double pjn = pjn0;
    if (col < n) {
      P[col] += h;
      rho[col] += drho[col] * h;
    } else if (col < 2 * n) {
      m[col - n] += h;
    } else {
      pjn += h;
    }
    return ideal_junction_residuals_and_jacobian(P, rho, drho, m, pjn, g).residual;
  };
  // Row scales: a pressure row by the largest port dynamic head, the mass row
  // by the throughput.
  double head = 1.0, mscale = 1e-12;
  for (std::size_t i = 0; i < n; ++i) {
    const double u = m0[i] / (rho0[i] * g.area[i]);
    head = std::max(head, 0.5 * rho0[i] * u * u);
    mscale = std::max(mscale, std::abs(m0[i]));
  }
  for (std::size_t c = 0; c < 2 * n + 1; ++c) {
    const double x = c < n ? P0[c] : (c < 2 * n ? m0[c - n] : pjn0);
    // Flow steps well inside the role blend (width 1e-4 of the throughput),
    // or central differences carry its curvature as truncation error. The
    // loss-factor term (#526) adds curvature: at 1e-7 / 1e-8 of the throughput
    // the worst excess was 8.8e-4 / 8.8e-6 -- the h^2 of truncation, not a
    // Jacobian error -- so the step is 1e-9.
    const double h = c < 2 * n && c >= n ? 1e-9 * std::max(mscale, 1e-6)
                                         : 1e-6 * std::max(std::abs(x), 1.0);
    const std::vector<double> rp = eval(c, h), rm = eval(c, -h);
    // Compare J_ij * |dx_j| for a unit change of the column's natural scale.
    const double xs = c < 2 * n && c >= n ? mscale : head;
    for (std::size_t row = 0; row <= n; ++row) {
      const double fd = (rp[row] - rm[row]) / (2.0 * h);
      const double rs = row < n ? head : mscale;
      EXPECT_NEAR(r.jacobian[row][c] * xs / rs, fd * xs / rs, 1e-6)
          << "row " << row << " col " << c;
    }
  }
}

}  // namespace

// Compressible states (rho moving with P), merge, split and a port reversing.
// Areas are shrunk so the ports run at 20-60 m/s: with heads of ~1 Pa on 2-bar
// statics the central differences would measure round-off, not the Jacobian.
TEST(IdealJunction, JacobianMatchesCentralDifferences) {
  const std::vector<double> P = {2.00e5, 2.01e5, 1.99e5};
  const std::vector<double> rho = {2.3, 2.31, 2.29};
  std::vector<double> drho(3);
  for (int i = 0; i < 3; ++i) drho[i] = rho[i] / P[i];
  auto realistic = [](Ports p) {
    for (double& a : p.g.area) a *= 0.02;
    return p;
  };

  const Ports merge = realistic(merge_tee(0.35, 2.0, 60.0 * M_PI / 180.0));
  expect_jacobian_matches_fd(P, rho, drho, merge.m, 2.005e5, merge.g);

  const Ports split = realistic(split_tee(0.6, 3.0, 90.0 * M_PI / 180.0));
  expect_jacobian_matches_fd(P, rho, drho, split.m, 2.005e5, split.g);

  // Through the role blend (width 1e-4 of the throughput): the branch carries
  // nothing, a fraction of the width either way, and well outside it.
  for (double s : {-1e-4, -1e-5, 0.0, 1e-5, 1e-4, 3e-3}) {
    const Ports p = realistic(split_tee(1.0 - s, 1.0, 45.0 * M_PI / 180.0));
    expect_jacobian_matches_fd(P, rho, drho, p.m, 2.005e5, p.g);
  }

  // Port signs map the connecting element's flow.
  Ports flipped = split;
  flipped.g.port_sign = {-1.0, 1.0, 1.0};
  flipped.m[0] = -flipped.m[0];
  expect_jacobian_matches_fd(P, rho, drho, flipped.m, 2.005e5, flipped.g);
}

// The all-zero state is an ordinary point: finite rows, P_i = p_jn, and a
// nonsingular pressure block.
TEST(IdealJunction, AllZeroFlowIsRegular) {
  const Ports p = split_tee(0.5, 1.0, M_PI / 2.0);
  const std::vector<double> P = {1e5, 1e5, 1e5}, rho = {1.2, 1.2, 1.2};
  const std::vector<double> drho = {1.2e-5, 1.2e-5, 1.2e-5};
  const std::vector<double> m = {0.0, 0.0, 0.0};
  const IdealJunctionResult r =
      ideal_junction_residuals_and_jacobian(P, rho, drho, m, 1e5, p.g);
  for (int i = 0; i < 4; ++i) {
    EXPECT_TRUE(std::isfinite(r.residual[i]));
    EXPECT_NEAR(r.residual[i], 0.0, 1e-9);
  }
  for (int i = 0; i < 3; ++i) EXPECT_DOUBLE_EQ(r.jacobian[i][i], 1.0);
  expect_jacobian_matches_fd(P, rho, drho, m, 1e5, p.g);
}

TEST(IdealJunction, RejectsInconsistentPorts) {
  IdealJunctionGeometry g;
  g.area = {1.0, 1.0};
  g.axis_rad = {0.0, M_PI};
  g.port_sign = {1.0, 1.0};
  const std::vector<double> two = {1.0, 1.0};
  const std::vector<double> three = {1.0, 1.0, 1.0};
  EXPECT_THROW(ideal_junction_residuals_and_jacobian(three, three, three, three, 0.0, g),
               std::invalid_argument);
  const std::vector<double> one = {1.0};
  g.area = {1.0};
  g.axis_rad = {0.0};
  g.port_sign = {1.0};
  EXPECT_THROW(ideal_junction_residuals_and_jacobian(one, one, one, one, 0.0, g),
               std::invalid_argument);
  (void)two;
}

// ---------------------------------------------------------------------------
// Receiver loss factor f (#526): Pt_j = Pt_mix,j - f_j L_j.
// ---------------------------------------------------------------------------

namespace {

std::vector<double> solved_statics_f(Ports p, const std::vector<double>& f) {
  p.g.loss_factor = f;
  return solved_statics(p);
}

double rel2(const Ports& p, std::size_t i, double Vx, double Vy) {
  const double u = p.m[i] / (kRho * p.g.area[i]);
  const double dx = u * std::cos(p.g.axis_rad[i]) - Vx;
  const double dy = u * std::sin(p.g.axis_rad[i]) - Vy;
  return dx * dx + dy * dy;
}

}  // namespace

// f = 1 is the closure itself: not "close", identical, in rows and Jacobian.
TEST(IdealJunctionLossFactor, UnitFactorIsBitIdentical) {
  const std::vector<double> P = {2.00e5, 2.01e5, 1.99e5};
  const std::vector<double> rho = {2.3, 2.31, 2.29};
  std::vector<double> drho(3);
  for (int i = 0; i < 3; ++i) drho[i] = rho[i] / P[i];
  for (Ports p : {merge_tee(0.35, 2.0, 1.0), split_tee(0.6, 3.0, 1.2),
                  split_tee(1.0 + 1e-5, 1.0, 0.8)}) {
    for (double& a : p.g.area) a *= 0.02;
    const IdealJunctionResult base =
        ideal_junction_residuals_and_jacobian(P, rho, drho, p.m, 2.005e5, p.g);
    p.g.loss_factor = {1.0, 1.0, 1.0};
    const IdealJunctionResult one =
        ideal_junction_residuals_and_jacobian(P, rho, drho, p.m, 2.005e5, p.g);
    for (std::size_t r = 0; r < base.residual.size(); ++r) {
      EXPECT_EQ(one.residual[r], base.residual[r]);
      for (std::size_t c = 0; c < base.jacobian[r].size(); ++c) {
        EXPECT_EQ(one.jacobian[r][c], base.jacobian[r][c]);
      }
    }
  }
}

// Split: one supplier, no mixing, so each receiver loses exactly f_j times its
// vector Borda head -- per port, with unequal factors.
TEST(IdealJunctionLossFactor, SplitLossScalesPerPort) {
  for (double th_deg : {45.0, 90.0}) {
    for (double qs : {0.2, 0.5, 0.8}) {
      for (const std::vector<double>& f :
           {std::vector<double>{1.0, 0.0, 0.0}, std::vector<double>{0.0, 0.5, 1.5},
            std::vector<double>{2.5, 2.0, 0.3}}) {
        const Ports p = split_tee(qs, 2.0, th_deg * M_PI / 180.0);
        const std::vector<double> P = solved_statics_f(p, f);
        const double u0 = 1.0 / kRho;
        const double head = 0.5 * kRho * u0 * u0;
        const double Pt0 = total_pressure(P[0], p.m[0], p.g.area[0]);
        for (std::size_t j : {1u, 2u}) {
          const double Ptj = total_pressure(P[j], p.m[j], p.g.area[j]);
          const double borda = 0.5 * kRho * rel2(p, j, u0, 0.0);
          EXPECT_NEAR((Pt0 - Ptj) / head, f[j] * borda / head, kSmoothTol)
              << "port " << j << " th " << th_deg << " qs " << qs << " f " << f[j];
        }
        // The supplier sits at p_jn whatever ITS factor: a port keeps its f when
        // it reverses, and the factor must act only while it receives.
        EXPECT_NEAR(P[0], 0.0, kSmoothTol * head);
      }
    }
  }
}

// Merge: the dissipation is the suppliers' mixing term plus the SCALED
// receiver term, so it stays >= 0 for every f >= 0 (f = 0: mixing only).
TEST(IdealJunctionLossFactor, MergeDissipationIsMixingPlusScaledReceiverLoss) {
  for (double psi : {1.0, 2.0, 4.0}) {
    for (double d_deg : {30.0, 90.0, 150.0}) {
      for (double q : {0.1, 0.5, 0.9}) {
        for (double fc : {0.0, 0.5, 1.0, 1.5}) {
          const Ports p = merge_tee(q, psi, d_deg * M_PI / 180.0);
          const std::vector<double> P = solved_statics_f(p, {fc, 0.3, 2.0});
          // V from the suppliers (ports 1 and 2), mass-weighted.
          double Vx = 0.0, Vy = 0.0;
          for (std::size_t k : {1u, 2u}) {
            const double u = p.m[k] / (kRho * p.g.area[k]);
            Vx += -p.m[k] * u * std::cos(p.g.axis_rad[k]);
            Vy += -p.m[k] * u * std::sin(p.g.axis_rad[k]);
          }
          double phi = 0.0;
          for (std::size_t i = 0; i < 3; ++i) {
            phi -= p.m[i] * total_pressure(P[i], p.m[i], p.g.area[i]);
          }
          const double mix =
              0.5 * kRho * (-p.m[1] * rel2(p, 1, Vx, Vy) - p.m[2] * rel2(p, 2, Vx, Vy));
          const double rec = fc * p.m[0] * 0.5 * kRho * rel2(p, 0, Vx, Vy);
          const double head = 0.5 * kRho / (kRho * kRho);
          EXPECT_NEAR(phi / head, (mix + rec) / head, 10.0 * kSmoothTol)
              << "psi " << psi << " d " << d_deg << " q " << q << " f " << fc;
          EXPECT_GE(phi, -1e-9);
          // Suppliers (factors 0.3 and 2.0, both with v_k != V) stay at p_jn.
          EXPECT_NEAR(P[1], 0.0, kSmoothTol * head);
          EXPECT_NEAR(P[2], 0.0, kSmoothTol * head);
        }
      }
    }
  }
}

// Jacobian with factors: merge, split, unequal per-port factors, and through
// the role blend where a port turns from receiver to supplier.
TEST(IdealJunctionLossFactor, JacobianMatchesCentralDifferences) {
  const std::vector<double> P = {2.00e5, 2.01e5, 1.99e5};
  const std::vector<double> rho = {2.3, 2.31, 2.29};
  std::vector<double> drho(3);
  for (int i = 0; i < 3; ++i) drho[i] = rho[i] / P[i];
  auto realistic = [](Ports p, const std::vector<double>& f) {
    for (double& a : p.g.area) a *= 0.02;
    p.g.loss_factor = f;
    return p;
  };
  for (const std::vector<double>& f :
       {std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{0.5, 0.5, 0.5},
        std::vector<double>{1.5, 1.5, 1.5}, std::vector<double>{0.2, 1.7, 0.6}}) {
    const Ports merge = realistic(merge_tee(0.35, 2.0, 60.0 * M_PI / 180.0), f);
    expect_jacobian_matches_fd(P, rho, drho, merge.m, 2.005e5, merge.g);
    const Ports split = realistic(split_tee(0.6, 3.0, 90.0 * M_PI / 180.0), f);
    expect_jacobian_matches_fd(P, rho, drho, split.m, 2.005e5, split.g);
    for (double s : {-1e-4, -1e-5, 0.0, 1e-5, 1e-4, 3e-3}) {
      const Ports p = realistic(split_tee(1.0 - s, 1.0, 45.0 * M_PI / 180.0), f);
      expect_jacobian_matches_fd(P, rho, drho, p.m, 2.005e5, p.g);
    }
  }
}

TEST(IdealJunctionLossFactor, RejectsInvalidFactors) {
  Ports p = split_tee(0.5, 1.0, M_PI / 2.0);
  const std::vector<double> P = {1e5, 1e5, 1e5}, rho = {1.2, 1.2, 1.2};
  const std::vector<double> z = {0.0, 0.0, 0.0};
  for (const std::vector<double>& f :
       {std::vector<double>{1.0, -0.1, 1.0}, std::vector<double>{1.0, NAN, 1.0},
        std::vector<double>{1.0, INFINITY, 1.0}, std::vector<double>{1.0, 1.0}}) {
    p.g.loss_factor = f;
    EXPECT_THROW(ideal_junction_residuals_and_jacobian(P, rho, z, p.m, 1e5, p.g),
                 std::invalid_argument);
  }
}
