#include "tee_junction.h"
#include "math_constants.h"
#include "solver_interface.h"
#include "composition.h"
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

using namespace combaero;
using namespace combaero::solver;

// Helper: central FD for scalar functions of q
static double fd_dK5(double q, double eps = 1e-7) {
    return (K5(q + eps) - K5(q - eps)) / (2.0 * eps);
}
static double fd_dK(double (*K)(double, double, double), double q, double psi,
                    double theta, double eps = 1e-7) {
    return (K(q + eps, psi, theta) - K(q - eps, psi, theta)) / (2.0 * eps);
}

// -----------------------------------------------------------------------------
// Test 1: K5 independence of theta and psi (Bassett Eq. 15)
// -----------------------------------------------------------------------------
TEST(TeeJunctionKCoefficients, K5IndependentOfThetaAndPsi) {
    // K5 = q^2 - 1.5*q + 0.5. At q=0.5: 0.25 - 0.75 + 0.5 = 0.0
    // Note: spec validation comment "== -0.125" contains an arithmetic error;
    // the correct value from Eq. 15 is 0.0.
    const double expected = 0.25 - 1.5 * 0.5 + 0.5;
    const double thetas[] = {30.0, 45.0, 60.0, 90.0};
    const double psis[] = {1.0, 2.0, 3.0};
    for (double deg : thetas) {
        for (double psi : psis) {
            (void)psi;
            double theta = deg * M_PI / 180.0;
            (void)theta;
            EXPECT_NEAR(K5(0.5), expected, 1e-12);
        }
    }
    // Boundary values
    EXPECT_NEAR(K5(0.0), 0.5, 1e-12);
    EXPECT_NEAR(K5(1.0), 0.0, 1e-12);
}

// -----------------------------------------------------------------------------
// Test 2: K6 at theta=90deg, psi=1 (Bassett Fig. 7a reference ~0.87)
// -----------------------------------------------------------------------------
TEST(TeeJunctionKCoefficients, K6AtNinetyDegreeEqualArea) {
    const double theta = M_PI / 2.0;
    const double psi = 1.0;
    const double c = std::cos(0.75 * theta);
    const double expected = 0.25 + 1.0 - 2.0 * 0.5 * psi * c;
    EXPECT_NEAR(K6(0.5, psi, theta), expected, 1e-10);
    EXPECT_NEAR(K6(0.5, psi, theta), 0.8673, 5e-4);
}

// -----------------------------------------------------------------------------
// Test 3: K12 at theta=90deg, psi=1 (Bassett Fig. 10b reference ~0.35)
// -----------------------------------------------------------------------------
TEST(TeeJunctionKCoefficients, K12AtNinetyDegreeEqualArea) {
    const double theta = M_PI / 2.0;
    const double psi = 1.0;
    const double c = std::cos(0.75 * theta);
    const double D = psi + 0.5 * c;
    const double N = 1.0 - 0.25 - 0.25 * psi * c;
    const double expected = (2.0 * psi / D) * N + 0.25 * psi * psi - 1.0;
    EXPECT_NEAR(K12(0.5, psi, theta), expected, 1e-10);
    EXPECT_NEAR(K12(0.5, psi, theta), 0.348, 5e-3);
}

// -----------------------------------------------------------------------------
// Test 4: Jacobian continuity through q=0 (C1 check via FD)
// -----------------------------------------------------------------------------
TEST(TeeJunctionBlend, JacobianContinuityThroughZero) {
    const double theta = M_PI / 2.0;
    const double psi = 1.0;
    const double blend_k = 30.0;
    const double eps_fd = 1e-6;
    const double tol = 1e-5;

    for (int i = -20; i <= 20; ++i) {
        double q = i * 0.01;

        double dKs_anal = merging_tee_dK_straight_dq(q, psi, theta, blend_k);
        double Ks_p = merging_tee_K_straight(q + eps_fd, psi, theta, blend_k);
        double Ks_m = merging_tee_K_straight(q - eps_fd, psi, theta, blend_k);
        EXPECT_NEAR(dKs_anal, (Ks_p - Ks_m) / (2.0 * eps_fd), tol)
            << "straight dK/dq at q=" << q;

        double dKb_anal = merging_tee_dK_branch_dq(q, psi, theta, blend_k);
        double Kb_p = merging_tee_K_branch(q + eps_fd, psi, theta, blend_k);
        double Kb_m = merging_tee_K_branch(q - eps_fd, psi, theta, blend_k);
        EXPECT_NEAR(dKb_anal, (Kb_p - Kb_m) / (2.0 * eps_fd), tol)
            << "branch dK/dq at q=" << q;
    }
}

// -----------------------------------------------------------------------------
// Test 5: Blend agreement between merging and branching at q=0
// -----------------------------------------------------------------------------
TEST(TeeJunctionBlend, AgreementAtZeroFlow) {
    const double tol = 1e-10;
    const double thetas[] = {M_PI / 6.0, M_PI / 4.0, M_PI / 3.0, M_PI / 2.0};
    const double psis[] = {0.5, 1.0, 2.0};

    for (double theta : thetas) {
        for (double psi : psis) {
            EXPECT_NEAR(merging_tee_K_straight(0.0, psi, theta),
                        branching_tee_K_straight(0.0, psi, theta), tol)
                << "K_straight psi=" << psi << " theta=" << theta;
            EXPECT_NEAR(merging_tee_K_branch(0.0, psi, theta),
                        branching_tee_K_branch(0.0, psi, theta), tol)
                << "K_branch psi=" << psi << " theta=" << theta;
        }
    }
}

// -----------------------------------------------------------------------------
// Test 6: Post-solve topology flag
// -----------------------------------------------------------------------------
TEST(TeeJunctionDiagnostics, TopologyValidFlag) {
    EXPECT_FALSE(tee_topology_valid(-0.15, 0.05));
    EXPECT_TRUE(tee_topology_valid(-0.15, 0.20));
    EXPECT_TRUE(tee_topology_valid(0.5));
    EXPECT_FALSE(tee_topology_valid(1.1));
    EXPECT_TRUE(tee_topology_valid(1.1, 0.15));
}

// -----------------------------------------------------------------------------
// Input validity check
// -----------------------------------------------------------------------------
TEST(TeeJunctionDiagnostics, InputValidityCheck) {
    auto ok = tee_check_inputs(0.5, 1.0, M_PI / 2.0);
    EXPECT_TRUE(ok.valid());
    EXPECT_EQ(ok.status(), CorrelationStatus::Valid);

    auto bad_q = tee_check_inputs(1.5, 1.0, M_PI / 2.0);
    EXPECT_FALSE(bad_q.q_in_range);
    EXPECT_EQ(bad_q.status(), CorrelationStatus::Extrapolated);

    auto bad_psi = tee_check_inputs(0.5, 0.01, M_PI / 2.0);
    EXPECT_FALSE(bad_psi.psi_valid);

    auto bad_theta = tee_check_inputs(0.5, 1.0, M_PI);
    EXPECT_FALSE(bad_theta.theta_valid);
}

// -----------------------------------------------------------------------------
// Robustness: soft bounds prevent NaN/Inf outside valid range
// -----------------------------------------------------------------------------
TEST(TeeJunctionBlend, NoNaNOutsideValidRange) {
    // Very small psi (below TEE_PSI_MIN)
    double Ks = merging_tee_K_straight(0.5, 1e-6, M_PI / 2.0);
    double Kb = merging_tee_K_branch(0.5, 1e-6, M_PI / 2.0);
    EXPECT_TRUE(std::isfinite(Ks)) << "K_straight finite for tiny psi";
    EXPECT_TRUE(std::isfinite(Kb)) << "K_branch finite for tiny psi";

    // theta outside valid range
    double Ks2 = merging_tee_K_straight(0.5, 1.0, M_PI);
    double Kb2 = merging_tee_K_branch(0.5, 1.0, M_PI);
    EXPECT_TRUE(std::isfinite(Ks2)) << "K_straight finite for theta=pi";
    EXPECT_TRUE(std::isfinite(Kb2)) << "K_branch finite for theta=pi";

    // Large |q| (flow reversal extrapolation)
    double Ks3 = merging_tee_K_straight(-5.0, 1.0, M_PI / 2.0);
    double Ks4 = merging_tee_K_straight(5.0, 1.0, M_PI / 2.0);
    EXPECT_TRUE(std::isfinite(Ks3));
    EXPECT_TRUE(std::isfinite(Ks4));
}

// -----------------------------------------------------------------------------
// Analytical dK/dq vs FD for all raw K functions at 5 test points
// -----------------------------------------------------------------------------
TEST(TeeJunctionKCoefficients, AnalyticalDerivativesMatchFD) {
    const double psi = 1.5;
    const double theta = M_PI / 3.0;
    const double eps = 1e-7;
    const double tol = 1e-5;
    const double test_q[] = {0.0, 0.25, 0.5, 0.75, 1.0};

    for (double q : test_q) {
        EXPECT_NEAR(dK5_dq(q), fd_dK5(q, eps), tol) << "dK5/dq at q=" << q;
        EXPECT_NEAR(dK6_dq(q, psi, theta), fd_dK(K6, q, psi, theta, eps), tol)
            << "dK6/dq at q=" << q;
        EXPECT_NEAR(dK1_dq(q, psi, theta), fd_dK(K1, q, psi, theta, eps), tol)
            << "dK1/dq at q=" << q;
        EXPECT_NEAR(dK3_dq(q, psi, theta), fd_dK(K3, q, psi, theta, eps), tol)
            << "dK3/dq at q=" << q;
        EXPECT_NEAR(dK4_dq(q, psi, theta), fd_dK(K4, q, psi, theta, eps), tol)
            << "dK4/dq at q=" << q;
        EXPECT_NEAR(dK11_dq(q, psi, theta), fd_dK(K11, q, psi, theta, eps), tol)
            << "dK11/dq at q=" << q;
        EXPECT_NEAR(dK12_dq(q, psi, theta), fd_dK(K12, q, psi, theta, eps), tol)
            << "dK12/dq at q=" << q;
        EXPECT_NEAR(dK7_dq(q, psi, theta), fd_dK(K7, q, psi, theta, eps), tol)
            << "dK7/dq at q=" << q;
        EXPECT_NEAR(dK8_dq(q, psi, theta), fd_dK(K8, q, psi, theta, eps), tol)
            << "dK8/dq at q=" << q;
        EXPECT_NEAR(dK9_dq(q, psi, theta), fd_dK(K9, q, psi, theta, eps), tol)
            << "dK9/dq at q=" << q;
        EXPECT_NEAR(dK10_dq(q, psi, theta), fd_dK(K10, q, psi, theta, eps), tol)
            << "dK10/dq at q=" << q;
    }
}

// Soft lower bound returns lo for x << lo, x for x >> lo, is smooth.
TEST(TeeJunctionBlend, SoftLowerBoundProperties) {
    // For large x: soft_lower(x, lo) ~ x
    EXPECT_NEAR(soft_lower(1.0, 0.05), 1.0, 1e-4);
    EXPECT_NEAR(soft_lower(10.0, 0.05), 10.0, 1e-6);

    // For x = lo: result > lo (smooth, above lo)
    double at_lo = soft_lower(TEE_PSI_MIN, TEE_PSI_MIN);
    EXPECT_GT(at_lo, TEE_PSI_MIN);
    EXPECT_LT(at_lo, TEE_PSI_MIN * 1.1);

    // For x << lo: result ~ lo
    EXPECT_NEAR(soft_lower(0.0, TEE_PSI_MIN), TEE_PSI_MIN, TEE_PSI_MIN * 0.05);
    EXPECT_NEAR(soft_lower(-1.0, TEE_PSI_MIN), TEE_PSI_MIN, TEE_PSI_MIN * 0.01);
}

// -----------------------------------------------------------------------------
// Solver interface: m_dot Jacobian verification via FD
// -----------------------------------------------------------------------------
static std::vector<double> dry_air_Y() {
    std::vector<double> Y(15, 0.0);
    Y[0] = 0.767; // N2
    Y[1] = 0.233; // O2
    return Y;
}

TEST(TeeJunctionSolverInterface, MergingTeeJacobiansMdot) {
    const double m_dot_com = 0.5;
    const double m_dot_branch = 0.2;
    const double dP0_s = 150.0;
    const double dP0_b = 200.0;
    const double P_static = 2e5;
    const double T = 600.0;
    auto Y = dry_air_Y();
    const double theta = M_PI / 2.0;
    const double psi = 1.0;
    const double F_C = 0.01;
    const double eps_m = 1e-5;

    auto res0 = merging_tee_residuals_and_jacobian(
        m_dot_com, m_dot_branch, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);

    // FD for m_dot_com
    auto r_mc_p = merging_tee_residuals_and_jacobian(
        m_dot_com + eps_m, m_dot_branch, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);
    auto r_mc_m = merging_tee_residuals_and_jacobian(
        m_dot_com - eps_m, m_dot_branch, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);
    double fd_Rs_mcom = (r_mc_p.R_straight - r_mc_m.R_straight) / (2.0 * eps_m);
    double fd_Rb_mcom = (r_mc_p.R_branch - r_mc_m.R_branch) / (2.0 * eps_m);

    EXPECT_NEAR(res0.dR_straight_d_mdot_com, fd_Rs_mcom,
                std::abs(fd_Rs_mcom) * 1e-4 + 0.5);
    EXPECT_NEAR(res0.dR_branch_d_mdot_com, fd_Rb_mcom,
                std::abs(fd_Rb_mcom) * 1e-4 + 0.5);

    // FD for m_dot_branch
    auto r_mb_p = merging_tee_residuals_and_jacobian(
        m_dot_com, m_dot_branch + eps_m, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);
    auto r_mb_m = merging_tee_residuals_and_jacobian(
        m_dot_com, m_dot_branch - eps_m, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);
    double fd_Rs_mb = (r_mb_p.R_straight - r_mb_m.R_straight) / (2.0 * eps_m);
    double fd_Rb_mb = (r_mb_p.R_branch - r_mb_m.R_branch) / (2.0 * eps_m);

    EXPECT_NEAR(res0.dR_straight_d_mdot_branch, fd_Rs_mb,
                std::abs(fd_Rs_mb) * 1e-4 + 0.5);
    EXPECT_NEAR(res0.dR_branch_d_mdot_branch, fd_Rb_mb,
                std::abs(fd_Rb_mb) * 1e-4 + 0.5);
}

TEST(TeeJunctionSolverInterface, BranchingTeeJacobiansMdot) {
    const double m_dot_com = 0.5;
    const double m_dot_branch = 0.15;
    const double dP0_s = 100.0;
    const double dP0_b = 180.0;
    const double P_static = 1.5e5;
    const double T = 500.0;
    auto Y = dry_air_Y();
    const double theta = M_PI / 3.0;
    const double psi = 1.2;
    const double F_C = 0.008;
    const double eps_m = 1e-5;

    auto res0 = branching_tee_residuals_and_jacobian(
        m_dot_com, m_dot_branch, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);

    auto r_p = branching_tee_residuals_and_jacobian(
        m_dot_com + eps_m, m_dot_branch, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);
    auto r_m = branching_tee_residuals_and_jacobian(
        m_dot_com - eps_m, m_dot_branch, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);
    double fd_Rs = (r_p.R_straight - r_m.R_straight) / (2.0 * eps_m);
    double fd_Rb = (r_p.R_branch - r_m.R_branch) / (2.0 * eps_m);

    EXPECT_NEAR(res0.dR_straight_d_mdot_com, fd_Rs,
                std::abs(fd_Rs) * 1e-4 + 0.5);
    EXPECT_NEAR(res0.dR_branch_d_mdot_com, fd_Rb,
                std::abs(fd_Rb) * 1e-4 + 0.5);

    // Gap 7: FD verification for m_dot_branch derivatives
    auto r_bp = branching_tee_residuals_and_jacobian(
        m_dot_com, m_dot_branch + eps_m, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);
    auto r_bm = branching_tee_residuals_and_jacobian(
        m_dot_com, m_dot_branch - eps_m, dP0_s, dP0_b, P_static, T, Y, theta, psi, F_C);
    double fd_Rs_mb = (r_bp.R_straight - r_bm.R_straight) / (2.0 * eps_m);
    double fd_Rb_mb = (r_bp.R_branch - r_bm.R_branch) / (2.0 * eps_m);

    EXPECT_NEAR(res0.dR_straight_d_mdot_branch, fd_Rs_mb,
                std::abs(fd_Rs_mb) * 1e-4 + 0.5);
    EXPECT_NEAR(res0.dR_branch_d_mdot_branch, fd_Rb_mb,
                std::abs(fd_Rb_mb) * 1e-4 + 0.5);
}

// -----------------------------------------------------------------------------
// Figure 5 validation: contraction coefficient mu (Bassett 2001, Eq. 22)
//
// mu = 1 / (1 + sqrt(1 + 1/(q*psi)^2 - 2*cos(3*theta/4)/(q*psi)))
//
// Figure 5 plots mu vs q for equal-area junctions (psi=1) at five branch
// angles. Reference data digitized from Bassett (2001), Fig. 5.
// Tolerance 0.006 reflects digitization uncertainty (~half a grid division).
//
// Note: K6 = q^2*psi^2*(1 - 1/mu)^2 reduces algebraically to Eq. 27, which
// is identical to the K6 formula in tee_junction.h (verified in test below).
// -----------------------------------------------------------------------------
namespace {
double mu_bassett(double q, double psi, double theta) {
    // Bassett 2001, Eq. 22
    const double c = std::cos(0.75 * theta);
    const double r = 1.0 / (q * psi);
    return 1.0 / (1.0 + std::sqrt(1.0 + r * r - 2.0 * c * r));
}
} // namespace

TEST(TeeJunctionKCoefficients, Figure5ContractionCoefficientMu) {
    // Digitized from Bassett (2001) Fig. 5, equal-area junctions (psi=1.0).
    // Rows: {q, mu_digitized, theta_degrees}
    // theta=0 and theta>90 are outside the validated range but plotted in Fig. 5.
    struct Point { double q, mu; int theta_deg; };
    const Point pts[] = {
        // theta = 0 deg (mu = q analytically for psi=1)
        {0.012, 0.011, 0}, {0.139, 0.137, 0}, {0.296, 0.299, 0},
        {0.491, 0.490, 0}, {0.656, 0.656, 0}, {0.842, 0.843, 0},
        {0.997, 1.000, 0},
        // theta = 30 deg
        {0.023, 0.021, 30}, {0.459, 0.432, 30}, {0.636, 0.569, 30},
        {0.776, 0.654, 30}, {0.939, 0.710, 30}, {0.997, 0.722, 30},
        // theta = 90 deg (within validated range)
        {0.030, 0.025, 90}, {0.227, 0.193, 90}, {0.313, 0.251, 90},
        {0.384, 0.292, 90}, {0.445, 0.320, 90}, {0.529, 0.362, 90},
        {0.598, 0.387, 90}, {0.700, 0.416, 90}, {0.813, 0.446, 90},
        {0.977, 0.469, 90},
        // theta = 120 deg (outside validated range, extrapolation)
        {0.026, 0.024, 120}, {0.170, 0.142, 120}, {0.271, 0.208, 120},
        {0.416, 0.278, 120}, {0.567, 0.332, 120}, {0.710, 0.367, 120},
        {0.846, 0.392, 120}, {0.999, 0.418, 120},
        // theta = 180 deg (outside validated range, extrapolation)
        {0.026, 0.025, 180}, {0.181, 0.137, 180}, {0.333, 0.210, 180},
        {0.470, 0.255, 180}, {0.638, 0.295, 180}, {0.781, 0.320, 180},
        {0.907, 0.340, 180}, {0.998, 0.353, 180},
    };
    const double tol = 0.006; // digitization uncertainty
    for (const auto& p : pts) {
        const double theta = p.theta_deg * M_PI / 180.0;
        const double mu = mu_bassett(p.q, 1.0, theta);
        EXPECT_NEAR(mu, p.mu, tol)
            << "mu mismatch at theta=" << p.theta_deg << " q=" << p.q;
    }
}

// Verify K6 = q^2*psi^2*(1 - 1/mu)^2 matches the closed-form K6 formula
// (Bassett Eq. 26 -> Eq. 27 algebraic identity). Tested at psi=1, 2, 3.
TEST(TeeJunctionKCoefficients, K6EqualsMusquaredFormula) {
    const double qs[] = {0.1, 0.3, 0.5, 0.7, 0.9};
    const double thetas[] = {M_PI / 6.0, M_PI / 4.0, M_PI / 3.0, M_PI / 2.0};
    const double psis[] = {1.0, 2.0, 3.0};
    for (double psi : psis) {
        for (double theta : thetas) {
            for (double q : qs) {
                const double mu = mu_bassett(q, psi, theta);
                const double K6_mu = q * q * psi * psi * (1.0 - 1.0 / mu) * (1.0 - 1.0 / mu);
                EXPECT_NEAR(K6(q, psi, theta), K6_mu, 1e-10)
                    << "K6 != q^2*psi^2*(1-1/mu)^2 at psi=" << psi
                    << " theta=" << theta << " q=" << q;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Bassett Eq. 32: K12(q) - K11(1-q) = 2q-1  for psi=1
// Derived from definition + pA=pB assumption: holds for any theta.
// (Note: K12 and K11 use complementary q values, not the same q.)
// -----------------------------------------------------------------------------
TEST(TeeJunctionKCoefficients, K12MinusK11ComplementIdentity) {
    const double qs[] = {0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9};
    const double thetas[] = {M_PI / 6.0, M_PI / 4.0, M_PI / 3.0, M_PI / 2.0};
    for (double theta : thetas) {
        for (double q : qs) {
            EXPECT_NEAR(K12(q, 1.0, theta) - K11(1.0 - q, 1.0, theta), 2.0 * q - 1.0, 1e-12)
                << "K12(q)-K11(1-q) != 2q-1 at theta=" << theta << " q=" << q;
        }
    }
}

// Same identity for K7/K8 (flow type 4, Bassett Section 4.2).
TEST(TeeJunctionKCoefficients, K7MinusK8ComplementIdentity) {
    const double qs[] = {0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9};
    const double thetas[] = {M_PI / 6.0, M_PI / 4.0, M_PI / 3.0, M_PI / 2.0};
    for (double theta : thetas) {
        for (double q : qs) {
            EXPECT_NEAR(K7(q, 1.0, theta) - K8(1.0 - q, 1.0, theta), 2.0 * q - 1.0, 1e-12)
                << "K7(q)-K8(1-q) != 2q-1 at theta=" << theta << " q=" << q;
        }
    }
}

// -----------------------------------------------------------------------------
// 90-degree symmetry identities (Bassett Table 2 / Section 4.1)
// At theta=pi/2: cos(3*(pi-theta)/4) = cos(3*theta/4) -> K1=K6, K3=K4.
// -----------------------------------------------------------------------------
TEST(TeeJunctionKCoefficients, NinetyDegreeSymmetryK1EqualsK6) {
    const double theta = M_PI / 2.0;
    const double qs[] = {0.1, 0.3, 0.5, 0.7, 0.9};
    const double psis[] = {0.5, 1.0, 2.0};
    for (double psi : psis) {
        for (double q : qs) {
            EXPECT_NEAR(K1(q, psi, theta), K6(q, psi, theta), 1e-14)
                << "K1 != K6 at 90deg, psi=" << psi << " q=" << q;
        }
    }
}

TEST(TeeJunctionKCoefficients, NinetyDegreeSymmetryK3EqualsK4) {
    const double theta = M_PI / 2.0;
    const double qs[] = {0.1, 0.3, 0.5, 0.7, 0.9};
    const double psis[] = {0.5, 1.0, 2.0};
    for (double psi : psis) {
        for (double q : qs) {
            EXPECT_NEAR(K3(q, psi, theta), K4(q, psi, theta), 1e-14)
                << "K3 != K4 at 90deg, psi=" << psi << " q=" << q;
        }
    }
}

// At theta=135deg K4 == K3 at theta=45deg (paper: "135deg curve represents K3 for 45deg").
TEST(TeeJunctionKCoefficients, K4At135EqualK3At45) {
    const double theta45  = M_PI / 4.0;
    const double theta135 = 3.0 * M_PI / 4.0;
    const double qs[] = {0.1, 0.3, 0.5, 0.7, 0.9};
    const double psis[] = {0.5, 1.0, 2.0};
    for (double psi : psis) {
        for (double q : qs) {
            EXPECT_NEAR(K4(q, psi, theta135), K3(q, psi, theta45), 1e-14)
                << "K4(135) != K3(45) at psi=" << psi << " q=" << q;
        }
    }
}

// K9 = K10 for equal-area junction at theta=90deg (Bassett Section 4.2).
TEST(TeeJunctionKCoefficients, K9EqualsK10AtNinetyDegEqualArea) {
    const double theta = M_PI / 2.0;
    const double qs[] = {0.1, 0.3, 0.5, 0.7, 0.9};
    for (double q : qs) {
        EXPECT_NEAR(K9(q, 1.0, theta), K10(q, 1.0, theta), 1e-14)
            << "K9 != K10 at 90deg, psi=1, q=" << q;
    }
}

// -----------------------------------------------------------------------------
// Figure 10c: K12 for branching tee at theta=45deg, various psi.
// Model line data digitized from Bassett (2001) Fig. 10c.
// Tolerances reflect digitization uncertainty (grows with |K12| magnitude).
// -----------------------------------------------------------------------------
TEST(TeeJunctionKCoefficients, Figure10cK12ModelLine) {
    const double theta = M_PI / 4.0;
    struct Point { double q, K12_ref; };

    // psi=1 model line (values in [-0.78, 0.27], tol=0.025)
    const Point psi1[] = {
        {0.0877, -0.76586}, {0.1039, -0.72264}, {0.1172, -0.6856},
        {0.1334, -0.64856}, {0.1496, -0.61152}, {0.1658, -0.57448},
        {0.1820, -0.53743}, {0.1938, -0.48805}, {0.2070, -0.47570},
        {0.2232, -0.44483}, {0.2395, -0.40779}, {0.2557, -0.37692},
        {0.2719, -0.33988}, {0.2881, -0.31518}, {0.3043, -0.28432},
        {0.3205, -0.25345}, {0.3367, -0.22258}, {0.3529, -0.19788},
        {0.3691, -0.17319}, {0.3853, -0.14232}, {0.4015, -0.11763},
        {0.4177, -0.09293}, {0.4339, -0.06824}, {0.4501, -0.04354},
        {0.4663, -0.01885}, {0.4825, -0.00033}, {0.4987,  0.01819},
        {0.5149,  0.03671}, {0.5312,  0.06141}, {0.5474,  0.07376},
        {0.5636,  0.09228}, {0.5798,  0.11080}, {0.5960,  0.12932},
        {0.6122,  0.14167}, {0.6284,  0.15401}, {0.6446,  0.17253},
        {0.6608,  0.18488}, {0.6770,  0.19106}, {0.6932,  0.20340},
        {0.7094,  0.21575}, {0.7256,  0.22810}, {0.7418,  0.23427},
        {0.7580,  0.24044}, {0.7742,  0.24662}, {0.7904,  0.25279},
        {0.8066,  0.25897}, {0.8229,  0.25897}, {0.8391,  0.26514},
        {0.8553,  0.26514}, {0.8715,  0.26514}, {0.8877,  0.27131},
        {0.9039,  0.26514}, {0.9201,  0.26514}, {0.9363,  0.26514},
        {0.9525,  0.25897}, {0.9687,  0.25897}, {0.9849,  0.25279},
        {0.0818, -0.77821}, {0.0656, -0.82142},
    };
    for (const auto& p : psi1)
        EXPECT_NEAR(K12(p.q, 1.0, theta), p.K12_ref, 0.025)
            << "psi=1 model line q=" << p.q;

    // psi=3 model line (values in [-0.43, 2.43], tol=0.10)
    const Point psi3[] = {
        {0.6520,  2.42592}, {0.6402,  2.33331}, {0.6210,  2.22219},
        {0.6137,  2.14810}, {0.5960,  2.04315}, {0.5886,  1.99376},
        {0.5754,  1.90116}, {0.5650,  1.82090}, {0.5444,  1.70977},
        {0.5370,  1.63569}, {0.5135,  1.51222}, {0.4987,  1.40726},
        {0.4958,  1.38257}, {0.5076,  1.45048}, {0.4708,  1.22823},
        {0.4428,  1.08623}, {0.4266,  0.98746}, {0.4133,  0.90720},
        {0.3941,  0.82694}, {0.4030,  0.78373}, {0.3838,  0.74668},
        {0.3765,  0.69112}, {0.3455,  0.54295}, {0.3308,  0.46270},
        {0.3146,  0.38244}, {0.2984,  0.28366}, {0.2792,  0.19723},
        {0.2660,  0.13549}, {0.2542,  0.05524}, {0.2262, -0.06206},
        {0.2188, -0.09911}, {0.1938, -0.17936}, {0.2070, -0.17319},
        {0.1790, -0.27814}, {0.1614, -0.35223}, {0.1437, -0.42631},
    };
    for (const auto& p : psi3)
        EXPECT_NEAR(K12(p.q, 3.0, theta), p.K12_ref, 0.10)
            << "psi=3 model line q=" << p.q;

    // psi=4 model line (values in [-0.69, 2.45], tol=0.09)
    const Point psi4[] = {
        {0.4619,  2.45061}, {0.4516,  2.33949}, {0.4383,  2.18514},
        {0.4310,  2.08019}, {0.4148,  1.91968}, {0.3971,  1.75916},
        {0.3956,  1.67890}, {0.3853,  1.64186}, {0.3765,  1.54308},
        {0.3617,  1.38874}, {0.3426,  1.23440}, {0.3264,  1.08623},
        {0.3175,  0.99363}, {0.2969,  0.81459}, {0.2851,  0.72816},
        {0.2763,  0.62321}, {0.2542,  0.47504}, {0.2453,  0.38861},
        {0.2262,  0.23427}, {0.1967,  0.11080}, {0.2011,  0.06141},
        {0.1953, -0.00650}, {0.1953,  0.01819}, {0.1732, -0.11763},
        {0.1525, -0.24110}, {0.1407, -0.33370}, {0.1157, -0.46952},
        {0.0965, -0.51274}, {0.1098, -0.49422}, {0.0907, -0.59917},
        {0.0833, -0.68560},
    };
    for (const auto& p : psi4)
        EXPECT_NEAR(K12(p.q, 4.0, theta), p.K12_ref, 0.09)
            << "psi=4 model line q=" << p.q;
}

// Bassett model accuracy at psi=2 (intermediate area ratio, theta=45 deg).
// The Bassett correlation systematically underpredicts K12 at high q for psi=2:
// errors grow from ~0.10 at q=0 to ~0.52 at q=1. This is the primary known
// limitation and the target for the Stigler vector-momentum upgrade.
// Tolerance 0.55 establishes the current Bassett baseline for future comparison.
TEST(TeeJunctionKCoefficients, Figure10cK12BassettLimitationPsi2) {
    const double theta = M_PI / 4.0;
    struct Point { double q, K12_exp; };
    const Point psi2[] = {
        {-0.0007, -0.8996}, {0.1988, -0.1758}, {0.2994,  0.0994},
        { 0.3999,  0.4249}, {0.5008,  0.6684}, {0.6000,  1.0030},
        { 0.6993,  1.4040}, {0.8006,  1.6315}, {0.8993,  2.1635},
        { 0.9975,  2.4141},
    };
    for (const auto& p : psi2)
        EXPECT_NEAR(K12(p.q, 2.0, theta), p.K12_exp, 0.55)
            << "Bassett psi=2 baseline q=" << p.q;
}

// Experimental scatter: physical measurements from Bassett (2001) Fig. 10c.
// psi=1, 3, 4 agree with the Bassett model within experimental scatter (tol=0.20).
TEST(TeeJunctionKCoefficients, Figure10cK12ExperimentalScatter) {
    const double theta = M_PI / 4.0;
    struct Point { double q, K12_exp; };

    const Point psi1[] = {
        {0.0996, -0.7180}, {0.2000, -0.4772}, {0.4005, -0.0212},
        {0.5009,  0.0894}, {0.6001,  0.2264}, {0.7000,  0.2745},
        {0.8011,  0.3166}, {0.9003,  0.3137}, {0.9994,  0.3018},
    };
    for (const auto& p : psi1)
        EXPECT_NEAR(K12(p.q, 1.0, theta), p.K12_exp, 0.20)
            << "psi=1 exp q=" << p.q;

    const Point psi3[] = {
        {0.0992, -0.5024}, {0.1989, -0.0052}, {0.2981,  0.3757},
        {0.3997,  0.8005}, {0.4998,  1.4596}, {0.5993,  2.0562},
    };
    for (const auto& p : psi3)
        EXPECT_NEAR(K12(p.q, 3.0, theta), p.K12_exp, 0.20)
            << "psi=3 exp q=" << p.q;

    const Point psi4[] = {
        {0.1002, -0.4012}, {0.2000,  0.0973},
        {0.2986,  0.8468}, {0.3993,  1.6536},
    };
    for (const auto& p : psi4)
        EXPECT_NEAR(K12(p.q, 4.0, theta), p.K12_exp, 0.20)
            << "psi=4 exp q=" << p.q;
}

// -----------------------------------------------------------------------------
// Unified0D compressible tee: Jacobian verification via central FD
// Convention matches Python caller: m_str = m_com - m_branch.
// Perturbing m_dot_com also adjusts str.m_dot by the same delta.
// Perturbing m_dot_branch also adjusts str.m_dot by -delta.
// All other Jacobian entries perturb only their own BranchInput field.
// -----------------------------------------------------------------------------
namespace {
BranchInput make_bi(double P_s, double Pt, double T, double m_dot,
                     double A, double theta, double gamma = 1.4,
                     double R_gas = 287.05) {
    BranchInput b{};
    b.P_static = P_s; b.Pt = Pt; b.T = T; b.m_dot = m_dot;
    b.A = A; b.theta = theta; b.gamma_eff = gamma; b.R_gas = R_gas;
    return b;
}
} // namespace

TEST(CompressibleTeeJacobian, BranchingFDAllFields) {
    const double eps = 1e-5;
    const double tol_rel = 2e-4;  // central FD tolerance
    const double tol_abs = 1.0;   // Pa absolute floor

    // Typical branching state (M_dat ~ 0.25 in common duct)
    BranchInput com = make_bi(1.95e5, 2.0e5, 400.0,  0.50, 0.01, 0.0);
    BranchInput str = make_bi(1.88e5, 1.93e5, 402.0, 0.30, 0.01, 0.0);
    BranchInput bra = make_bi(1.82e5, 1.90e5, 398.0, 0.20, 0.01, M_PI / 2.0);

    CompressibleTeeResult res0 = compressible_branching_tee_rj(com, str, bra);

    auto check = [&](double anal, double fd_val, const char* name) {
        double scale = std::max(std::abs(fd_val), std::abs(anal)) + tol_abs;
        EXPECT_NEAR(anal, fd_val, scale * tol_rel) << name;
    };

    // dR/dPt_com
    { auto cp = com; cp.Pt += eps; auto cm = com; cm.Pt -= eps;
      check(res0.dR0_dPt_com,
            (compressible_branching_tee_rj(cp,str,bra).R_0 - compressible_branching_tee_rj(cm,str,bra).R_0)/(2*eps),
            "dR0_dPt_com");
      check(res0.dR1_dPt_com,
            (compressible_branching_tee_rj(cp,str,bra).R_1 - compressible_branching_tee_rj(cm,str,bra).R_1)/(2*eps),
            "dR1_dPt_com"); }

    // dR/dP_static_com
    { auto cp = com; cp.P_static += eps; auto cm = com; cm.P_static -= eps;
      check(res0.dR0_dP_com,
            (compressible_branching_tee_rj(cp,str,bra).R_0 - compressible_branching_tee_rj(cm,str,bra).R_0)/(2*eps),
            "dR0_dP_com");
      check(res0.dR1_dP_com,
            (compressible_branching_tee_rj(cp,str,bra).R_1 - compressible_branching_tee_rj(cm,str,bra).R_1)/(2*eps),
            "dR1_dP_com"); }

    // dR/dT_com
    { auto cp = com; cp.T += eps; auto cm = com; cm.T -= eps;
      check(res0.dR0_dT_com,
            (compressible_branching_tee_rj(cp,str,bra).R_0 - compressible_branching_tee_rj(cm,str,bra).R_0)/(2*eps),
            "dR0_dT_com");
      check(res0.dR1_dT_com,
            (compressible_branching_tee_rj(cp,str,bra).R_1 - compressible_branching_tee_rj(cm,str,bra).R_1)/(2*eps),
            "dR1_dT_com"); }

    // dR/dPt_str
    { auto sp = str; sp.Pt += eps; auto sm = str; sm.Pt -= eps;
      check(res0.dR0_dPt_str,
            (compressible_branching_tee_rj(com,sp,bra).R_0 - compressible_branching_tee_rj(com,sm,bra).R_0)/(2*eps),
            "dR0_dPt_str"); }

    // dR/dPt_bra
    { auto bp = bra; bp.Pt += eps; auto bm = bra; bm.Pt -= eps;
      check(res0.dR1_dPt_bra,
            (compressible_branching_tee_rj(com,str,bp).R_1 - compressible_branching_tee_rj(com,str,bm).R_1)/(2*eps),
            "dR1_dPt_bra"); }

    // dR/dm_dot_com: caller convention m_str = m_com - m_branch, so str.m_dot changes too
    { auto cp = com; cp.m_dot += eps; auto sp = str; sp.m_dot += eps;
      auto cm = com; cm.m_dot -= eps; auto sm = str; sm.m_dot -= eps;
      check(res0.dR0_dmdot_com,
            (compressible_branching_tee_rj(cp,sp,bra).R_0 - compressible_branching_tee_rj(cm,sm,bra).R_0)/(2*eps),
            "dR0_dmdot_com");
      check(res0.dR1_dmdot_com,
            (compressible_branching_tee_rj(cp,sp,bra).R_1 - compressible_branching_tee_rj(cm,sm,bra).R_1)/(2*eps),
            "dR1_dmdot_com"); }

    // dR/dm_dot_branch: bra.m_dot up, str.m_dot down (dm_str/dm_branch = -1)
    { auto bp = bra; bp.m_dot += eps; auto sp = str; sp.m_dot -= eps;
      auto bm = bra; bm.m_dot -= eps; auto sm = str; sm.m_dot += eps;
      check(res0.dR0_dmdot_branch,
            (compressible_branching_tee_rj(com,sp,bp).R_0 - compressible_branching_tee_rj(com,sm,bm).R_0)/(2*eps),
            "dR0_dmdot_branch");
      check(res0.dR1_dmdot_branch,
            (compressible_branching_tee_rj(com,sp,bp).R_1 - compressible_branching_tee_rj(com,sm,bm).R_1)/(2*eps),
            "dR1_dmdot_branch"); }
}

TEST(CompressibleTeeJacobian, MergingFDAllFields) {
    const double eps = 1e-5;
    const double tol_rel = 2e-4;
    const double tol_abs = 1.0;

    // Typical merging state (M_dat ~ 0.25 in combined supplier stream)
    BranchInput com = make_bi(1.88e5, 1.95e5, 401.0, -0.50, 0.01, 0.0);
    BranchInput str = make_bi(1.95e5, 2.0e5,  400.0,  0.30, 0.01, 0.0);
    BranchInput bra = make_bi(1.93e5, 1.98e5, 402.0,  0.20, 0.01, M_PI / 2.0);

    CompressibleTeeResult res0 = compressible_merging_tee_rj(com, str, bra);

    auto check = [&](double anal, double fd_val, const char* name) {
        double scale = std::max(std::abs(fd_val), std::abs(anal)) + tol_abs;
        EXPECT_NEAR(anal, fd_val, scale * tol_rel) << name;
    };

    // dR/dPt_com
    { auto cp = com; cp.Pt += eps; auto cm = com; cm.Pt -= eps;
      check(res0.dR0_dPt_com,
            (compressible_merging_tee_rj(cp,str,bra).R_0 - compressible_merging_tee_rj(cm,str,bra).R_0)/(2*eps),
            "dR0_dPt_com");
      check(res0.dR1_dPt_com,
            (compressible_merging_tee_rj(cp,str,bra).R_1 - compressible_merging_tee_rj(cm,str,bra).R_1)/(2*eps),
            "dR1_dPt_com"); }

    // dR/dPt_str
    { auto sp = str; sp.Pt += eps; auto sm = str; sm.Pt -= eps;
      check(res0.dR0_dPt_str,
            (compressible_merging_tee_rj(com,sp,bra).R_0 - compressible_merging_tee_rj(com,sm,bra).R_0)/(2*eps),
            "dR0_dPt_str"); }

    // dR/dPt_bra
    { auto bp = bra; bp.Pt += eps; auto bm = bra; bm.Pt -= eps;
      check(res0.dR1_dPt_bra,
            (compressible_merging_tee_rj(com,str,bp).R_1 - compressible_merging_tee_rj(com,str,bm).R_1)/(2*eps),
            "dR1_dPt_bra"); }

    // dR/dP_static_str
    { auto sp = str; sp.P_static += eps; auto sm = str; sm.P_static -= eps;
      check(res0.dR0_dP_str,
            (compressible_merging_tee_rj(com,sp,bra).R_0 - compressible_merging_tee_rj(com,sm,bra).R_0)/(2*eps),
            "dR0_dP_str");
      check(res0.dR1_dP_str,
            (compressible_merging_tee_rj(com,sp,bra).R_1 - compressible_merging_tee_rj(com,sm,bra).R_1)/(2*eps),
            "dR1_dP_str"); }

    // dR/dT_str
    { auto sp = str; sp.T += eps; auto sm = str; sm.T -= eps;
      check(res0.dR0_dT_str,
            (compressible_merging_tee_rj(com,sp,bra).R_0 - compressible_merging_tee_rj(com,sm,bra).R_0)/(2*eps),
            "dR0_dT_str");
      check(res0.dR1_dT_str,
            (compressible_merging_tee_rj(com,sp,bra).R_1 - compressible_merging_tee_rj(com,sm,bra).R_1)/(2*eps),
            "dR1_dT_str"); }

    // dR/dP_static_bra
    { auto bp = bra; bp.P_static += eps; auto bm = bra; bm.P_static -= eps;
      check(res0.dR0_dP_bra,
            (compressible_merging_tee_rj(com,str,bp).R_0 - compressible_merging_tee_rj(com,str,bm).R_0)/(2*eps),
            "dR0_dP_bra");
      check(res0.dR1_dP_bra,
            (compressible_merging_tee_rj(com,str,bp).R_1 - compressible_merging_tee_rj(com,str,bm).R_1)/(2*eps),
            "dR1_dP_bra"); }

    // dR/dT_bra
    { auto bp = bra; bp.T += eps; auto bm = bra; bm.T -= eps;
      check(res0.dR0_dT_bra,
            (compressible_merging_tee_rj(com,str,bp).R_0 - compressible_merging_tee_rj(com,str,bm).R_0)/(2*eps),
            "dR0_dT_bra");
      check(res0.dR1_dT_bra,
            (compressible_merging_tee_rj(com,str,bp).R_1 - compressible_merging_tee_rj(com,str,bm).R_1)/(2*eps),
            "dR1_dT_bra"); }

    // dR/dm_dot_com: m_str = m_com - m_branch, dm_str/dm_com = 1
    { auto cp = com; cp.m_dot += eps; auto sp = str; sp.m_dot += eps;
      auto cm = com; cm.m_dot -= eps; auto sm = str; sm.m_dot -= eps;
      check(res0.dR0_dmdot_com,
            (compressible_merging_tee_rj(cp,sp,bra).R_0 - compressible_merging_tee_rj(cm,sm,bra).R_0)/(2*eps),
            "dR0_dmdot_com");
      check(res0.dR1_dmdot_com,
            (compressible_merging_tee_rj(cp,sp,bra).R_1 - compressible_merging_tee_rj(cm,sm,bra).R_1)/(2*eps),
            "dR1_dmdot_com"); }

    // dR/dm_dot_branch: bra.m_dot up, str.m_dot down (dm_str/dm_branch = -1)
    { auto bp = bra; bp.m_dot += eps; auto sp = str; sp.m_dot -= eps;
      auto bm = bra; bm.m_dot -= eps; auto sm = str; sm.m_dot += eps;
      check(res0.dR0_dmdot_branch,
            (compressible_merging_tee_rj(com,sp,bp).R_0 - compressible_merging_tee_rj(com,sm,bm).R_0)/(2*eps),
            "dR0_dmdot_branch");
      check(res0.dR1_dmdot_branch,
            (compressible_merging_tee_rj(com,sp,bp).R_1 - compressible_merging_tee_rj(com,sm,bm).R_1)/(2*eps),
            "dR1_dmdot_branch"); }
}

TEST(CompressibleTeeJacobian, ZeroFlowIsRegularised) {
    BranchInput zero = make_bi(1e5, 1e5, 400.0, 1e-15, 0.01, 0.0);
    auto res_b = compressible_branching_tee_rj(zero, zero, zero);
    auto res_m = compressible_merging_tee_rj(zero, zero, zero);
    for (double v : {res_b.R_0, res_b.R_1, res_b.dR0_dmdot_com, res_b.dR1_dmdot_com,
                     res_m.R_0, res_m.R_1, res_m.dR0_dmdot_com, res_m.dR1_dmdot_com})
        EXPECT_FALSE(std::isnan(v)) << "NaN in zero-flow result";
}

// Zero flow: regularisation must prevent NaN in all outputs
TEST(TeeJunctionSolverInterface, ZeroFlowIsRegularised) {
    auto Y = dry_air_Y();
    auto res = merging_tee_residuals_and_jacobian(
        0.0, 0.0, 50.0, 50.0, 1e5, 400.0, Y, M_PI / 2.0, 1.0, 0.01);
    EXPECT_FALSE(std::isnan(res.R_straight));
    EXPECT_FALSE(std::isnan(res.R_branch));
    EXPECT_FALSE(std::isnan(res.dR_straight_d_mdot_com));
    EXPECT_FALSE(std::isnan(res.dR_branch_d_mdot_com));
    EXPECT_FALSE(std::isnan(res.dR_straight_d_mdot_branch));
    EXPECT_FALSE(std::isnan(res.dR_branch_d_mdot_branch));
}
