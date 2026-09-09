#pragma once

// -----------------------------------------------------------------------------
// Mynard & Valen-Sendstad Unified0D junction loss closure, on dual numbers.
//
// C++ port of python/combaero/network/_mynard2010.py::junction_loss_coefficient
// -- see that module for the derivation, the paper equation numbers and the
// provenance of every constant. This header carries only what the port itself
// needs to state; it deliberately does not restate the physics.
//
// Kept in exact correspondence with the Python so the two can be compared
// row-for-row against validation/junction/data/mynard_golden_data.h (see
// tests/test_mynard_junction.cpp).
//
// Reference
// ---------
// [Mynard & Valen-Sendstad 2015] J.P. Mynard, K. Valen-Sendstad,
//     "A unified method for estimating pressure losses at vascular
//     junctions", International Journal for Numerical Methods in
//     Biomedical Engineering 31(7) (2015) e02717.
//
// Scope: three ports. The Python computes C for general N and K only for
// N <= 3; MultiPortChamberElement refuses N > 3 at construction, so nothing that ships
// exercises the general-N path and this port does not carry it.
//
// BRANCH-ON-PRIMAL. Six decisions are taken on primal values rather than on
// duals, each marked at its site below:
//
//   1. supplier vs collector, from sign(Q)
//   2. the second-quadrant flip of the pseudocollector angle
//   3. the pseudosupplier direction flip, from sign of a Q-weighted moment
//   4. which multiple of 2*pi the two angle wraps remove
//   5. which port is "common" for the K normalisation
//   6. whether a collector continues straight (the dividing-streamline term)
//
// Within a branch the derivative is exact. At a switching surface it is
// one-sided, which is the correct thing to hand Newton and is far better than
// a zero column. Note that (1) is a genuine seam in the residual, not only in
// its derivative: crossing it changes which formula is evaluated. That is a
// property of the model, is present identically in the Python, and is
// discussed in validation/junction/MPCE_CPP_PORT_DESIGN.md.
// -----------------------------------------------------------------------------

#include <array>
#include <cmath>
#include <cstddef>

#include "dual_number.h"
#include "math_constants.h"

namespace combaero::solver {

// Mynard Eq 36's CFD-fitted energy-transfer coefficients. Gated by eta_scale,
// whose production default is 0.0 -- see MultiPortChamberElement.DEFAULT_ETA_SCALE for
// why the term is off and what measurement retired it.
constexpr double kMynardEtaA0 = 0.8;
constexpr double kMynardEtaA1 = -0.2;

// Regulariser on the collector loss coefficient, (1 - exp(-FlowRatio/tau)).
// From the reference Matlab, not from the paper: it keeps C finite as a
// collector's flow ratio goes to zero.
constexpr double kFlowRatioDamping = 0.02;

// Dividing-streamline pressure recovery on a collector that continues
// straight, and the angle within which "continues straight" is judged.
constexpr double kDividingStreamlineRecovery = 0.5;
constexpr double kCollinearTolRad = 1.0e-6;

// A port is excluded from both masks only at exactly zero; the Python uses
// Q < 0.0 for the collector test, so a zero-flow port counts as a supplier.
// Reproduced exactly -- the element snaps such ports before calling in, and
// the two classifications have to agree.
constexpr int kMynardPorts = 3;

template <int M> struct MynardResult {
  // Per branch; supplier entries are zero, as in the Python.
  std::array<DualN<M>, kMynardPorts> C{};
  // Per NON-COMMON port, in ascending port order. Which ports those are
  // depends on the flow pattern: the collectors when there is one supplier,
  // the suppliers when there is one collector.
  std::array<DualN<M>, kMynardPorts - 1> K{};
  std::array<int, kMynardPorts - 1> k_port{};   // port index each K belongs to
  std::array<DualN<M>, kMynardPorts> flow_ratio{}; // per collector; else zero
  int n_supplier = 0;
  int n_collector = 0;
  bool valid = false; // false when the split has no supplier or no collector
};

// U: per-port velocity, POSITIVE INTO the junction (supplier), negative out.
// A: per-port area. theta: per-port angle from an arbitrary reference, radians.
//
// Returns valid=false rather than throwing when the flow pattern has no
// supplier or no collector; the caller owns that decision, exactly as the
// Python's caller does.
template <int M>
MynardResult<M> mynard_junction_loss_coefficient(
    const std::array<DualN<M>, kMynardPorts>& U,
    const std::array<double, kMynardPorts>& A,
    const std::array<double, kMynardPorts>& theta_in,
    double joining_etransfer_alpha = 0.0, double eta_scale = 0.0) {
  using D = DualN<M>;
  MynardResult<M> out;

  // Angles are primal throughout. They come from the declared geometry and
  // from means over it; only WHICH angles are selected depends on the flow,
  // and that is decision (1) below.
  std::array<double, kMynardPorts> theta{};
  for (int i = 0; i < kMynardPorts; ++i) {
    double t = std::fmod(theta_in[i] + M_PI, 2.0 * M_PI);
    if (t < 0.0) t += 2.0 * M_PI;
    theta[i] = t - M_PI;
  }

  std::array<D, kMynardPorts> Q{};
  for (int i = 0; i < kMynardPorts; ++i) Q[i] = U[i] * A[i];

  // (1) BRANCH-ON-PRIMAL: supplier/collector from the sign of Q.
  std::array<bool, kMynardPorts> is_collector{};
  for (int i = 0; i < kMynardPorts; ++i) {
    is_collector[i] = Q[i].v < 0.0;
    if (is_collector[i]) {
      ++out.n_collector;
    } else {
      ++out.n_supplier;
    }
  }
  if (out.n_supplier == 0 || out.n_collector == 0) return out;

  D qtot = D::constant(0.0);
  for (int i = 0; i < kMynardPorts; ++i) {
    if (!is_collector[i]) qtot = qtot + Q[i];
  }
  for (int i = 0; i < kMynardPorts; ++i) {
    if (is_collector[i]) out.flow_ratio[i] = (0.0 - Q[i]) / qtot;
  }

  // Reorient so the pseudocollector sits at zero.
  double pseudo_col_angle = 0.0;
  for (int i = 0; i < kMynardPorts; ++i) {
    if (is_collector[i]) pseudo_col_angle += theta[i];
  }
  pseudo_col_angle /= static_cast<double>(out.n_collector);

  D sin_moment = D::constant(0.0);
  D cos_moment = D::constant(0.0);
  for (int i = 0; i < kMynardPorts; ++i) {
    if (is_collector[i]) continue;
    sin_moment = sin_moment + Q[i] * std::sin(theta[i]);
    cos_moment = cos_moment + Q[i] * std::cos(theta[i]);
  }
  D pseudo_sup_initial = datan2(sin_moment, cos_moment);

  // (2) BRANCH-ON-PRIMAL: put the pseudosupplier in the second quadrant
  // relative to the collector.
  if (std::abs(pseudo_sup_initial.v - pseudo_col_angle) < M_PI / 2.0) {
    pseudo_col_angle += M_PI;
  }
  for (int i = 0; i < kMynardPorts; ++i) {
    // (4) BRANCH-ON-PRIMAL: the wrap's multiple of 2*pi. Primal-only here.
    double t = std::fmod(theta[i] - pseudo_col_angle + M_PI, 2.0 * M_PI);
    if (t < 0.0) t += 2.0 * M_PI;
    theta[i] = t - M_PI;
  }

  // (3) BRANCH-ON-PRIMAL: the pseudosupplier direction. The moment is a dual,
  // but only its SIGN is used, so the flip itself contributes no derivative.
  D dir_moment = D::constant(0.0);
  for (int i = 0; i < kMynardPorts; ++i) {
    if (!is_collector[i]) dir_moment = dir_moment + Q[i] * std::sin(theta[i]);
  }
  // np.sign of a mean; the division by the count cannot change the sign.
  if (dir_moment.v < 0.0) {
    for (int i = 0; i < kMynardPorts; ++i) theta[i] = -theta[i];
  }

  D abs_sin_moment = D::constant(0.0);
  D abs_cos_moment = D::constant(0.0);
  for (int i = 0; i < kMynardPorts; ++i) {
    if (is_collector[i]) continue;
    abs_sin_moment = abs_sin_moment + Q[i] * std::sin(std::abs(theta[i]));
    abs_cos_moment = abs_cos_moment + Q[i] * std::cos(std::abs(theta[i]));
  }
  D pseudo_sup_angle = datan2(abs_sin_moment, abs_cos_moment);

  // Mynard's empirical energy transfer, per collector. eta_scale = 0 zeroes
  // it without changing the arithmetic, matching the Python.
  std::array<D, kMynardPorts> etransfer{};
  for (int i = 0; i < kMynardPorts; ++i) {
    if (!is_collector[i]) continue;
    double sgn = (theta[i] < 0.0) ? -1.0 : ((theta[i] > 0.0) ? 1.0 : 0.0);
    D eta = ((M_PI - pseudo_sup_angle) * (kMynardEtaA0 * sgn) + kMynardEtaA1) * eta_scale;
    etransfer[i] = eta * (1.0 - out.flow_ratio[i]);
  }

  // combaero's joining-side asymmetry correction (not in the paper). Vanishes
  // at equal supplier areas; provenance on MultiPortChamberElement.
  if (joining_etransfer_alpha != 0.0 && out.n_supplier >= 2) {
    double a_max = 0.0;
    double a_min = 0.0;
    bool first = true;
    for (int i = 0; i < kMynardPorts; ++i) {
      if (is_collector[i]) continue;
      if (first) {
        a_max = A[i];
        a_min = A[i];
        first = false;
      } else {
        a_max = std::max(a_max, A[i]);
        a_min = std::min(a_min, A[i]);
      }
    }
    if (a_max + a_min > 0.0) {
      double asym = joining_etransfer_alpha * (a_max - a_min) / (a_max + a_min);
      for (int i = 0; i < kMynardPorts; ++i) {
        if (is_collector[i]) etransfer[i] = etransfer[i] + asym;
      }
    }
  }

  D uq = D::constant(0.0);
  for (int i = 0; i < kMynardPorts; ++i) {
    if (!is_collector[i]) uq = uq + U[i] * Q[i];
  }
  D pseudo_velocity_avg = uq / qtot;

  for (int i = 0; i < kMynardPorts; ++i) {
    if (!is_collector[i]) continue;
    D tot_pseudo_area = qtot / ((1.0 - etransfer[i]) * pseudo_velocity_avg);
    D area_ratio = tot_pseudo_area / A[i];
    // (4) BRANCH-ON-PRIMAL again, this time on a dual: the wrap passes the
    // partials through and shifts only the value.
    D phi = dwrap_to_2pi(pseudo_sup_angle - theta[i]);
    D damping = 1.0 - dexp((0.0 - out.flow_ratio[i]) / kFlowRatioDamping);
    out.C[i] = damping * (1.0 - (1.0 / (area_ratio * out.flow_ratio[i])) *
                                    dcos((M_PI - phi) * 0.75));
  }

  // (5) BRANCH-ON-PRIMAL: the common port. One collector means the collector
  // normalises; otherwise the (single) supplier does. The K entries are then
  // indexed by the OTHER side, which is why K is per non-common port.
  const bool collector_is_common = (out.n_collector == 1);
  int common = -1;
  for (int i = 0; i < kMynardPorts; ++i) {
    if (is_collector[i] == collector_is_common) {
      common = i;
      break;
    }
  }
  D u_com = U[common];

  int k = 0;
  for (int i = 0; i < kMynardPorts; ++i) {
    if (is_collector[i] == collector_is_common) continue;
    // The Python broadcasts U[Ci]**2 / Ucom**2 against U[Si]**2 / U[Ci]**2,
    // so the collector term is the common one when there are two suppliers
    // and the per-entry one when there are two collectors. Written out here
    // rather than relying on a broadcast the reader has to reconstruct.
    const int col = collector_is_common ? common : i;
    const int sup = collector_is_common ? i : common;
    D ratio_col = (U[col] * U[col]) / (u_com * u_com);
    D ratio_sup = (U[sup] * U[sup]) / (U[col] * U[col]);
    out.K[k] = ratio_col * (out.C[col] * 2.0 + ratio_sup - 1.0);
    out.k_port[k] = i;
    ++k;
  }

  // Dividing-streamline recovery, single-supplier case only. The Python's
  // `phi` here is the collector's, recomputed for the test.
  if (out.n_supplier == 1) {
    k = 0;
    for (int i = 0; i < kMynardPorts; ++i) {
      if (is_collector[i] == collector_is_common) continue;
      D phi = dwrap_to_2pi(pseudo_sup_angle - theta[i]);
      // (6) BRANCH-ON-PRIMAL: with one supplier the pseudosupplier angle is
      // that supplier's own, so this test is on geometry and introduces no
      // discontinuity in the residual -- see the Python's note.
      if (std::abs(std::abs(phi.v) - M_PI) < kCollinearTolRad) {
        out.K[k] = out.K[k] - (1.0 - out.flow_ratio[i]) * kDividingStreamlineRecovery;
      }
      ++k;
    }
  }

  out.valid = true;
  return out;
}

} // namespace combaero::solver
