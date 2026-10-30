#pragma once

// -----------------------------------------------------------------------------
// Jobson (1955) orifice closure, one templated body (#854).
//
// The network kernel (solver_interface.cpp, DualN forward mode for its
// Jacobian) and the plain-double API in orifice.h (Table 2, expansibility)
// evaluate the same function, so the reproduction tests exercise the code the
// solver runs. See orifice.h for the physics and the beta extension.
//
// Scalar inputs are DualN<N>: p0, rho0 the upstream stagnation state; p1 the
// approach STATIC pressure (p1 = p0 for a reservoir), whose density follows
// the isentrope through (p0, rho0); p the back pressure, p < p1; Ci the
// incompressible contraction; n the isentropic index. All in SI.
// -----------------------------------------------------------------------------

#include <cmath>

#include "dual_number.h"

namespace combaero::orifice_jobson {

using solver::DualN;

// x^y with both arguments carrying derivatives: exp(y ln x), x > 0.
template <int N> DualN<N> pow_dd(const DualN<N>& x, const DualN<N>& y) {
  return solver::dexp(y * solver::dlog(x));
}

// Ci from the element's incompressible coefficient in ISO's form
// (m = Cd E A sqrt(2 rho1 dp)) and the velocity coefficient Cv. The actual
// flow is Cv times the theoretical one and the approach velocity is the
// actual flow's, so the incompressible limit reads
// Cv Ci / sqrt(1 - (Cv Ci)^2 beta^4) = Cd E, E = 1/sqrt(1 - beta^4).
template <int N> DualN<N> ci_from_cd(const DualN<N>& Cd, double beta, double Cv) {
  if (!(beta > 0.0 && beta < 1.0)) {
    return Cd / Cv;
  }
  const double b4 = beta * beta * beta * beta;
  const double E = 1.0 / std::sqrt(1.0 - b4);
  const DualN<N> k = Cd * E;
  return k / solver::dsqrt(1.0 + k * k * b4) / Cv;
}

template <int N> struct Flux {
  DualN<N> g;   // m / (Cv A): the jet's mass flow per bore area [kg/(s m^2)]
  DualN<N> C;   // contraction coefficient, jet area / bore area [-]
  DualN<N> G;   // isentropic mass flux of the jet section, K_N sqrt(p0 rho0)
  bool choked = false;
  // The momentum quadratic's discriminant fell to zero (held there): the
  // closure has no real root, only reachable far outside the paper's Ci.
  bool degenerate = false;
};

// The momentum balance, solved for the mass flux g = C G:
//   subcritical  (p1 - p) + f_eff g^2 / rho1 = g u                 (Eq. 11)
//   choked       (p1 - p) - (pc - p) C + f_eff g^2 / rho1 = g uc   (Eq. 15)
// with f_eff = f(Ci, beta, Cv) + Cv beta^2 (the approach momentum g V1,
// V1 = Cv g beta^2 / rho1 the ACTUAL approach velocity that sets the pipe's
// static p1, moves to the left), so f_eff = 1/Ci - 1/(2 Ci^2) + Cv^2
// beta^4 / 2 and both read f_eff g^2/rho1 - B g + (p1 - p) = 0 with
// B = u, or uc + (pc - p) / Gc choked. B, G and therefore g are C1 through
// the critical ratio: dG/dp = 0 there and dB/dp = -1/Gc on both sides.
// Jobson's root (the minus sign), written 2c / (B + sqrt(B^2 - 4 f_eff c /
// rho1)) so it stays finite for f_eff <= 0 (Ci <= 0.5). Branch on the primal
// pressure ratio; one-sided derivatives exactly at r_c, where they agree.
template <int N>
Flux<N> flux(const DualN<N>& p0, const DualN<N>& rho0, const DualN<N>& p1, const DualN<N>& p,
             const DualN<N>& Ci, const DualN<N>& n, double beta, double Cv = 1.0) {
  const double b2 = (beta > 0.0 && beta < 1.0) ? beta * beta : 0.0;
  const double b4 = b2 * b2;
  // f(Ci, beta, Cv) + Cv beta^2 = 1/Ci - 1/(2 Ci^2) + Cv^2 beta^4 / 2.
  const DualN<N> f_eff = 1.0 / Ci - 0.5 / (Ci * Ci) + 0.5 * Cv * Cv * b4;
  const DualN<N> rho1 = rho0 * pow_dd(p1 / p0, 1.0 / n);
  const DualN<N> rc = pow_dd(2.0 / (n + 1.0), n / (n - 1.0));
  const DualN<N> r = p / p0;

  Flux<N> out;
  out.choked = r.v < rc.v;
  const DualN<N> re = out.choked ? rc : r;
  // Jet velocity and density on the isentrope from the TOTAL state to re.
  const DualN<N> u =
      solver::dsqrt(2.0 * n / (n - 1.0) * (p0 / rho0) * (1.0 - pow_dd(re, (n - 1.0) / n)));
  out.G = rho0 * pow_dd(re, 1.0 / n) * u;
  const DualN<N> B = out.choked ? u + (re * p0 - p) / out.G : u;
  const DualN<N> c = p1 - p;
  const DualN<N> disc = B * B - 4.0 * f_eff * c / rho1;
  out.degenerate = !(disc.v > 0.0);
  const DualN<N> root = out.degenerate ? DualN<N>::constant(0.0) : solver::dsqrt(disc);
  out.g = 2.0 * c / (B + root);
  out.C = out.g / out.G;
  return out;
}

}  // namespace combaero::orifice_jobson
