#include "../include/fanno_mach.h"

#include "../include/solver_interface.h"
#include "../include/stagnation.h"
#include "../include/thermo.h"
#include "../include/transport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace combaero {

namespace {

double sound_speed_sq(double T, const std::vector<double>& X) {
  return isentropic_expansion_coefficient(T, X) * specific_gas_constant(X) * T;
}

// d(a^2)/dT analytically: a^2 = gamma R T, gamma = cp/cv with cv = cp - R_u
// (molar), so dgamma/dT = -R_u cp' / cv^2. Exact, so the length integrand
// carries no finite-difference noise and the quadrature can meet its
// tolerance.
double sound_speed_sq_dT(double T, const std::vector<double>& X) {
  const double cp_mol = cp(T, X);
  const double cv_mol = cp_mol - thermo::R_GAS;
  const double gamma = cp_mol / cv_mol;
  const double dgamma = -thermo::R_GAS * dcp_dT(T, X) / (cv_mol * cv_mol);
  return specific_gas_constant(X) * (gamma + T * dgamma);
}

// Static temperature at Mach M from the stagnation temperature:
// h(T) + M^2 a(T)^2 / 2 = h(Tt). F(T) rises with T, F(Tt) >= 0, so a
// safeguarded Newton on [T_lo, Tt] cannot leave the physical branch.
double static_T(double Tt, double M, const std::vector<double>& X) {
  const double h0 = h_mass(Tt, X);
  const double R = specific_gas_constant(X);
  auto F = [&](double T) {
    const double a2 = isentropic_expansion_coefficient(T, X) * R * T;
    return h_mass(T, X) + 0.5 * M * M * a2 - h0;
  };
  double lo = std::max(20.0, 0.05 * Tt);
  double hi = Tt;
  if (M <= 0.0) {
    return Tt;
  }
  // Ideal-gas guess at constant gamma, then Newton with the bracket.
  const double g = isentropic_expansion_coefficient(Tt, X);
  double T = Tt / (1.0 + 0.5 * (g - 1.0) * M * M);
  for (int it = 0; it < 100; ++it) {
    const double f = F(T);
    if (f > 0.0) {
      hi = T;
    } else {
      lo = T;
    }
    const double fp = cp_mass(T, X) + 0.5 * M * M * sound_speed_sq_dT(T, X);
    double T_new = (fp > 0.0) ? T - f / fp : 0.5 * (lo + hi);
    if (!(T_new > lo && T_new < hi)) {
      T_new = 0.5 * (lo + hi);
    }
    if (std::abs(T_new - T) <= 1e-14 * T) {
      return T_new;
    }
    T = T_new;
  }
  return T;
}

double friction(double Re, double e_D, const std::string& model, double f_mult) {
  if (model == "fixed") {
    return f_mult;
  }
  return std::get<0>(solver::friction_and_jacobian(model, Re, e_D).result) * f_mult;
}

// Composite 12-point Gauss-Legendre on kFannoPanels geometric panels. A FIXED
// rule, not an adaptive one: the length is then a smooth function of G and the
// Mach limits, so differences of it (the network Jacobians) carry no noise from
// a refinement pattern that switches as the arguments move. Measured against a
// 64-panel reference over inlet Mach 0.01-0.99: <= 2e-10 at constant f; with a
// Re-dependent f the floor is ~1e-8, set by kinks in mu(T), not by the rule.
constexpr int kFannoPanels = 4;
constexpr std::array<double, 6> kGlX = {0.12523340851146894, 0.3678314989981802,
                                        0.5873179542866175,  0.7699026741943047,
                                        0.9041172563704748,  0.9815606342467192};
constexpr std::array<double, 6> kGlW = {0.24914704581340297, 0.2334925365383549,
                                        0.2031674267230659,  0.16007832854334636,
                                        0.10693932599531876, 0.04717533638651128};

template <typename Fn>
double composite_gl(const Fn& f, double a, double b) {
  // Geometric panel edges: the interval in s = 1/M^2 spans orders of magnitude
  // at low Mach.
  const double ratio = std::pow(b / a, 1.0 / kFannoPanels);
  double total = 0.0;
  double lo = a;
  for (int k = 0; k < kFannoPanels; ++k) {
    const double hi = (k == kFannoPanels - 1) ? b : lo * ratio;
    const double c = 0.5 * (lo + hi);
    const double h = 0.5 * (hi - lo);
    double acc = 0.0;
    for (std::size_t j = 0; j < kGlX.size(); ++j) {
      acc += kGlW[j] * (f(c - h * kGlX[j]) + f(c + h * kGlX[j]));
    }
    total += acc * h;
    lo = hi;
  }
  return total;
}

// Brent's method for a sign change of fn on [a, b].
template <typename Fn>
double brent(const Fn& fn, double a, double b, double fa, double fb, double rel_tol) {
  double c = a, fc = fa, d = b - a, e = d;
  for (int it = 0; it < 200; ++it) {
    if ((fb > 0.0) == (fc > 0.0)) {
      c = a;
      fc = fa;
      d = e = b - a;
    }
    if (std::abs(fc) < std::abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    const double tol = 2.0 * 1e-16 * std::abs(b) + 0.5 * rel_tol * std::abs(b);
    const double m = 0.5 * (c - b);
    if (std::abs(m) <= tol || fb == 0.0) {
      return b;
    }
    if (std::abs(e) >= tol && std::abs(fa) > std::abs(fb)) {
      double p, q, r;
      const double s = fb / fa;
      if (a == c) {
        p = 2.0 * m * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) {
        q = -q;
      } else {
        p = -p;
      }
      if (2.0 * p < std::min(3.0 * m * q - std::abs(tol * q), std::abs(e * q))) {
        e = d;
        d = p / q;
      } else {
        d = m;
        e = m;
      }
    } else {
      d = m;
      e = m;
    }
    a = b;
    fa = fb;
    b += (std::abs(d) > tol) ? d : (m > 0.0 ? tol : -tol);
    fb = fn(b);
  }
  return b;
}

// Isentropic static pressure at T from the stagnation state: for an ideal gas
// at fixed composition s(T, P) = s(T, Pref) - R ln(P / Pref), so equal entropy
// gives P in closed form.
double isentropic_P(double Pt, double Tt, double T, const std::vector<double>& X) {
  const double R = specific_gas_constant(X);
  return Pt * std::exp((s_mass(T, X, 101325.0) - s_mass(Tt, X, 101325.0)) / R);
}

double isentropic_flux(double Pt, double Tt, double M, const std::vector<double>& X) {
  if (M <= 0.0) {
    return 0.0;
  }
  const double T = static_T(Tt, M, X);
  const double P = isentropic_P(Pt, Tt, T, X);
  return density(T, P, X) * M * std::sqrt(sound_speed_sq(T, X));
}

void check_duct(double L, double D) {
  if (!(L > 0.0) || !(D > 0.0)) {
    throw std::invalid_argument("fanno (Mach): L and D must be positive");
  }
}

}  // namespace

FannoMachState fanno_state_at_mach(double G, double Tt, double M,
                                   const std::vector<double>& X) {
  if (!(M > 0.0) || !(G > 0.0) || !(Tt > 0.0)) {
    throw std::invalid_argument("fanno_state_at_mach: G, Tt and M must be positive");
  }
  FannoMachState s;
  s.M = M;
  s.T = static_T(Tt, M, X);
  s.a = std::sqrt(sound_speed_sq(s.T, X));
  s.u = M * s.a;
  s.rho = G / s.u;
  s.P = s.rho * specific_gas_constant(X) * s.T;
  return s;
}

double fanno_dx_dmach(double G, double Tt, double M, const std::vector<double>& X,
                      double D, double roughness, const std::string& friction_model,
                      double f_multiplier) {
  const FannoMachState s = fanno_state_at_mach(G, Tt, M, X);
  const double R = specific_gas_constant(X);
  const double a2 = s.a * s.a;
  const double da2_dT = sound_speed_sq_dT(s.T, X);
  // Energy h(T) + M^2 a^2 / 2 = h0, differentiated at fixed h0.
  const double dT_dM = -M * a2 / (cp_mass(s.T, X) + 0.5 * M * M * da2_dT);
  const double du_dM = s.a + M * da2_dT / (2.0 * s.a) * dT_dM;
  const double drho_dM = -G * du_dM / (s.u * s.u);
  const double dP_dM = R * (drho_dM * s.T + s.rho * dT_dM);
  const double Re = G * D / viscosity(s.T, s.P, X);
  const double f = friction(Re, roughness / D, friction_model, f_multiplier);
  return -(dP_dM + G * du_dM) * 2.0 * D / (f * G * s.u);
}

double fanno_length_between(double G, double Tt, double M1, double M2,
                            const std::vector<double>& X, double D, double roughness,
                            const std::string& friction_model, double f_multiplier) {
  if (!(M1 > 0.0) || M2 < M1 || M2 > 1.0) {
    throw std::invalid_argument("fanno_length_between: need 0 < M1 <= M2 <= 1");
  }
  if (M2 == M1) {
    return 0.0;
  }
  // In s = 1/M^2, dx/ds = (M^3 / 2) dx/dM up to sign: bounded at low Mach,
  // where dx/dM itself grows as 1/M^3.
  auto integrand = [&](double s) {
    const double M = 1.0 / std::sqrt(s);
    return 0.5 * M * M * M *
           fanno_dx_dmach(G, Tt, M, X, D, roughness, friction_model, f_multiplier);
  };
  return composite_gl(integrand, 1.0 / (M2 * M2), 1.0 / (M1 * M1));
}

double fanno_sonic_mass_flux(double Pt, double Tt, const std::vector<double>& X) {
  return isentropic_flux(Pt, Tt, 1.0, X);
}

double fanno_inlet_mach(double Pt, double Tt, double G, const std::vector<double>& X) {
  if (G <= 0.0) {
    return 0.0;
  }
  const double G_max = fanno_sonic_mass_flux(Pt, Tt, X);
  if (G >= G_max) {
    return 1.0;
  }
  // Newton on the isentropic flux, which rises monotonically to G_max at
  // M = 1, kept inside the bracket [lo, hi] (bisection when a step leaves it).
  // Starts from the low-Mach limit G ~ G_max M (rho and a near stagnation).
  auto fn = [&](double M) { return isentropic_flux(Pt, Tt, M, X) - G; };
  double lo = 0.0, hi = 1.0;
  double M = std::min(G / G_max, 0.9);
  for (int it = 0; it < 60; ++it) {
    const double f = fn(M);
    if (f > 0.0) {
      hi = M;
    } else {
      lo = M;
    }
    const double h = 1e-7 * std::max(M, 1e-3);
    const double fp = (fn(std::min(M + h, 1.0)) - fn(M - h)) / (std::min(M + h, 1.0) - (M - h));
    double M_new = (fp > 0.0) ? M - f / fp : 0.5 * (lo + hi);
    if (!(M_new > lo && M_new < hi)) {
      M_new = 0.5 * (lo + hi);
    }
    if (std::abs(M_new - M) <= 1e-15 * M + 1e-300) {
      return M_new;
    }
    M = M_new;
  }
  return M;
}

FannoDuctResult fanno_duct(double Pt, double Tt, double G, const std::vector<double>& X,
                           double L, double D, double roughness,
                           const std::string& friction_model, double f_multiplier) {
  check_duct(L, D);
  if (!(G > 0.0)) {
    throw std::invalid_argument("fanno_duct: G must be positive");
  }
  FannoDuctResult out;
  out.G = G;
  const double M_in = fanno_inlet_mach(Pt, Tt, G, X);
  out.inlet = fanno_state_at_mach(G, Tt, M_in, X);
  if (M_in >= 1.0) {
    out.choked = true;
    out.L_star = 0.0;
    out.exit = out.inlet;
    out.Pt_exit = P0_from_static(out.exit.P, out.exit.T, 1.0, X);
    return out;
  }
  out.L_star = fanno_length_between(G, Tt, M_in, 1.0, X, D, roughness, friction_model,
                                    f_multiplier);
  double M_exit = 1.0;
  if (out.L_star <= L) {
    out.choked = out.L_star < L;
  } else {
    // Exit Mach: the length from M_in reaches L. Newton on the cumulative
    // length, whose slope is dx/dM exactly, inside a bracket.
    auto length_to = [&](double Me) {
      return fanno_length_between(G, Tt, M_in, Me, X, D, roughness, friction_model,
                                  f_multiplier) -
             L;
    };
    M_exit = brent(length_to, M_in, 1.0, -L, out.L_star - L, kFannoMachRelTol);
  }
  out.exit = fanno_state_at_mach(G, Tt, M_exit, X);
  out.Pt_exit = P0_from_static(out.exit.P, out.exit.T, M_exit, X, 1e-13, 100);
  return out;
}

double fanno_choked_mass_flux(double Pt, double Tt, const std::vector<double>& X,
                              double L, double D, double roughness,
                              const std::string& friction_model, double f_multiplier) {
  check_duct(L, D);
  const double G_max = fanno_sonic_mass_flux(Pt, Tt, X);
  // L*(G) falls monotonically from infinity (G -> 0) to 0 (G -> G_max).
  auto excess = [&](double logG) {
    const double G = std::exp(logG);
    const double M_in = fanno_inlet_mach(Pt, Tt, G, X);
    if (M_in >= 1.0) {
      return -L;
    }
    return fanno_length_between(G, Tt, M_in, 1.0, X, D, roughness, friction_model,
                                f_multiplier) -
           L;
  };
  double hi = std::log(G_max * (1.0 - 1e-12));
  double lo = hi - 1.0;
  double f_lo = excess(lo);
  while (f_lo <= 0.0 && lo > hi - 60.0) {
    lo -= 2.0;
    f_lo = excess(lo);
  }
  const double f_hi = excess(hi);
  if (f_hi >= 0.0) {
    return G_max;
  }
  return std::exp(brent(excess, lo, hi, f_lo, f_hi, 1e-14));
}

namespace {

// Exit pressure per unit mass flux: at fixed (M, Tt) the station state's P is
// G R T / (M a), linear in G, and so is its stagnation pressure (P times the
// isentropic ratio, which depends on T and M only). So for a target pressure
// the flux at any exit Mach is explicit: G = P_target / exit_per_flux(M).
double exit_per_flux(double Me, double Tt, const std::vector<double>& X, bool exit_total) {
  const FannoMachState s = fanno_state_at_mach(1.0, Tt, Me, X);
  if (!exit_total) {
    return s.P;
  }
  const double R = specific_gas_constant(X);
  return s.P * std::exp((s_mass(Tt, X, 101325.0) - s_mass(s.T, X, 101325.0)) / R);
}

struct FlowProblem {
  double Pt0, Tt0, P, L, D, rough, fmult;
  const std::vector<double>* X;
  const std::string* model;
  bool exit_total;
  double Me_guess;
};

// Length the flux G(Me) = P / exit_per_flux(Me) needs to reach exit Mach Me,
// less L. Falls monotonically with Me: positive for a slow exit (a long duct
// would be needed), negative once G passes the choked flux. A flux the inlet
// cannot supply (sonic inlet) needs no length at all: -L.
double psi(double Me, double P, double Pt0, double Tt0, const FlowProblem& pb) {
  const std::vector<double>& X = *pb.X;
  const double G = P / exit_per_flux(Me, Tt0, X, pb.exit_total);
  if (G >= fanno_sonic_mass_flux(Pt0, Tt0, X)) {
    return -pb.L;
  }
  const double Mi = fanno_inlet_mach(Pt0, Tt0, G, X);
  if (Mi >= Me) {
    return -pb.L;
  }
  return fanno_length_between(G, Tt0, Mi, Me, X, pb.D, pb.rough, *pb.model, pb.fmult) - pb.L;
}

// The raw (unregularised) solve at drive Pt0 - P > 0.
FannoChannelFlow raw_flow(const FlowProblem& pb) {
  const std::vector<double>& X = *pb.X;
  FannoChannelFlow out;
  const double choked_test = psi(1.0, pb.P, pb.Pt0, pb.Tt0, pb);
  if (choked_test >= 0.0) {
    // Choked: even a sonic exit at this back pressure would need more duct
    // than there is, so the exit sits at sonic above P and the flux is the
    // choked one, set by the inlet state alone.
    out.choked = true;
    out.M_exit = 1.0;
    out.G = fanno_choked_mass_flux(pb.Pt0, pb.Tt0, X, pb.L, pb.D, pb.rough, *pb.model, pb.fmult);
    out.M_in = fanno_inlet_mach(pb.Pt0, pb.Tt0, out.G, X);
    // L*(G; Pt0, Tt0) = L, differentiated implicitly.
    auto Lstar = [&](double G, double Pt0, double Tt0) {
      const double Mi = fanno_inlet_mach(Pt0, Tt0, G, X);
      return fanno_length_between(G, Tt0, Mi, 1.0, X, pb.D, pb.rough, *pb.model, pb.fmult);
    };
    const double hG = 1e-6 * out.G, hP = 1e-6 * pb.Pt0, hT = 1e-6 * pb.Tt0;
    const double LG = (Lstar(out.G + hG, pb.Pt0, pb.Tt0) - Lstar(out.G - hG, pb.Pt0, pb.Tt0)) / (2 * hG);
    const double LP = (Lstar(out.G, pb.Pt0 + hP, pb.Tt0) - Lstar(out.G, pb.Pt0 - hP, pb.Tt0)) / (2 * hP);
    const double LT = (Lstar(out.G, pb.Pt0, pb.Tt0 + hT) - Lstar(out.G, pb.Pt0, pb.Tt0 - hT)) / (2 * hT);
    out.dG_dPt0 = -LP / LG;
    out.dG_dTt0 = -LT / LG;
    out.dG_dP_target = 0.0;
    return out;
  }
  // Unchoked: the exit Mach is the root of psi, which falls monotonically in
  // Me. A caller's guess (the previous solve) gives a tight bracket, widened
  // geometrically until it holds a sign change; without one the bracket is
  // grown down from Me = 1 by halving.
  auto f = [&](double Me) { return psi(Me, pb.P, pb.Pt0, pb.Tt0, pb); };
  double lo, hi, f_lo, f_hi;
  if (pb.Me_guess > 0.0 && pb.Me_guess < 1.0) {
    double w = 0.01;
    lo = pb.Me_guess * (1.0 - w);
    hi = std::min(1.0, pb.Me_guess * (1.0 + w));
    f_lo = f(lo);
    f_hi = f(hi);
    while ((f_lo <= 0.0 || f_hi > 0.0) && w < 1.0) {
      w *= 4.0;
      if (f_lo <= 0.0) {
        hi = lo;
        f_hi = f_lo;
        lo = std::max(pb.Me_guess * (1.0 - w), 1e-6);
        f_lo = f(lo);
      } else {
        lo = hi;
        f_lo = f_hi;
        hi = std::min(1.0, pb.Me_guess * (1.0 + w));
        f_hi = f(hi);
      }
    }
  } else {
    hi = 1.0;
    f_hi = choked_test;
    lo = 0.5;
    f_lo = f(lo);
  }
  while (f_lo <= 0.0 && lo > 1e-6) {
    hi = lo;
    f_hi = f_lo;
    lo *= 0.5;
    f_lo = f(lo);
  }
  const double Me = brent(f, lo, hi, f_lo, f_hi, 1e-13);
  const double phi = exit_per_flux(Me, pb.Tt0, X, pb.exit_total);
  out.M_exit = Me;
  out.G = pb.P / phi;
  out.M_in = fanno_inlet_mach(pb.Pt0, pb.Tt0, out.G, X);
  // psi(Me; P, Pt0, Tt0) = 0, differentiated implicitly. Its slope in Me stays
  // finite at sonic -- dx/dM -> 0 there, but G(Me) still moves -- so this is
  // well conditioned right up to the choke, where it meets the choked branch.
  //
  // psi(Me, P) = Lam(G, Me) with G = P / phi(Me, Tt0) and
  // Lam = (length from M_in(G, Pt0, Tt0) to Me) - L, so most partials are
  // exact: d Lam/d Me = dx/dM at Me (the upper limit), and Pt0 acts only
  // through M_in, with d Lam/d M_in = -dx/dM at M_in. Only d Lam/d G and
  // d Lam/d Tt0 are differenced, each on smooth fixed-rule integrals.
  const double G = out.G;
  const double Mi = out.M_in;
  const std::string& model = *pb.model;
  auto Lam = [&](double Gx, double Ptx, double Ttx) {
    const double Mix = fanno_inlet_mach(Ptx, Ttx, Gx, X);
    return fanno_length_between(Gx, Ttx, Mix, Me, X, pb.D, pb.rough, model, pb.fmult) - pb.L;
  };
  const double hG = 1e-6 * G;
  const double hT = 1e-6 * pb.Tt0;
  const double Lam_G = (Lam(G + hG, pb.Pt0, pb.Tt0) - Lam(G - hG, pb.Pt0, pb.Tt0)) / (2.0 * hG);
  const double Lam_T = (Lam(G, pb.Pt0, pb.Tt0 + hT) - Lam(G, pb.Pt0, pb.Tt0 - hT)) / (2.0 * hT);
  const double Lam_Me = fanno_dx_dmach(G, pb.Tt0, Me, X, pb.D, pb.rough, model, pb.fmult);
  // Isentropic inlet: flux(Pt0, Tt0, M_in) = G with flux proportional to Pt0
  // at fixed (Tt0, M), so d M_in/d Pt0 = -(G / Pt0) / (d flux/d M).
  const double hMi = 1e-7 * std::max(Mi, 1e-4);
  const double flux_M = (isentropic_flux(pb.Pt0, pb.Tt0, Mi + hMi, X) -
                         isentropic_flux(pb.Pt0, pb.Tt0, Mi - hMi, X)) /
                        (2.0 * hMi);
  const double dMi_dPt = -(G / pb.Pt0) / flux_M;
  const double Lam_Pt = -fanno_dx_dmach(G, pb.Tt0, Mi, X, pb.D, pb.rough, model, pb.fmult) * dMi_dPt;
  // G = P / phi(Me, Tt0)
  const double hM = 1e-6 * Me;
  const double Mp = std::min(Me + hM, 1.0);
  const double Mm = Mp - 2.0 * hM;
  const double phi_M = (exit_per_flux(Mp, pb.Tt0, X, pb.exit_total) -
                        exit_per_flux(Mm, pb.Tt0, X, pb.exit_total)) / (Mp - Mm);
  const double phi_T = (exit_per_flux(Me, pb.Tt0 + hT, X, pb.exit_total) -
                        exit_per_flux(Me, pb.Tt0 - hT, X, pb.exit_total)) / (2.0 * hT);
  const double G_P = 1.0 / phi;
  const double G_M = -pb.P * phi_M / (phi * phi);
  const double G_T = -pb.P * phi_T / (phi * phi);
  const double psi_M = Lam_G * G_M + Lam_Me;
  const double psi_P = Lam_G * G_P;
  const double psi_Pt = Lam_Pt;
  const double psi_T = Lam_G * G_T + Lam_T;
  const double dMe_dP = -psi_P / psi_M;
  const double dMe_dPt = -psi_Pt / psi_M;
  const double dMe_dT = -psi_T / psi_M;
  out.dG_dP_target = G_P + G_M * dMe_dP;
  out.dG_dPt0 = G_M * dMe_dPt;
  out.dG_dTt0 = G_T + G_M * dMe_dT;
  return out;
}

}  // namespace

FannoChannelFlow fanno_channel_flow(double Pt0, double Tt0, const std::vector<double>& X,
                                    double P_target, bool exit_total, double L, double D,
                                    double roughness, const std::string& friction_model,
                                    double f_multiplier, double M_exit_guess) {
  check_duct(L, D);
  FannoChannelFlow out;
  const double drive = Pt0 - P_target;
  if (!(drive > 0.0) || !(Pt0 > 0.0) || !(Tt0 > 0.0)) {
    return out;
  }
  FlowProblem pb{Pt0, Tt0, P_target, L, D, roughness, f_multiplier, &X, &friction_model, exit_total,
                 M_exit_guess};
  if (drive >= kFannoFlowDriveFloor) {
    return raw_flow(pb);
  }
  // Below the floor: G = G0 h(t), t = drive / floor, with G0 the flux at the
  // floor and h the odd cubic matching the raw law's value AND slope there:
  // sigma = (dG/d drive) floor / G0 is the local exponent (1/2 for a sqrt law;
  // larger at low Re, where f rises as the flow falls), and
  //   h(t) = (3 - sigma)/2 t + (sigma - 1)/2 t^3,   h(1) = 1, h'(1) = sigma,
  // finite and positive at t = 0 and monotone for sigma < 3.
  pb.P = Pt0 - kFannoFlowDriveFloor;
  const FannoChannelFlow at = raw_flow(pb);
  const double sigma =
      (at.G > 0.0) ? std::clamp(-at.dG_dP_target * kFannoFlowDriveFloor / at.G, 0.0, 2.9) : 0.5;
  const double t = drive / kFannoFlowDriveFloor;
  const double h = 0.5 * (3.0 - sigma) * t + 0.5 * (sigma - 1.0) * t * t * t;
  const double dh = 0.5 * (3.0 - sigma) + 1.5 * (sigma - 1.0) * t * t;
  out = at;
  out.G = at.G * h;
  // d/dP_target through t only; d/dPt0 through t and through G0, whose own
  // target Pt0 - floor moves with Pt0. sigma is held fixed: its own
  // variation is second order within the floor.
  out.dG_dP_target = -at.G * dh / kFannoFlowDriveFloor;
  out.dG_dPt0 = at.G * dh / kFannoFlowDriveFloor + h * (at.dG_dPt0 + at.dG_dP_target);
  out.dG_dTt0 = h * at.dG_dTt0;
  return out;
}

}  // namespace combaero
