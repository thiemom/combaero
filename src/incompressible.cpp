#include "incompressible.h"

#include "friction.h"
#include "geometry.h"
#include "math_constants.h" // Provides M_PI for MSVC
#include "thermo.h"
#include "transport.h"
#include "input_validation.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace combaero {

// -------------------------------------------------------------
// Bernoulli equation
// -------------------------------------------------------------

double bernoulli_P2(double P1, double v1, double v2, double rho, double dz,
                    double g) {
  // P1 + 0.5*rho*v1^2 + rho*g*h1 = P2 + 0.5*rho*v2^2 + rho*g*h2
  // P2 = P1 + 0.5*rho*(v1^2 - v2^2) - rho*g*dz
  return P1 + 0.5 * rho * (v1 * v1 - v2 * v2) - rho * g * dz;
}

double bernoulli_v2(double P1, double P2, double v1, double rho, double dz,
                    double g) {
  // P1 + 0.5*rho*v1^2 = P2 + 0.5*rho*v2^2 + rho*g*dz
  // v2^2 = v1^2 + 2(P1 - P2)/rho - 2g*dz
  double v2_sq = v1 * v1 + 2.0 * (P1 - P2) / rho - 2.0 * g * dz;
  if (!validation::is_finite_nonnegative(v2_sq)) {
    throw std::invalid_argument(
        "bernoulli_v2: insufficient energy for flow (v2^2 < 0)");
  }
  return std::sqrt(v2_sq);
}

// -------------------------------------------------------------
// Orifice / restriction flow
// -------------------------------------------------------------

double orifice_mdot(double P1, double P2, double A, double Cd, double rho) {
  if (!validation::is_finite(P1) || !validation::is_finite(P2) || P1 < P2) {
    throw std::invalid_argument("orifice_mdot: P1 must be >= P2");
  }
  if (!validation::is_finite_positive(A) || !validation::is_finite_positive(Cd) || !validation::is_finite_positive(rho)) {
    throw std::invalid_argument("orifice_mdot: A, Cd, rho must be positive");
  }
  // m_dot = Cd * A * sqrt(2 * rho * Delta_P)
  double dP = P1 - P2;
  return Cd * A * std::sqrt(2.0 * rho * dP);
}

double orifice_Q(double P1, double P2, double A, double Cd, double rho) {
  // Q = m_dot / rho
  return orifice_mdot(P1, P2, A, Cd, rho) / rho;
}

double orifice_velocity(double P1, double P2, double rho) {
  if (!validation::is_finite(P1) || !validation::is_finite(P2) || P1 < P2) {
    throw std::invalid_argument("orifice_velocity: P1 must be >= P2");
  }
  if (!validation::is_finite_positive(rho)) {
    throw std::invalid_argument("orifice_velocity: rho must be positive");
  }
  // v = sqrt(2 * Delta_P / rho)
  return std::sqrt(2.0 * (P1 - P2) / rho);
}

double orifice_area(double mdot, double P1, double P2, double Cd, double rho) {
  if (!validation::is_finite(P1) || !validation::is_finite(P2) || P1 <= P2) {
    throw std::invalid_argument("orifice_area: P1 must be > P2");
  }
  if (!validation::is_finite_positive(mdot) || !validation::is_finite_positive(Cd) || !validation::is_finite_positive(rho)) {
    throw std::invalid_argument("orifice_area: mdot, Cd, rho must be positive");
  }
  // A = m_dot / (Cd * sqrt(2 * rho * Delta_P))
  double dP = P1 - P2;
  return mdot / (Cd * std::sqrt(2.0 * rho * dP));
}

double orifice_dP(double mdot, double A, double Cd, double rho) {
  if (!validation::is_finite_nonnegative(mdot) || !validation::is_finite_positive(A) || !validation::is_finite_positive(Cd) || !validation::is_finite_positive(rho)) {
    throw std::invalid_argument("orifice_dP: invalid parameters");
  }
  // Delta_P = (m_dot / (Cd * A))^2 / (2 * rho)
  double term = mdot / (Cd * A);
  return term * term / (2.0 * rho);
}

// -------------------------------------------------------------
// Pipe flow (Darcy-Weisbach)
// -------------------------------------------------------------

double channel_dP(double v, double L, double D, double f, double rho) {
  if (!validation::is_finite_positive(D) || !validation::is_finite_nonnegative(L) || !validation::is_finite_nonnegative(f) || !validation::is_finite_positive(rho)) {
    throw std::invalid_argument("channel_dP: invalid parameters");
  }
  // Delta_P = f * (L/D) * (rho * v^2 / 2)
  return f * (L / D) * 0.5 * rho * v * v;
}

double channel_dP_mdot(double mdot, double L, double D, double f, double rho) {
  double v = channel_velocity(mdot, D, rho);
  return channel_dP(v, L, D, f, rho);
}

double channel_velocity(double mdot, double D, double rho) {
  if (!validation::is_finite_positive(D) || !validation::is_finite_positive(rho)) {
    throw std::invalid_argument("channel_velocity: D and rho must be positive");
  }
  // v = m_dot / (rho * A) where A = pi*D^2/4
  double A = M_PI * D * D / 4.0;
  return mdot / (rho * A);
}

double channel_mdot(double v, double D, double rho) {
  if (!validation::is_finite_positive(D) || !validation::is_finite_positive(rho)) {
    throw std::invalid_argument("channel_mdot: D and rho must be positive");
  }
  // m_dot = rho * v * A where A = pi*D^2/4
  double A = M_PI * D * D / 4.0;
  return rho * v * A;
}

// -------------------------------------------------------------
// Hydraulic utilities
// -------------------------------------------------------------

double dynamic_pressure(double v, double rho) {
  // q = 0.5 * rho * v^2
  return 0.5 * rho * v * v;
}

double velocity_from_q(double q, double rho) {
  if (!validation::is_finite_nonnegative(q) || !validation::is_finite_positive(rho)) {
    throw std::invalid_argument("velocity_from_q: q >= 0 and rho > 0 required");
  }
  // v = sqrt(2 * q / rho)
  return std::sqrt(2.0 * q / rho);
}

// -------------------------------------------------------------
// Pressure loss coefficient (zeta / K)
// -------------------------------------------------------------

double pressure_loss(double v, double rho, double zeta) {
  if (!validation::is_finite_positive(rho) || !validation::is_finite_nonnegative(zeta)) {
    throw std::invalid_argument(
        "pressure_loss: rho > 0 and zeta >= 0 required");
  }
  // Delta_P = zeta * 0.5*rho*v^2
  return zeta * 0.5 * rho * v * v;
}

double velocity_from_pressure_loss(double dP, double rho, double zeta) {
  if (!validation::is_finite_nonnegative(dP) || !validation::is_finite_positive(rho) || !validation::is_finite_positive(zeta)) {
    throw std::invalid_argument(
        "velocity_from_pressure_loss: dP >= 0, rho > 0, zeta > 0 required");
  }
  // v = sqrt(2*Delta_P / (zeta*rho))
  return std::sqrt(2.0 * dP / (zeta * rho));
}

double zeta_from_Cd(double Cd) {
  if (!validation::is_finite_positive(Cd)) {
    throw std::invalid_argument("zeta_from_Cd: Cd must be positive");
  }
  // zeta = 1 / Cd^2
  return 1.0 / (Cd * Cd);
}

double Cd_from_zeta(double zeta) {
  if (!validation::is_finite_positive(zeta)) {
    throw std::invalid_argument("Cd_from_zeta: zeta must be positive");
  }
  // Cd = 1 / sqrt(zeta)
  return 1.0 / std::sqrt(zeta);
}

// -------------------------------------------------------------
// High-level thermo-aware API
// -------------------------------------------------------------

// Helper: dispatch friction factor from correlation name
static double friction_from_correlation(double Re, double e_D,
                                        const std::string &correlation) {
  if (correlation == "haaland")
    return friction_haaland(Re, e_D);
  if (correlation == "serghides")
    return friction_serghides(Re, e_D);
  if (correlation == "colebrook")
    return friction_colebrook(Re, e_D);
  if (correlation == "petukhov")
    return friction_petukhov(Re);
  throw std::invalid_argument(
      "incompressible: unknown friction correlation '" + correlation +
      "'. "
      "Valid options: 'haaland', 'serghides', 'colebrook', 'petukhov'");
}

IncompressibleFlowSolution orifice_flow_thermo(double T, double P,
                                               const std::vector<double> &X,
                                               double P_back, double A,
                                               double Cd) {
  if (!validation::is_finite_positive(T))
    throw std::invalid_argument("orifice_flow: T must be positive");
  if (!validation::is_finite_positive(P))
    throw std::invalid_argument("orifice_flow: P must be positive");
  if (!validation::is_finite_nonnegative(P_back))
    throw std::invalid_argument("orifice_flow: P_back must be non-negative");
  if (P < P_back)
    throw std::invalid_argument("orifice_flow: P must be >= P_back");
  if (!validation::is_finite_positive(A))
    throw std::invalid_argument("orifice_flow: A must be positive");
  if (!validation::is_finite_positive(Cd) || Cd > 1.0)
    throw std::invalid_argument("orifice_flow: Cd must be in (0, 1]");

  const double rho = density(T, P, X);
  const double dP = P - P_back;
  const double v = (dP > 0.0) ? std::sqrt(2.0 * dP / rho) : 0.0;
  const double mdot = Cd * A * rho * v;

  IncompressibleFlowSolution sol;
  sol.mdot = mdot;
  sol.v = Cd * v; // effective throat velocity (Cd * ideal velocity)
  sol.dP = dP;
  sol.Re = 0.0; // not meaningful for an orifice (no length scale)
  sol.rho = rho;
  sol.f = Cd; // store Cd in f field for orifice results
  return sol;
}

// Helper to build the symmetric structural profile array over length L
static std::vector<IncompressibleStation>
build_incompressible_profile(double L, double P_in, double dP, double T,
                             double rho, double v, double h_val,
                             std::size_t n_steps) {
  std::vector<IncompressibleStation> profile;
  if (n_steps < 2)
    n_steps = 2; // Need at least start and end

  profile.reserve(n_steps);

  // Incompressible properties (v, rho, T, h, M) are constant throughout.
  // Pressure drops linearly due to uniform friction assumption.
  const double n_steps_d = static_cast<double>(n_steps);
  const double dx = L / (n_steps_d - 1.0);
  const double dP_step = dP / (n_steps_d - 1.0);

  // Since this is a low-Mach assumption, M approx 0.
  const double M = 0.0;

  for (std::size_t i = 0; i < n_steps; ++i) {
    IncompressibleStation st;
    st.x = i * dx;
    st.P = P_in - (i * dP_step);
    st.T = T;
    st.rho = rho;
    st.v = v;
    st.M = M;
    st.h = h_val;
    profile.push_back(st);
  }
  return profile;
}

IncompressibleFlowSolution
channel_flow(double T, double P, const std::vector<double> &X, double u,
             double L, double D, double f, std::size_t n_steps,
             bool store_profile) {
  if (!validation::is_finite_positive(T))
    throw std::invalid_argument("channel_flow: T must be positive");
  if (!validation::is_finite_positive(P))
    throw std::invalid_argument("channel_flow: P must be positive");
  if (!validation::is_finite_nonnegative(u))
    throw std::invalid_argument("channel_flow: u must be non-negative");
  if (!validation::is_finite_nonnegative(L))
    throw std::invalid_argument("channel_flow: L must be non-negative");
  if (!validation::is_finite_positive(D))
    throw std::invalid_argument("channel_flow: D must be positive");
  if (!validation::is_finite_nonnegative(f))
    throw std::invalid_argument("channel_flow: f must be non-negative");

  const double rho = density(T, P, X);
  const double mu = viscosity(T, P, X);
  const double Re = (mu > 0.0) ? rho * u * D / mu : 0.0;
  const double dP = channel_dP(u, L, D, f, rho);
  const double A = M_PI * D * D / 4.0;
  const double mdot = rho * u * A;

  const double h_val = h_mass(T, X); // [J/kg], as IncompressibleStation::h

  IncompressibleFlowSolution sol;
  sol.inlet.set_TPX(T, P, X);
  sol.outlet.set_TPX(T, P - dP, X); // Assume adiabatic for now
  sol.mdot = mdot;
  sol.v = u;
  sol.dP = dP;
  sol.Re = Re;
  sol.rho = rho;
  sol.f = f;

  if (store_profile) {
    sol.profile =
        build_incompressible_profile(L, P, dP, T, rho, u, h_val, n_steps);
  }
  return sol;
}

IncompressibleFlowSolution
channel_flow_rough(double T, double P, const std::vector<double> &X, double u,
                   double L, double D, double roughness,
                   const std::string &correlation, std::size_t n_steps,
                   bool store_profile) {
  if (!validation::is_finite_positive(T))
    throw std::invalid_argument("channel_flow_rough: T must be positive");
  if (!validation::is_finite_positive(P))
    throw std::invalid_argument("channel_flow_rough: P must be positive");
  if (!validation::is_finite_nonnegative(u))
    throw std::invalid_argument("channel_flow_rough: u must be non-negative");
  if (!validation::is_finite_nonnegative(L))
    throw std::invalid_argument("channel_flow_rough: L must be non-negative");
  if (!validation::is_finite_positive(D))
    throw std::invalid_argument("channel_flow_rough: D must be positive");
  if (!validation::is_finite_nonnegative(roughness))
    throw std::invalid_argument(
        "channel_flow_rough: roughness must be non-negative");

  const double rho = density(T, P, X);
  const double mu = viscosity(T, P, X);
  const double Re = (mu > 0.0) ? rho * u * D / mu : 0.0;
  const double e_D = roughness / D;
  const double f = friction_from_correlation(Re, e_D, correlation);
  const double dP = channel_dP(u, L, D, f, rho);
  const double A = M_PI * D * D / 4.0;
  const double mdot = rho * u * A;

  const double h_val = h_mass(T, X); // [J/kg], as IncompressibleStation::h

  IncompressibleFlowSolution sol;
  sol.inlet.set_TPX(T, P, X);
  sol.outlet.set_TPX(T, P - dP, X); // Assume adiabatic for now
  sol.mdot = mdot;
  sol.v = u;
  sol.dP = dP;
  sol.Re = Re;
  sol.rho = rho;
  sol.f = f;

  if (store_profile) {
    sol.profile =
        build_incompressible_profile(L, P, dP, T, rho, u, h_val, n_steps);
  }
  return sol;
}

// Orifice overload with user-supplied Cd correlation.
// Re is computed from the ideal throat velocity (Cd=1) so the correlation
// receives a meaningful Reynolds number without requiring iteration.
IncompressibleFlowSolution
orifice_flow_thermo(double T, double P, const std::vector<double> &X,
                    double P_back, double A, const IncompressibleCdFn &cd_fn) {
  if (!validation::is_finite_positive(T))
    throw std::invalid_argument("orifice_flow: T must be positive");
  if (!validation::is_finite_positive(P))
    throw std::invalid_argument("orifice_flow: P must be positive");
  if (!validation::is_finite_nonnegative(P_back))
    throw std::invalid_argument("orifice_flow: P_back must be non-negative");
  if (P < P_back)
    throw std::invalid_argument("orifice_flow: P must be >= P_back");
  if (!validation::is_finite_positive(A))
    throw std::invalid_argument("orifice_flow: A must be positive");

  const double rho = density(T, P, X);
  const double mu = viscosity(T, P, X);
  const double dP = P - P_back;
  const double v_ideal = (dP > 0.0) ? std::sqrt(2.0 * dP / rho) : 0.0;

  // Characteristic length: equivalent diameter from area
  const double D_eq = std::sqrt(4.0 * A / M_PI);
  const double Re = (mu > 0.0) ? rho * v_ideal * D_eq / mu : 0.0;

  const double Cd = cd_fn(T, P, X, Re);
  if (!validation::is_finite_positive(Cd) || Cd > 1.0)
    throw std::runtime_error(
        "orifice_flow: Cd correlation returned value outside (0, 1]");

  const double mdot = Cd * A * rho * v_ideal;

  IncompressibleFlowSolution sol;
  sol.mdot = mdot;
  sol.v = Cd * v_ideal;
  sol.dP = dP;
  sol.Re = Re;
  sol.rho = rho;
  sol.f = Cd;
  return sol;
}

// Channel overload with additional K-loss correlation (bends, fittings, etc.).
IncompressibleFlowSolution
channel_flow_rough(double T, double P, const std::vector<double> &X, double u,
                   double L, double D, double roughness,
                   const std::string &correlation,
                   const IncompressibleKLossFn &k_loss_fn, std::size_t n_steps,
                   bool store_profile) {
  // Compute base friction solution
  IncompressibleFlowSolution sol = channel_flow_rough(
      T, P, X, u, L, D, roughness, correlation, n_steps, store_profile);

  // Evaluate additional K-loss at the computed Re (inlet state, loop-free)
  const double K = k_loss_fn(T, P, X, sol.Re);
  const double dP_extra = K * 0.5 * sol.rho * u * u;
  sol.dP += dP_extra;
  // The outlet and the profile carry the total drop, not the friction drop
  // alone (they were left at P - dP_friction).
  sol.outlet.set_TPX(T, P - sol.dP, X);
  if (store_profile) {
    sol.profile = build_incompressible_profile(L, P, sol.dP, T, sol.rho, u,
                                               h_mass(T, X), n_steps);
  }
  return sol;
}

std::tuple<double, double, double>
channel_pressure_drop(double T, double P, const std::vector<double> &X, double v,
                      double D, double L, double roughness,
                      const std::string &correlation) {
  if (!validation::is_finite_positive(T))
    throw std::invalid_argument("channel_pressure_drop: T must be positive");
  if (!validation::is_finite_positive(P))
    throw std::invalid_argument("channel_pressure_drop: P must be positive");
  if (!validation::is_finite_positive(D))
    throw std::invalid_argument("channel_pressure_drop: D must be positive");
  if (!validation::is_finite_nonnegative(L))
    throw std::invalid_argument("channel_pressure_drop: L must be non-negative");
  if (!validation::is_finite_nonnegative(v))
    throw std::invalid_argument("channel_pressure_drop: v must be non-negative");
  if (!validation::is_finite_nonnegative(roughness))
    throw std::invalid_argument(
        "channel_pressure_drop: roughness must be non-negative");

  const double rho = density(T, P, X);
  const double Re = reynolds(T, P, X, v, D);
  const double e_D = roughness / D;
  const double f = friction_from_correlation(Re, e_D, correlation);
  const double dP = channel_dP(v, L, D, f, rho);

  return std::make_tuple(dP, Re, f);
}

} // namespace combaero
