#include "../include/state.h"
#include "../include/thermo.h"
#include "../include/transport.h"
#include <algorithm>
#include <limits>
#include <stdexcept>

// -------------------------------------------------------------
// State property getters
// -------------------------------------------------------------

double State::mw() const { return mwmix(X); }
double State::cp() const { return ::cp(T, X); }
double State::h() const { return ::h(T, X); }
double State::s() const { return ::s(T, X, P); }
double State::cv() const { return ::cv(T, X); }
double State::u() const { return ::u(T, X); }
double State::rho() const { return ::density(T, P, X); }
double State::R() const { return ::specific_gas_constant(X); }
double State::gamma() const { return ::isentropic_expansion_coefficient(T, X); }
double State::a() const { return ::speed_of_sound(T, X); }

// -------------------------------------------------------------
// Transport properties
// -------------------------------------------------------------

double State::mu() const { return ::viscosity(T, P, X); }
double State::k() const { return ::thermal_conductivity(T, P, X); }
double State::nu() const { return ::kinematic_viscosity(T, P, X); }
double State::Pr() const { return ::prandtl(T, P, X); }
double State::alpha() const { return ::thermal_diffusivity(T, P, X); }

// -------------------------------------------------------------
// State Setters
// -------------------------------------------------------------

State &State::set_TPX(double T_new, double P_new,
                      const std::vector<double> &X_new) {
  T = T_new;
  P = P_new;
  X = X_new;
  return *this;
}

State &State::set_TPY(double T_new, double P_new,
                      const std::vector<double> &Y_new) {
  T = T_new;
  P = P_new;
  X = mass_to_mole(Y_new);
  return *this;
}

State &State::set_TP(double T_new, double P_new) {
  T = T_new;
  P = P_new;
  return *this;
}

State &State::set_DP_mass(double rho_new, double P_new) {
  if (rho_new <= 0)
    throw std::invalid_argument("Density must be positive");
  P = P_new;
  T = P_new / (rho_new * R());
  return *this;
}

State &State::set_HP_mass(double h_new, double P_new) {
  P = P_new;
  T = calc_T_from_h_mass(h_new, X, T);
  return *this;
}

State &State::set_SP_mass(double s_new, double P_new) {
  P = P_new;
  T = calc_T_from_s_mass(s_new, P, X, T);
  return *this;
}

State &State::set_UV_mass(double u_new, double v_new) {
  if (v_new <= 0)
    throw std::invalid_argument("Specific volume must be positive");
  T = calc_T_from_u_mass(u_new, X, T);
  P = R() * T / v_new;
  return *this;
}

State &State::set_SV_mass(double s_new, double v_new) {
  if (v_new <= 0)
    throw std::invalid_argument("Specific volume must be positive");
  T = calc_T_from_sv_mass(s_new, v_new, X, T);
  P = R() * T / v_new;
  return *this;
}

State &State::set_PV_mass(double P_new, double v_new) {
  if (v_new <= 0)
    throw std::invalid_argument("Specific volume must be positive");
  P = P_new;
  T = P_new * v_new / R();
  return *this;
}

State &State::set_UP_mass(double u_new, double P_new) {
  P = P_new;
  T = calc_T_from_u_mass(u_new, X, T);
  return *this;
}

State &State::set_VH_mass(double v_new, double h_new) {
  if (v_new <= 0)
    throw std::invalid_argument("Specific volume must be positive");
  T = calc_T_from_h_mass(h_new, X, T);
  P = R() * T / v_new;
  return *this;
}

State &State::set_SH_mass(double s_new, double h_new) {
  T = calc_T_from_sh_mass(s_new, h_new, X, T);
  // Find pressure such that s_mass(T, P) = s_new
  // s_mass(T, P) = s_molar(T, P) / MW
  // s(T, P) = sum(X_i * s^0_i(T)) - R * ln(P/P0) - R * sum(X_i ln X_i)
  // We want P such that s_mass(T, P) = s_new
  // MW * s_new = s^0_mix(T) - R * ln(P/P0) - R * sum(X_i ln X_i)
  // s(T, P0) = s^0_mix(T) - R * sum(X_i ln X_i)
  // mw * s_new = s(T, P0) - R * ln(P/P0)
  // ln(P/P0) = (s(T, P0) - mw * s_new) / R

  double MW = mw();                            // g/mol
  double s_target_molar = s_new * MW / 1000.0; // J/(mol K)

  double P0 = 101325.0;
  double s_at_P0 = ::s(T, X, P0, P0);
  double ln_P_P0 = (s_at_P0 - s_target_molar) / combaero::thermo::R_GAS;

  P = P0 * std::exp(ln_P_P0);
  return *this;
}

// -------------------------------------------------------------
// Stream mixing
// -------------------------------------------------------------

Stream mix(const std::vector<Stream> &streams, double P_out) {
  if (streams.empty()) {
    throw std::runtime_error("mix: at least one stream required");
  }

  if (streams.size() == 1) {
    Stream out = streams[0];
    if (P_out > 0.0) {
      out.state.P = P_out;
    }
    return out;
  }

  const std::size_t n_species = streams[0].X().size();

  // Validate all streams have same species count
  for (const auto &s : streams) {
    if (s.X().size() != n_species) {
      throw std::runtime_error(
          "mix: all streams must have same number of species");
    }
  }

  // Determine output pressure (minimum if not specified)
  if (P_out < 0.0) {
    P_out = std::numeric_limits<double>::max();
    for (const auto &s : streams) {
      P_out = std::min(P_out, s.P());
    }
  }

  // Mass balance: total mass flow
  double mdot_total = 0.0;
  for (const auto &s : streams) {
    mdot_total += s.mdot;
  }

  if (mdot_total <= 0.0) {
    throw std::runtime_error("mix: total mass flow must be positive");
  }

  // Species balance: accumulate molar flows
  // n_dot_k = sum_i (mdot_i / MW_i) * X_i[k]
  std::vector<double> n_dot(n_species, 0.0);
  double n_dot_total = 0.0;

  for (const auto &s : streams) {
    double MW_i = s.mw();                    // g/mol
    double n_dot_i = s.mdot / (MW_i * 1e-3); // mol/s (MW in kg/mol)
    for (std::size_t k = 0; k < n_species; ++k) {
      n_dot[k] += n_dot_i * s.X()[k];
    }
    n_dot_total += n_dot_i;
  }

  // Convert to mole fractions
  std::vector<double> X_out(n_species);
  for (std::size_t k = 0; k < n_species; ++k) {
    X_out[k] = n_dot[k] / n_dot_total;
  }

  // Enthalpy balance: H_out * mdot_out = sum_i (H_i * mdot_i)
  // where H is specific enthalpy [J/kg]
  double H_dot_total = 0.0; // [W] = [J/s]
  for (const auto &s : streams) {
    double h_molar = s.h();                  // J/mol
    double MW_i = s.mw();                    // g/mol
    double h_mass = h_molar / (MW_i * 1e-3); // J/kg
    H_dot_total += h_mass * s.mdot;
  }

  double h_out_mass = H_dot_total / mdot_total; // J/kg

  // Convert to molar enthalpy for T solver
  double MW_out = mwmix(X_out);                      // g/mol
  double h_out_molar = h_out_mass * (MW_out * 1e-3); // J/mol

  // Solve for T_out: h(T_out, X_out) = h_out_molar
  // Use mass-weighted average T as initial guess
  double T_guess = 0.0;
  for (const auto &s : streams) {
    T_guess += s.T() * s.mdot;
  }
  T_guess /= mdot_total;

  double T_out = calc_T_from_h(h_out_molar, X_out, T_guess);

  // Build output stream
  Stream out;
  out.state.T = T_out;
  out.state.P = P_out;
  out.state.X = X_out;
  out.mdot = mdot_total;

  return out;
}
