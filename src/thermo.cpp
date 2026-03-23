#include "../include/thermo.h"
#include "../include/composition.h"
#include "../include/humidair.h"
#include "../include/thermo_transport_data.h"
#include "../include/transport.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace combaero {

double J_per_mol_to_J_per_kg(double value, double molar_mass) {
  return value * 1000.0 /
         molar_mass; // molar_mass is in g/mol, convert to kg/mol
}

// -------------------------------------------------------------
// Species metadata and lookup functions
// -------------------------------------------------------------

std::string species_name(std::size_t species_index) {
    if (species_index >= species_names.size()) {
        throw std::out_of_range("species_name: index out of bounds");
    }
    return species_names[species_index];
}

std::size_t species_index_from_name(const std::string& name) {
    auto it = species_index.find(name);
    if (it == species_index.end()) {
        throw std::invalid_argument("species_index_from_name: species not found: " + name);
    }
    return static_cast<std::size_t>(it->second);
}

std::size_t num_species() {
    return species_names.size();
}

double species_molar_mass(std::size_t species_index) {
    if (species_index >= molar_masses.size()) {
        throw std::out_of_range("species_molar_mass: index out of bounds");
    }
    return molar_masses[species_index];
}

double species_molar_mass_from_name(const std::string& name) {
    return species_molar_mass(species_index_from_name(name));
}


// -------------------------------------------------------------
// NASA 9 Polynomial evaluation
// -------------------------------------------------------------

namespace {
const NASA9_Interval& get_interval(std::size_t species_idx, double T) {
    const auto& nasa = nasa_coeffs[species_idx];
    for (const auto& interval : nasa.intervals) {
        if (T >= interval.T_min && T <= interval.T_max) {
            return interval;
        }
    }
    // Fallback to nearest if out of range
    if (T < nasa.intervals.front().T_min) return nasa.intervals.front();
    return nasa.intervals.back();
}
}

double cp_R(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("cp_R: species index out of bounds");
  }
  const auto& interval = get_interval(species_idx, T);
  const auto& a = interval.coeffs;
  // Cp/R = a1/T^2 + a2/T + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
  return a[0]/(T*T) + a[1]/T + a[2] + a[3]*T + a[4]*T*T + a[5]*T*T*T + a[6]*T*T*T*T;
}

double h_RT(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("h_RT: species index out of bounds");
  }
  const auto& interval = get_interval(species_idx, T);
  const auto& a = interval.coeffs;
  // H/RT = -a1/T^2 + a2*ln(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + b1/T
  return -a[0]/(T*T) + a[1]*std::log(T)/T + a[2] + a[3]*T/2.0 + a[4]*T*T/3.0 +
         a[5]*T*T*T/4.0 + a[6]*T*T*T*T/5.0 + a[8]/T;
}

double s_R(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("s_R: species index out of bounds");
  }
  const auto& interval = get_interval(species_idx, T);
  const auto& a = interval.coeffs;
  // S/R = -a1/(2*T^2) - a2/T + a3*ln(T) + a4*T + a5*T^2/2 + a6*T^3/3 + a7*T^4/4 + a9
  return -a[0]/(2.0*T*T) - a[1]/T + a[2]*std::log(T) + a[3]*T + a[4]*T*T/2.0 +
         a[5]*T*T*T/3.0 + a[6]*T*T*T*T/4.0 + a[9];
}

double g_over_RT(std::size_t species_idx, double T) {
  return h_RT(species_idx, T) - s_R(species_idx, T);
}

// -------------------------------------------------------------
// Dimensional property evaluation
// -------------------------------------------------------------

double cp_species(std::size_t species_idx, double T) {
  return cp_R(species_idx, T) * thermo::R_GAS;
}

double h_species(std::size_t species_idx, double T) {
  return h_RT(species_idx, T) * thermo::R_GAS * T;
}

double s_species(std::size_t species_idx, double T) {
  return s_R(species_idx, T) * thermo::R_GAS;
}

// -------------------------------------------------------------
// Mixture property evaluation
// -------------------------------------------------------------

double cp(double T, const std::vector<double> &X) {
  double cp_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    cp_mix += X[i] * cp_species(i, T);
  }
  return cp_mix;
}

double h(double T, const std::vector<double> &X) {
  double h_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    h_mix += X[i] * h_species(i, T);
  }
  return h_mix;
}

double s(double T, const std::vector<double> &X, double P, double P_ref) {
  double s_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    if (X[i] > 1e-15) {
      s_mix += X[i] * (s_species(i, T) - thermo::R_GAS * std::log(X[i] * P / P_ref));
    }
  }
  return s_mix;
}

double cv(double T, const std::vector<double> &X) {
  return cp(T, X) - thermo::R_GAS;
}

double u(double T, const std::vector<double> &X) {
  return h(T, X) - thermo::R_GAS * T;
}

double density(double T, double P, const std::vector<double> &X) {
  return (P * mwmix(X)) / (thermo::R_GAS * 1000.0 * T);
}

double molar_volume(double T, double P) {
  return (thermo::R_GAS * T) / P;
}

double specific_gas_constant(const std::vector<double> &X) {
  return (thermo::R_GAS * 1000.0) / mwmix(X);
}

double isentropic_expansion_coefficient(double T, const std::vector<double> &X) {
  return cp(T, X) / cv(T, X);
}

double speed_of_sound(double T, const std::vector<double> &X) {
  double gamma = isentropic_expansion_coefficient(T, X);
  double R_spec = specific_gas_constant(X);
  return std::sqrt(gamma * R_spec * T);
}

// -------------------------------------------------------------
// Mass-specific thermodynamic properties
// -------------------------------------------------------------

double cp_mass(double T, const std::vector<double> &X) {
  return cp(T, X) / (mwmix(X) / 1000.0);
}

double cv_mass(double T, const std::vector<double> &X) {
  return cv(T, X) / (mwmix(X) / 1000.0);
}

double h_mass(double T, const std::vector<double> &X) {
  return h(T, X) / (mwmix(X) / 1000.0);
}

double s_mass(double T, const std::vector<double> &X, double P, double P_ref) {
  return s(T, X, P, P_ref) / (mwmix(X) / 1000.0);
}

double u_mass(double T, const std::vector<double> &X) {
  return u(T, X) / (mwmix(X) / 1000.0);
}

// -------------------------------------------------------------
// Air Properties Bundle
// -------------------------------------------------------------

AirProperties air_properties(double T, double P, double humidity) {
  if (T <= 0) {
    throw std::invalid_argument("air_properties: temperature must be positive");
  }
  if (P <= 0) {
    throw std::invalid_argument("air_properties: pressure must be positive");
  }
  if (humidity < 0.0 || humidity > 1.0) {
    throw std::invalid_argument(
        "air_properties: humidity must be in range [0, 1]");
  }

  std::vector<double> X;
  if (humidity > 0.0) {
    X = humid_air_composition(T, P, humidity);
  } else {
    X = dry_air();
  }
  return transport_state(T, P, X);
}

// -------------------------------------------------------------
// Thermodynamic State Bundle
// -------------------------------------------------------------

ThermoState thermo_state(double T, double P, const std::vector<double> &X, double P_ref) {
  if (T <= 0.0) throw std::invalid_argument("thermo_state: temperature must be positive");
  if (P <= 0.0) throw std::invalid_argument("thermo_state: pressure must be positive");
  if (X.empty()) throw std::invalid_argument("thermo_state: composition must not be empty");
  if (P_ref <= 0.0) throw std::invalid_argument("thermo_state: reference pressure must be positive");
  ThermoState ts;
  ts.T = T;
  ts.P = P;
  ts.mw = mwmix(X);
  double mw_kg = ts.mw / 1000.0;
  ts.cp_mole = cp(T, X);
  ts.h_mole = h(T, X);
  ts.s_mole = s(T, X, P, P_ref);
  ts.cv_mole = cv(T, X);
  ts.u_mole = u(T, X);
  ts.rho = density(T, P, X);
  ts.gamma = ts.cp_mole / ts.cv_mole;
  ts.a = std::sqrt(ts.gamma * (thermo::R_GAS / mw_kg) * T);
  ts.cp = ts.cp_mole / mw_kg;
  ts.cv = ts.cv_mole / mw_kg;
  ts.h = ts.h_mole / mw_kg;
  ts.s = ts.s_mole / mw_kg;
  ts.u = ts.u_mole / mw_kg;
  return ts;
}

// -------------------------------------------------------------
// Complete State Bundle
// -------------------------------------------------------------

CompleteState complete_state(double T, double P, const std::vector<double> &X, double P_ref) {
  CompleteState cs;
  cs.thermo = thermo_state(T, P, X, P_ref);
  cs.transport = transport_state(T, P, X);
  return cs;
}

// -------------------------------------------------------------
// Derivatives wrt temperature
// -------------------------------------------------------------

double dh_dT(double T, const std::vector<double> &X) { return cp(T, X); }

double ds_dT(double T, const std::vector<double> &X) { return cp(T, X) / T; }

double dcp_dT(double T, const std::vector<double> &X) {
  double dcp_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    const auto& interval = get_interval(i, T);
    const auto& a = interval.coeffs;
    // d(Cp/R)/dT = -2*a1/T^3 - a2/T^2 + a4 + 2*a5*T + 3*a6*T^2 + 4*a7*T^3
    double dcp_R_dT = -2.0*a[0]/(T*T*T) - a[1]/(T*T) + a[3] + 2.0*a[4]*T +
                      3.0*a[5]*T*T + 4.0*a[6]*T*T*T;
    dcp_mix += X[i] * dcp_R_dT * thermo::R_GAS;
  }
  return dcp_mix;
}

double dg_over_RT_dT(double T, const std::vector<double> &X) {
  return -h(T, X) / (thermo::R_GAS * T * T);
}

// -------------------------------------------------------------
// Inverse solvers (T from properties)
// -------------------------------------------------------------

double calc_T_from_h(double h_target, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  double T = T_guess;
  for (std::size_t i = 0; i < max_iter; ++i) {
    double f = h(T, X) - h_target;
    double df = cp(T, X);
    double dT = f / df;
    T -= dT;
    if (std::abs(dT) < tol)
      return T;
  }
  return T;
}

double calc_T_from_s(double s_target, double P, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  double T = T_guess;
  for (std::size_t i = 0; i < max_iter; ++i) {
    double f = s(T, X, P) - s_target;
    double df = cp(T, X) / T;
    double dT = f / df;
    T -= dT;
    if (std::abs(dT) < tol)
      return T;
  }
  return T;
}

double calc_T_from_cp(double cp_target, const std::vector<double> &X,
                      double T_guess, double tol, std::size_t max_iter) {
  double T = T_guess;
  for (std::size_t i = 0; i < max_iter; ++i) {
    double cp_val = cp(T, X);
    double dcp_dT_val = dcp_dT(T, X);
    if (std::abs(dcp_dT_val) < 1.0e-10) {
      T *= (cp_val < cp_target) ? 1.1 : 0.9;
      if (T < 200.0) T = 200.0;
      if (T > 6000.0) T = 6000.0;
      continue;
    }
    double delta_T = (cp_target - cp_val) / dcp_dT_val;
    double max_step = 0.1 * T;
    if (std::abs(delta_T) > max_step)
      delta_T = (delta_T > 0) ? max_step : -max_step;
    T += delta_T;
    if (T < 200.0) T = 200.0;
    if (T > 6000.0) T = 6000.0;
    if (std::abs(delta_T) < tol || std::abs(cp_val - cp_target) < tol)
      return T;
  }
  return T;
}

double calc_T_from_u(double u_target, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  double T = T_guess;
  for (std::size_t i = 0; i < max_iter; ++i) {
    double f = u(T, X) - u_target;
    double df = cv(T, X);
    double dT = f / df;
    T -= dT;
    if (std::abs(dT) < tol)
      return T;
  }
  return T;
}

double calc_T_from_h_mass(double h_mass_target, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  return calc_T_from_h(h_mass_target * (mwmix(X) / 1000.0), X, T_guess, tol, max_iter);
}

double calc_T_from_s_mass(double s_mass_target, double P, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  return calc_T_from_s(s_mass_target * (mwmix(X) / 1000.0), P, X, T_guess, tol, max_iter);
}

double calc_T_from_u_mass(double u_mass_target, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  return calc_T_from_u(u_mass_target * (mwmix(X) / 1000.0), X, T_guess, tol, max_iter);
}

double calc_T_from_sv_mass(double s_mass_target, double v_mass_target, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
  double T = T_guess;
  double mw_kg = mwmix(X) / 1000.0;
  for (std::size_t i = 0; i < max_iter; ++i) {
    double rho = 1.0 / v_mass_target;
    double P = (rho * thermo::R_GAS * T) / mw_kg;
    double f = s_mass(T, X, P) - s_mass_target;
    // ds/dT = cp/T - R/P * dP/dT = cp/T - R/P * (rho R / mw) = cp/T - R/P * (P/T) = cp/T - R/T = cv/T
    double df = cv_mass(T, X) / T;
    double dT = f / df;
    T -= dT;
    if (std::abs(dT) < tol) return T;
  }
  return T;
}

double calc_T_from_sh_mass([[maybe_unused]] double s_mass_target, double h_mass_target, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
  // h = h(T) only for ideal gas. Find T from h first.
  double T = calc_T_from_h_mass(h_mass_target, X, T_guess, tol, max_iter);
  // Then P is determined by s(T, P) = s_target
  return T;
}

// -------------------------------------------------------------
// State-based overloads
// -------------------------------------------------------------

double mwmix(const State &s) { return combaero::mwmix(s.X); }
double cp(const State &s) { return combaero::cp(s.T, s.X); }
double h(const State &s) { return combaero::h(s.T, s.X); }
double s(const State &s) { return combaero::s(s.T, s.X, s.P); }
double cv(const State &s) { return combaero::cv(s.T, s.X); }
double u(const State &s) { return combaero::u(s.T, s.X); }
double density(const State &s) { return combaero::density(s.T, s.P, s.X); }
double specific_gas_constant(const State &s) { return combaero::specific_gas_constant(s.X); }
double isentropic_expansion_coefficient(const State &s) { return combaero::isentropic_expansion_coefficient(s.T, s.X); }
double speed_of_sound(const State &s) { return combaero::speed_of_sound(s.T, s.X); }

} // namespace combaero
