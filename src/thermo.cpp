#include "../include/thermo.h"
#include "../include/composition.h"
#include "../include/humidair.h"
#include "../include/thermo_c2_refit_data.h"
#include "../include/thermo_transport_data.h"
#include "../include/transport.h"
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
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

std::string species_name(std::size_t idx) {
    if (idx >= species_names.size()) {
        throw std::out_of_range("species_name: index out of bounds");
    }
    return species_names[idx];
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

double species_molar_mass(std::size_t idx) {
    if (idx >= molar_masses.size()) {
        throw std::out_of_range("species_molar_mass: index out of bounds");
    }
    return molar_masses[idx];
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

// --- the C2 re-fit (#549) ----------------------------------------------------

std::atomic<int> g_thermo_model{static_cast<int>(ThermoModel::NASA9)};

// The bridge covering T for this species, or nullptr (NASA9 model, or T
// outside every window).
const C2Bridge* active_bridge(std::size_t species_idx, double T) {
  if (g_thermo_model.load(std::memory_order_relaxed) !=
      static_cast<int>(ThermoModel::NASA9_C2)) {
    return nullptr;
  }
  for (const auto& b : c2_bridges[species_idx]) {
    if (std::abs(T - b.T_join) < b.half_width) {
      return &b;
    }
  }
  return nullptr;
}

double bridge_poly(const C2Bridge& b, double t) {
  double acc = 0.0;
  for (std::size_t k = b.c.size(); k-- > 0;) {
    acc = acc * t + b.c[k];
  }
  return acc;
}

double bridge_poly_dt(const C2Bridge& b, double t) {
  double acc = 0.0;
  for (std::size_t k = b.c.size(); k-- > 1;) {
    acc = acc * t + static_cast<double>(k) * b.c[k];
  }
  return acc;
}

// Antiderivative in t of the bridge polynomial.
double bridge_poly_int(const C2Bridge& b, double t) {
  double acc = 0.0;
  for (std::size_t k = b.c.size(); k-- > 0;) {
    acc = acc * t + b.c[k] / static_cast<double>(k + 1);
  }
  return acc * t;
}

// The NASA-9 interval that ends at the bridge's join (its left side).
const NASA9_Interval& left_interval(std::size_t species_idx, double T_join) {
  for (const auto& iv : nasa_coeffs[species_idx].intervals) {
    if (iv.T_max == T_join) {
      return iv;
    }
  }
  return nasa_coeffs[species_idx].intervals.front();
}

// NASA-9 H/R [K] and S/R of one interval (the published closed forms).
double nasa_H_R(const std::array<double, 10>& a, double T) {
  return -a[0] / T + a[1] * std::log(T) + a[2] * T + a[3] * T * T / 2.0 +
         a[4] * T * T * T / 3.0 + a[5] * T * T * T * T / 4.0 +
         a[6] * T * T * T * T * T / 5.0 + a[8];
}

double nasa_S_R(const std::array<double, 10>& a, double T) {
  return -a[0] / (2.0 * T * T) - a[1] / T + a[2] * std::log(T) + a[3] * T +
         a[4] * T * T / 2.0 + a[5] * T * T * T / 3.0 +
         a[6] * T * T * T * T / 4.0 + a[9];
}

// H/R [K] inside a bridge: NASA-9 at the window's left edge plus the
// bridge's integral of cp/R dT, closed form.
double bridge_H_R(std::size_t species_idx, const C2Bridge& b, double T) {
  const double w = b.half_width;
  const double t = (T - b.T_join) / w;
  const auto& L = left_interval(species_idx, b.T_join).coeffs;
  return nasa_H_R(L, b.T_join - w) +
         w * (bridge_poly_int(b, t) - bridge_poly_int(b, -1.0));
}

// S/R inside a bridge: NASA-9 at the left edge plus the integral of cp/(R T)
// dT over [-1, t], by 8-point Gauss-Legendre. The integrand is a polynomial
// over (T_join/w + t), analytic with its pole at t = -20: the rule is exact
// to rounding.
double bridge_S_R(std::size_t species_idx, const C2Bridge& b, double T) {
  static const double xg[4] = {0.1834346424956498, 0.5255324099163290,
                               0.7966664774136267, 0.9602898564975363};
  static const double wg[4] = {0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763};
  const double w = b.half_width;
  const double t = (T - b.T_join) / w;
  const double c = 0.5 * (t - 1.0), r = 0.5 * (t + 1.0);
  auto f = [&](double tau) { return bridge_poly(b, tau) / (b.T_join + w * tau); };
  double acc = 0.0;
  for (int i = 0; i < 4; ++i) {
    acc += wg[i] * (f(c - r * xg[i]) + f(c + r * xg[i]));
  }
  const auto& L = left_interval(species_idx, b.T_join).coeffs;
  return nasa_S_R(L, b.T_join - w) + w * r * acc;
}
}  // namespace

void set_thermo_model(ThermoModel model) {
  g_thermo_model.store(static_cast<int>(model), std::memory_order_relaxed);
}

ThermoModel thermo_model() {
  return static_cast<ThermoModel>(g_thermo_model.load(std::memory_order_relaxed));
}

double cp_R(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("cp_R: species index out of bounds");
  }
  double T_clamped = std::max(T, 10.0); const auto& interval = get_interval(species_idx, T_clamped); T = T_clamped;
  if (const C2Bridge* b = active_bridge(species_idx, T)) {
    return bridge_poly(*b, (T - b->T_join) / b->half_width);
  }
  const auto& a = interval.coeffs;
  // Cp/R = a1/T^2 + a2/T + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
  return a[0]/(T*T) + a[1]/T + a[2] + a[3]*T + a[4]*T*T + a[5]*T*T*T + a[6]*T*T*T*T;
}

double h_RT(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("h_RT: species index out of bounds");
  }
  double T_clamped = std::max(T, 10.0); const auto& interval = get_interval(species_idx, T_clamped); T = T_clamped;
  if (const C2Bridge* b = active_bridge(species_idx, T)) {
    return bridge_H_R(species_idx, *b, T) / T;
  }
  const auto& a = interval.coeffs;
  // H/RT = -a1/T^2 + a2*ln(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + b1/T
  return -a[0]/(T*T) + a[1]*std::log(T)/T + a[2] + a[3]*T/2.0 + a[4]*T*T/3.0 +
         a[5]*T*T*T/4.0 + a[6]*T*T*T*T/5.0 + a[8]/T;
}

double s_R(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("s_R: species index out of bounds");
  }
  double T_clamped = std::max(T, 10.0); const auto& interval = get_interval(species_idx, T_clamped); T = T_clamped;
  if (const C2Bridge* b = active_bridge(species_idx, T)) {
    return bridge_S_R(species_idx, *b, T);
  }
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

double dcp_R_dT(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("dcp_R_dT: species index out of bounds");
  }
  if (T < 10.0) {
    return 0.0;  // cp_R holds its 10 K value below the clamp
  }
  if (const C2Bridge* b = active_bridge(species_idx, T)) {
    return bridge_poly_dt(*b, (T - b->T_join) / b->half_width) / b->half_width;
  }
  const auto& a = get_interval(species_idx, T).coeffs;
  // d(Cp/R)/dT = -2*a1/T^3 - a2/T^2 + a4 + 2*a5*T + 3*a6*T^2 + 4*a7*T^3
  return -2.0*a[0]/(T*T*T) - a[1]/(T*T) + a[3] + 2.0*a[4]*T +
         3.0*a[5]*T*T + 4.0*a[6]*T*T*T;
}

double dcp_dT(double T, const std::vector<double> &X) {
  double dcp_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    if (const C2Bridge* b = active_bridge(i, T)) {
      dcp_mix += X[i] * thermo::R_GAS *
                 bridge_poly_dt(*b, (T - b->T_join) / b->half_width) / b->half_width;
      continue;
    }
    const auto& a = get_interval(i, T).coeffs;
    double dcp_R = -2.0*a[0]/(T*T*T) - a[1]/(T*T) + a[3] + 2.0*a[4]*T +
                   3.0*a[5]*T*T + 4.0*a[6]*T*T*T;
    dcp_mix += X[i] * dcp_R * thermo::R_GAS;
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
