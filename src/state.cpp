#include "../include/state.h"
#include "../include/thermo.h"
#include "../include/transport.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <cstdio>

namespace combaero {

// -------------------------------------------------------------
// State property getters
// -------------------------------------------------------------

double State::mw() const { return mwmix(X); }
double State::cp() const { return cp_mass(T, X); }
double State::h() const { return h_mass(T, X); }
double State::s() const { return s_mass(T, X, P); }
double State::cv() const { return cv_mass(T, X); }
double State::u() const { return u_mass(T, X); }
double State::rho() const { return density(T, P, X); }
double State::R() const { return 8.314462618 / (mw() / 1000.0); }
double State::gamma() const { return cp() / cv(); }
double State::a() const { return speed_of_sound(T, X); }

double State::mu() const { return viscosity(T, P, X); }
double State::k() const { return thermal_conductivity(T, P, X); }
double State::nu() const { return mu() / rho(); }
double State::Pr() const { return mu() * (cp()) / k(); }
double State::alpha() const { return k() / (rho() * cp()); }

// -------------------------------------------------------------
// Joint Setters
// -------------------------------------------------------------

State &State::set_TPX(double T_new, double P_new, const std::vector<double> &X_new) {
  T = T_new;
  P = P_new;
  set_X(X_new);
  return *this;
}

State &State::set_TPY(double T_new, double P_new, const std::vector<double> &Y_new) {
  T = T_new;
  P = P_new;
  set_Y(Y_new);
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
  // rho = P / (R * T) => T = P / (rho * R)
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
  // h = u + Pv = u + RT => find T such that u(T) + RT = h_target
  // This is slightly more complex, but for ideal gas we can iterate
  double T_guess = T;
  for (int i = 0; i < 20; ++i) {
    double f = u_mass(T_guess, X) + R() * T_guess - h_new;
    double df = cv_mass(T_guess, X) + R();
    double dT = f / df;
    T_guess -= dT;
    if (std::abs(dT) < 1e-4)
      break;
  }
  T = T_guess;
  P = R() * T / v_new;
  return *this;
}

State &State::set_SH_mass(double s_new, double h_new) {
  T = calc_T_from_sh_mass(s_new, h_new, X, T);
  // Iteratively find P that matches s(T, P) = s_new
  // s(T, P) = s0(T) - R ln(P/P_ref)
  // s_new - s0(T) = -R ln(P/P_ref) => ln(P/P_ref) = (s0(T) - s_new)/R
  // P = P_ref * exp((s0(T) - s_new)/R)
  double s0 = s_mass(T, X, 101325.0); // entropy at P_ref
  P = 101325.0 * std::exp((s0 - s_new) / R());
  return *this;
}

// -------------------------------------------------------------
// Joint Getters
// -------------------------------------------------------------

std::tuple<double, double, std::vector<double>> State::TPX() const {
  return {T, P, X};
}
std::tuple<double, double, std::vector<double>> State::TPY() const {
  return {T, P, Y};
}
std::tuple<double, double> State::TP() const { return {T, P}; }
std::tuple<double, double> State::DP_mass() const { return {rho(), P}; }
std::tuple<double, double> State::HP_mass() const { return {h(), P}; }
std::tuple<double, double> State::SP_mass() const { return {s(), P}; }
std::tuple<double, double> State::UV_mass() const { return {u(), 1.0 / rho()}; }
std::tuple<double, double> State::SV_mass() const { return {s(), 1.0 / rho()}; }
std::tuple<double, double> State::PV_mass() const { return {P, 1.0 / rho()}; }
std::tuple<double, double> State::UP_mass() const { return {u(), P}; }
std::tuple<double, double> State::VH_mass() const {
  return {1.0 / rho(), h()};
}
std::tuple<double, double> State::SH_mass() const { return {s(), h()}; }

// -------------------------------------------------------------
// I/O
// -------------------------------------------------------------

void State::print() const {
  std::printf("State: T=%.2f K, P=%.1f Pa, rho=%.3f kg/m3\n", T, P, rho());
}

// -------------------------------------------------------------
// Stream implementation
// -------------------------------------------------------------

Stream mix(const std::vector<Stream> &streams, double P_out, double delta_h) {
  if (streams.empty()) {
    return Stream();
  }

  double mdot_tot = 0.0;
  double H_tot = 0.0;
  std::size_t n_spec = streams[0].state.X.size();
  std::vector<double> Y_mix(n_spec, 0.0);
  double P_min = streams[0].state.P;

  for (const auto &s : streams) {
    mdot_tot += s.mdot;
    H_tot += s.mdot * s.state.h();
    if (s.state.P < P_min)
      P_min = s.state.P;
    for (std::size_t k = 0; k < n_spec; ++k) {
      Y_mix[k] += s.mdot * s.state.Y[k];
    }
  }

  Stream res;
  res.mdot = mdot_tot;
  if (mdot_tot > 1e-12) {
    for (std::size_t k = 0; k < n_spec; ++k) {
      Y_mix[k] /= mdot_tot;
    }
    double h_mix = H_tot / mdot_tot + delta_h;
    double P_res = (P_out < 0.0) ? P_min : P_out;
    res.state.set_Y(Y_mix);
    res.state.set_HP_mass(h_mix, P_res);
  } else {
    res.state = streams[0].state;
  }
  return res;
}

} // namespace combaero
