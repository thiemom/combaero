#include "liner_radiation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "dual_number.h"
#include "thermo_transport_data.h"

namespace combaero {

namespace {

// Gas radiation's 0.5 (1 + eps_w) wall-absorptivity factor (Eq. 8.5).
constexpr double kHalf = 0.5;
// Eq. 8.6 exponent on Tg/Tw, and the 2.5 = 4 - 1.5 it leaves in Eq. 8.7.
constexpr double kAbsorptivityExponent = 1.5;
constexpr double kR1Exponent = 2.5;

// Surface-temperature fixed point (wall_radiation_sides): Newton tolerance
// relative to the temperature, iteration cap, and the positive floor that
// keeps a Newton trial inside sqrt's domain (exact above it).
constexpr double kWallRadRelTol = 1e-13;
constexpr int kWallRadMaxIter = 100;
constexpr double kWallRadTFloor = 1e-3;

void require_finite(double v, const char *what) {
  if (!std::isfinite(v)) {
    throw std::invalid_argument(std::string(what) + " must be finite");
  }
}

void require_emissivity(double eps, const char *what) {
  if (!(eps > 0.0 && eps <= 1.0)) {
    throw std::invalid_argument(std::string(what) +
                                ": emissivity must be in (0, 1]");
  }
}

using solver::DualN;

// h_r of one face (see the header), as a DualN of whatever is seeded.
template <int N>
DualN<N> radiative_htc(int kind, const DualN<N> &Ts, const DualN<N> &Tsrc,
                       const DualN<N> &k) {
  if (kind == static_cast<int>(WallRadiationKind::Gas)) {
    const DualN<N> s = solver::dsqrt(solver::dpos_floor(Tsrc, kWallRadTFloor));
    const DualN<N> t = solver::dsqrt(solver::dpos_floor(Ts, kWallRadTFloor));
    // (s^5 - t^5) / (s^2 - t^2) = (s^4 + s^3 t + s^2 t^2 + s t^3 + t^4)/(s + t)
    const DualN<N> s2 = s * s;
    const DualN<N> t2 = t * t;
    const DualN<N> num = s2 * s2 + s2 * s * t + s2 * t2 + s * t2 * t + t2 * t2;
    const DualN<N> d25 = num / (s + t);
    return k * (s2 * s) * d25;
  }
  if (kind == static_cast<int>(WallRadiationKind::Surface)) {
    return k * (Ts + Tsrc) * (Ts * Ts + Tsrc * Tsrc);
  }
  return DualN<N>::constant(0.0);
}

template <int N> struct WallRadEval {
  DualN<N> H_a, T_a, H_b, T_b, hr_a, hr_b, q, Ts_a_new, Ts_b_new;
};

// One evaluation of the effective sides and the profile at given surface
// temperatures: p = [T_aw_a, h_a, T_aw_b, h_b, R_wall, T_src_a, k_a, T_src_b,
// k_b].
template <int N>
WallRadEval<N> wall_rad_eval(const DualN<N> &Ts_a, const DualN<N> &Ts_b,
                             const std::array<DualN<N>, WALL_RAD_N_PARAMS> &p,
                             int kind_a, int kind_b) {
  WallRadEval<N> e;
  e.hr_a = radiative_htc<N>(kind_a, Ts_a, p[5], p[6]);
  e.hr_b = radiative_htc<N>(kind_b, Ts_b, p[7], p[8]);
  e.H_a = p[1] + e.hr_a;
  e.H_b = p[3] + e.hr_b;
  e.T_a = (p[1] * p[0] + e.hr_a * p[5]) / e.H_a;
  e.T_b = (p[3] * p[2] + e.hr_b * p[7]) / e.H_b;
  const DualN<N> R_tot = 1.0 / e.H_a + p[4] + 1.0 / e.H_b;
  e.q = (e.T_a - e.T_b) / R_tot;
  e.Ts_a_new = e.T_a - e.q / e.H_a;
  e.Ts_b_new = e.T_b + e.q / e.H_b;
  return e;
}

} // namespace

GasEmissivityResult gas_emissivity_luminous(double T_g, double P, double q,
                                            double l_b, double L) {
  require_finite(T_g, "T_g");
  require_finite(P, "P");
  require_finite(q, "q");
  require_finite(l_b, "l_b");
  require_finite(L, "L");
  if (!(T_g > 0.0)) {
    throw std::invalid_argument("gas_emissivity_luminous: T_g must be > 0");
  }
  if (P < 0.0 || q < 0.0 || l_b < 0.0 || L < 0.0) {
    throw std::invalid_argument(
        "gas_emissivity_luminous: P, q, l_b and L must be >= 0");
  }
  GasEmissivityResult r;
  const double s = std::sqrt(q * l_b);
  // E = base P L s, base = 290 / 1000 Tg^-1.5 (P in Pa -> kPa).
  const double base = REEVES_EMISSIVITY_COEFF / LINER_RADIATION_PA_PER_KPA *
                      std::pow(T_g, -kAbsorptivityExponent);
  const double E = base * P * L * s;
  const double ex = std::exp(-E);
  r.eps = 1.0 - ex;
  r.deps_dT = ex * (-kAbsorptivityExponent * E / T_g);
  r.deps_dP = ex * base * L * s;
  r.deps_dL = ex * base * P * s;
  if (s > 0.0) {
    r.deps_dq = ex * base * P * L * 0.5 * l_b / s;
    r.deps_dl_b = ex * base * P * L * 0.5 * q / s;
  }
  return r;
}

std::tuple<double, double> luminosity_factor_hydrogen(double H_mass_percent) {
  require_finite(H_mass_percent, "H_mass_percent");
  if (!(H_mass_percent > 0.0)) {
    throw std::invalid_argument(
        "luminosity_factor_hydrogen: H_mass_percent must be > 0");
  }
  const double H = H_mass_percent;
  return {LUMINOSITY_HYDROGEN_COEFF / (H * H),
          -2.0 * LUMINOSITY_HYDROGEN_COEFF / (H * H * H)};
}

LinerGasRadiationResult liner_gas_radiation(double T_g, double T_w,
                                            double eps_g, double eps_w) {
  require_finite(T_g, "T_g");
  require_finite(T_w, "T_w");
  require_finite(eps_g, "eps_g");
  if (!(T_g > 0.0) || T_w < 0.0) {
    throw std::invalid_argument(
        "liner_gas_radiation: T_g must be > 0 and T_w >= 0");
  }
  if (eps_g < 0.0 || eps_g > 1.0) {
    throw std::invalid_argument("liner_gas_radiation: eps_g must be in [0, 1]");
  }
  require_emissivity(eps_w, "liner_gas_radiation eps_w");
  const double c = kHalf * STEFAN_BOLTZMANN;
  const double Tg15 = std::pow(T_g, kAbsorptivityExponent);
  const double Tw15 = std::pow(T_w, kAbsorptivityExponent);
  const double Tg25 = Tg15 * T_g;
  const double Tw25 = Tw15 * T_w;
  const double diff = Tg25 - Tw25;
  LinerGasRadiationResult r;
  const double pre = c * (1.0 + eps_w) * eps_g;
  r.R1 = pre * Tg15 * diff;
  r.dR1_dT_g = pre * (kAbsorptivityExponent * std::sqrt(T_g) * diff +
                      kR1Exponent * Tg15 * Tg15);
  r.dR1_dT_w = -pre * kR1Exponent * Tg15 * Tw15;
  r.dR1_deps_g = c * (1.0 + eps_w) * Tg15 * diff;
  r.dR1_deps_w = c * eps_g * Tg15 * diff;
  return r;
}

double liner_casing_exchange_factor(double eps_w, double eps_c,
                                    double area_ratio) {
  require_emissivity(eps_w, "liner_casing_radiation eps_w");
  require_emissivity(eps_c, "liner_casing_radiation eps_c");
  require_finite(area_ratio, "area_ratio");
  if (!(area_ratio > 0.0)) {
    throw std::invalid_argument(
        "liner_casing_radiation: area_ratio (A_w/A_c) must be > 0");
  }
  return STEFAN_BOLTZMANN * eps_w * eps_c /
         (eps_c + eps_w * (1.0 - eps_c) * area_ratio);
}

std::tuple<double, double, double>
liner_casing_radiation(double T_w, double T_c, double eps_w, double eps_c,
                       double area_ratio) {
  require_finite(T_w, "T_w");
  require_finite(T_c, "T_c");
  if (T_w < 0.0 || T_c < 0.0) {
    throw std::invalid_argument(
        "liner_casing_radiation: temperatures must be >= 0");
  }
  const double F = liner_casing_exchange_factor(eps_w, eps_c, area_ratio);
  const double Tw2 = T_w * T_w;
  const double Tc2 = T_c * T_c;
  return {F * (Tw2 * Tw2 - Tc2 * Tc2), 4.0 * F * Tw2 * T_w,
          -4.0 * F * Tc2 * T_c};
}

double mean_beam_length(double volume, double area, double coeff) {
  require_finite(volume, "volume");
  require_finite(area, "area");
  require_finite(coeff, "coeff");
  if (volume < 0.0 || !(area > 0.0) || !(coeff > 0.0)) {
    throw std::invalid_argument(
        "mean_beam_length: needs volume >= 0, area > 0 and coeff > 0");
  }
  return coeff * volume / area;
}

std::tuple<double, std::vector<double>>
elemental_fuel_air_ratio(const std::vector<double> &Y) {
  const std::size_t n = molar_masses.size();
  if (Y.size() != n) {
    throw std::invalid_argument(
        "elemental_fuel_air_ratio: Y needs one mass fraction per species");
  }
  std::vector<double> a(n);
  double w = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    const Molecular_Structure &ms = molecular_structures[k];
    a[k] = (static_cast<double>(ms.C) * ATOMIC_MASS_C +
            static_cast<double>(ms.H) * ATOMIC_MASS_H) /
           molar_masses[k];
    w += Y[k] * a[k];
  }
  if (!(w < 1.0)) {
    throw std::invalid_argument(
        "elemental_fuel_air_ratio: the mixture holds no air (all C and H)");
  }
  const double one_m = 1.0 - w;
  std::vector<double> dq(n);
  for (std::size_t k = 0; k < n; ++k) {
    dq[k] = a[k] / (one_m * one_m);
  }
  return {w / one_m, dq};
}

WallRadiationSides wall_radiation_sides(double T_aw_a, double h_a,
                                        double T_aw_b, double h_b,
                                        double R_wall, int kind_a,
                                        double T_src_a, double k_a, int kind_b,
                                        double T_src_b, double k_b) {
  const std::array<double, WALL_RAD_N_PARAMS> pv{
      T_aw_a, h_a, T_aw_b, h_b, R_wall, T_src_a, k_a, T_src_b, k_b};
  for (double v : pv) {
    require_finite(v, "wall_radiation_sides input");
  }
  if (!(h_a > 0.0) || !(h_b > 0.0) || R_wall < 0.0 || k_a < 0.0 ||
      k_b < 0.0) {
    throw std::invalid_argument(
        "wall_radiation_sides: needs h > 0, R_wall >= 0 and k >= 0");
  }
  for (int kind : {kind_a, kind_b}) {
    if (kind < 0 || kind > 2) {
      throw std::invalid_argument("wall_radiation_sides: kind must be 0, 1 or 2");
    }
  }
  if ((kind_a == 1 && !(T_src_a > 0.0)) || (kind_b == 1 && !(T_src_b > 0.0)) ||
      (kind_a == 2 && T_src_a < 0.0) || (kind_b == 2 && T_src_b < 0.0)) {
    throw std::invalid_argument(
        "wall_radiation_sides: a gas temperature must be > 0, a surface's >= 0");
  }

  // Newton on G(T) = T - Phi(T) in the two surface temperatures, from the
  // convective profile (no radiation).
  using D2 = DualN<2>;
  std::array<D2, WALL_RAD_N_PARAMS> p2;
  for (int i = 0; i < WALL_RAD_N_PARAMS; ++i) {
    p2[i] = D2::constant(pv[i]);
  }
  double Ta = 0.0;
  double Tb = 0.0;
  {
    const double R_tot = 1.0 / h_a + R_wall + 1.0 / h_b;
    const double q0 = (T_aw_a - T_aw_b) / R_tot;
    Ta = T_aw_a - q0 / h_a;
    Tb = T_aw_b + q0 / h_b;
  }
  int it = 0;
  bool converged = (kind_a == 0 && kind_b == 0);
  for (; !converged && it < kWallRadMaxIter; ++it) {
    const auto e = wall_rad_eval<2>(D2::seed(Ta, 0), D2::seed(Tb, 1), p2,
                                    kind_a, kind_b);
    const double g0 = Ta - e.Ts_a_new.v;
    const double g1 = Tb - e.Ts_b_new.v;
    const double j00 = 1.0 - e.Ts_a_new.d[0];
    const double j01 = -e.Ts_a_new.d[1];
    const double j10 = -e.Ts_b_new.d[0];
    const double j11 = 1.0 - e.Ts_b_new.d[1];
    const double det = j00 * j11 - j01 * j10;
    if (!(std::abs(det) > 0.0) || !std::isfinite(det)) {
      throw std::runtime_error(
          "wall_radiation_sides: singular surface-temperature Jacobian");
    }
    double dTa = -(j11 * g0 - j01 * g1) / det;
    double dTb = -(-j10 * g0 + j00 * g1) / det;
    // Keep the trial positive: at most halve either temperature.
    double lam = 1.0;
    while ((Ta + lam * dTa < 0.5 * Ta || Tb + lam * dTb < 0.5 * Tb) &&
           lam > 1e-6) {
      lam *= 0.5;
    }
    dTa *= lam;
    dTb *= lam;
    Ta += dTa;
    Tb += dTb;
    if (std::abs(dTa) <= kWallRadRelTol * std::max(1.0, std::abs(Ta)) &&
        std::abs(dTb) <= kWallRadRelTol * std::max(1.0, std::abs(Tb))) {
      converged = true;
    }
  }
  if (!converged || !(Ta > 0.0) || !(Tb > 0.0)) {
    throw std::runtime_error(
        "wall_radiation_sides: the surface temperatures did not converge");
  }

  // Total derivatives: seed the two surface temperatures (0, 1) and the
  // parameters (2..10); dT/dp = (I - dPhi/dT)^-1 dPhi/dp.
  constexpr int NP = WALL_RAD_N_PARAMS;
  constexpr int ND = NP + 2;
  using DN = DualN<ND>;
  std::array<DN, NP> pn;
  for (int i = 0; i < NP; ++i) {
    pn[i] = DN::seed(pv[i], i + 2);
  }
  const auto e = wall_rad_eval<ND>(DN::seed(Ta, 0), DN::seed(Tb, 1), pn, kind_a,
                                   kind_b);
  const double m00 = 1.0 - e.Ts_a_new.d[0];
  const double m01 = -e.Ts_a_new.d[1];
  const double m10 = -e.Ts_b_new.d[0];
  const double m11 = 1.0 - e.Ts_b_new.d[1];
  const double det = m00 * m11 - m01 * m10;
  if (!(std::abs(det) > 0.0) || !std::isfinite(det)) {
    throw std::runtime_error(
        "wall_radiation_sides: the surface temperatures are not isolated");
  }
  std::array<double, NP> dTa{};
  std::array<double, NP> dTb{};
  for (int j = 0; j < NP; ++j) {
    const double ba = e.Ts_a_new.d[j + 2];
    const double bb = e.Ts_b_new.d[j + 2];
    dTa[j] = (m11 * ba - m01 * bb) / det;
    dTb[j] = (-m10 * ba + m00 * bb) / det;
  }
  auto total = [&](const DN &y) {
    std::array<double, NP> out{};
    for (int j = 0; j < NP; ++j) {
      out[j] = y.d[j + 2] + y.d[0] * dTa[j] + y.d[1] * dTb[j];
    }
    return out;
  };

  WallRadiationSides r;
  r.H_a = e.H_a.v;
  r.T_a = e.T_a.v;
  r.H_b = e.H_b.v;
  r.T_b = e.T_b.v;
  r.T_s_a = Ta;
  r.T_s_b = Tb;
  r.h_r_a = e.hr_a.v;
  r.h_r_b = e.hr_b.v;
  r.q = e.q.v;
  r.iterations = it;
  r.d[0] = total(e.H_a);
  r.d[1] = total(e.T_a);
  r.d[2] = total(e.H_b);
  r.d[3] = total(e.T_b);
  const DN Ts_a = DN::seed(Ta, 0);
  const DN Ts_b = DN::seed(Tb, 1);
  const DN qa = e.hr_a * (pn[5] - Ts_a);
  const DN qb = e.hr_b * (Ts_b - pn[7]);
  r.q_rad_a = qa.v;
  r.q_rad_b = qb.v;
  r.dq_rad[0] = total(qa);
  r.dq_rad[1] = total(qb);
  return r;
}

WallAdiabaticFace wall_adiabatic_face(double T_aw, double h, int kind,
                                      double T_src, double k) {
  for (double v : {T_aw, h, T_src, k}) {
    require_finite(v, "wall_adiabatic_face input");
  }
  if (!(h > 0.0) || k < 0.0 || !(T_aw > 0.0)) {
    throw std::invalid_argument(
        "wall_adiabatic_face: needs h > 0, k >= 0 and T_aw > 0");
  }
  if (kind < 0 || kind > 2) {
    throw std::invalid_argument("wall_adiabatic_face: kind must be 0, 1 or 2");
  }
  if ((kind == 1 && !(T_src > 0.0)) || (kind == 2 && T_src < 0.0)) {
    throw std::invalid_argument(
        "wall_adiabatic_face: a gas temperature must be > 0, a surface's >= 0");
  }
  WallAdiabaticFace r;
  if (kind == 0 || k == 0.0) {
    r.T_s = T_aw;
    r.dT_s = {1.0, 0.0, 0.0, 0.0};
    return r;
  }

  // f(T) = h (T_aw - T) + h_r(T) (T_src - T), seeded in T (0) and the
  // parameters [T_aw, h, T_src, k] (1..4).
  constexpr int ND = WALL_ADIABATIC_N_PARAMS + 1;
  using DN = DualN<ND>;
  const DN p_Taw = DN::seed(T_aw, 1);
  const DN p_h = DN::seed(h, 2);
  const DN p_src = DN::seed(T_src, 3);
  const DN p_k = DN::seed(k, 4);
  auto balance = [&](double T) {
    const DN Ts = DN::seed(T, 0);
    const DN hr = radiative_htc<ND>(kind, Ts, p_src, p_k);
    return p_h * (p_Taw - Ts) + hr * (p_src - Ts);
  };

  // f falls strictly with T: the root lies between T_aw and T_src.
  double lo = std::max(std::min(T_aw, T_src), kWallRadTFloor);
  double hi = std::max(T_aw, T_src);
  double T = T_aw;
  int it = 0;
  bool converged = false;
  for (; it < kWallRadMaxIter; ++it) {
    const DN f = balance(T);
    if (f.v > 0.0) {
      lo = std::max(lo, T);
    } else if (f.v < 0.0) {
      hi = std::min(hi, T);
    } else {
      converged = true;
      break;
    }
    const double fT = f.d[0];
    double T_new = (fT < 0.0) ? T - f.v / fT : 0.5 * (lo + hi);
    if (!(T_new > lo && T_new < hi)) {
      T_new = 0.5 * (lo + hi); // Newton left the bracket: bisect
    }
    const double step = T_new - T;
    T = T_new;
    if (std::abs(step) <= kWallRadRelTol * std::max(1.0, std::abs(T)) ||
        hi - lo <= kWallRadRelTol * std::max(1.0, std::abs(T))) {
      converged = true;
      ++it;
      break;
    }
  }
  if (!converged || !(T > 0.0)) {
    throw std::runtime_error(
        "wall_adiabatic_face: the surface temperature did not converge");
  }

  const DN f = balance(T);
  const double fT = f.d[0];
  if (!(fT < 0.0) || !std::isfinite(fT)) {
    throw std::runtime_error(
        "wall_adiabatic_face: the surface temperature is not isolated");
  }
  r.T_s = T;
  r.iterations = it;
  for (int j = 0; j < WALL_ADIABATIC_N_PARAMS; ++j) {
    r.dT_s[static_cast<std::size_t>(j)] = -f.d[j + 1] / fT;
  }
  const DN Ts = DN::constant(T);
  const DN hr = radiative_htc<ND>(kind, Ts, p_src, p_k);
  r.h_r = hr.v;
  r.q_rad = (hr * (p_src - Ts)).v;
  return r;
}

} // namespace combaero
