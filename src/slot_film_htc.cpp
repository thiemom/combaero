#include "slot_film_htc.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "composition.h"
#include "thermo.h"
#include "transport.h"

namespace combaero {

namespace {

// ln(10): the blend half-width is stated in decades of x/s.
constexpr double kLn10 = 2.302585092994045684;
// g/mol per kg/mol: cp_mass = cp_molar / (W / 1000).
constexpr double kGramsPerKg = 1000.0;

// 8-point Gauss-Legendre on [-1, 1], on kBlendPanels panels in ln(x/s): the
// blend's integrand exp(g(u) + u) is smooth, so this is exact to round-off.
constexpr std::array<double, 8> kGLNode = {
    -0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
    -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
    0.7966664774136267,  0.9602898564975363};
constexpr std::array<double, 8> kGLWeight = {
    0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
    0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
    0.2223810344533745, 0.1012285362903763};
constexpr int kBlendPanels = 4;

void require(bool ok, const char *fn, const std::string &what) {
  if (!ok) {
    throw std::invalid_argument(std::string(fn) + ": " + what);
  }
}

void check_span(const char *fn, double a, double b) {
  require(std::isfinite(a) && std::isfinite(b), fn, "x/s must be finite");
  require(a >= 0.0, fn, "x_start_over_s must be >= 0");
  require(b >= a, fn, "x_end_over_s must be >= x_start_over_s");
  require(b > 0.0, fn,
          "a point value at the slot lip (x/s = 0) is infinite; give a span "
          "or x/s > 0");
}

// Seban's f(x/s) = St_s Re_s^0.3: curve B, the C1 blend, curve C.
struct SebanCurve {
  double lnB, lnC, uL, uR;

  SebanCurve(double C_near, double C_far) {
    lnB = std::log(C_near);
    lnC = std::log(C_far);
    // Crossing: ln C_B - 0.3 u = ln C_C - 0.6 u.
    const double u_star =
        (lnC - lnB) / (seban1960::x_exponent_far - seban1960::x_exponent_near);
    const double d = seban1960::blend_half_width_decades * kLn10;
    uL = u_star - d;
    uR = u_star + d;
  }

  double g(double u) const {
    const double gB = lnB - seban1960::x_exponent_near * u;
    const double gC = lnC - seban1960::x_exponent_far * u;
    if (u <= uL) {
      return gB;
    }
    if (u >= uR) {
      return gC;
    }
    const double t = (u - uL) / (uR - uL);
    return gB + (gC - gB) * t * t * (3.0 - 2.0 * t);
  }

  double f(double X) const { return std::exp(g(std::log(X))); }

  // Integral of f over [a, b], 0 <= a <= b.
  double integral(double a, double b) const {
    if (!(b > a)) {
      return 0.0;
    }
    const double XL = std::exp(uL), XR = std::exp(uR);
    double sum = 0.0;
    const double pB = 1.0 - seban1960::x_exponent_near;
    const double pC = 1.0 - seban1960::x_exponent_far;
    const double hiB = std::min(b, XL);
    if (a < hiB) {
      sum += std::exp(lnB) * (std::pow(hiB, pB) - std::pow(a, pB)) / pB;
    }
    const double loC = std::max(a, XR);
    if (loC < b) {
      sum += std::exp(lnC) * (std::pow(b, pC) - std::pow(loC, pC)) / pC;
    }
    const double lo = std::max(a, XL), hi = std::min(b, XR);
    if (lo < hi) {
      const double u0 = std::log(lo), u1 = std::log(hi);
      const double hu = (u1 - u0) / kBlendPanels;
      for (int p = 0; p < kBlendPanels; ++p) {
        const double mid = u0 + (p + 0.5) * hu;
        for (std::size_t q = 0; q < kGLNode.size(); ++q) {
          const double u = mid + 0.5 * hu * kGLNode[q];
          sum += 0.5 * hu * kGLWeight[q] * std::exp(g(u) + u);
        }
      }
    }
    return sum;
  }

  int branch(double a, double b) const {
    const double XL = std::exp(uL), XR = std::exp(uR);
    if (b <= XL) {
      return SLOT_FILM_BRANCH_NEAR;
    }
    if (a >= XR) {
      return SLOT_FILM_BRANCH_FAR;
    }
    return SLOT_FILM_BRANCH_BLEND;
  }
};

// The in-box part of the span: its upper end and the share it is of the span.
struct InBox {
  double hi = 0.0;
  double fraction = 0.0;
};

InBox in_box_part(double a, double b, double X_max) {
  InBox r;
  if (a > X_max || (a == X_max && b > a)) {
    return r;
  }
  r.hi = std::min(b, X_max);
  r.fraction = (b > a) ? (r.hi - a) / (b - a) : 1.0;
  return r;
}

// w(M): 0 below M_handback, 1 above M_handback + band, C1 smoothstep between.
void handback_weight(double M, double M_handback, double &w, double &dw) {
  const double t = (M - M_handback) / SLOT_FILM_HTC_HANDBACK_BAND;
  if (t <= 0.0) {
    w = 0.0;
    dw = 0.0;
  } else if (t >= 1.0) {
    w = 1.0;
    dw = 0.0;
  } else {
    w = t * t * (3.0 - 2.0 * t);
    dw = 6.0 * t * (1.0 - t) / SLOT_FILM_HTC_HANDBACK_BAND;
  }
}

struct SlotProps {
  double mu, dmu_dT, k, dk_dT, cp, dcp_dT;
};

SlotProps slot_props(double T, double P, const std::vector<double> &X) {
  const TransportDerivatives tr = transport_and_dT(T, P, X);
  const double W = mwmix(X) / kGramsPerKg;
  return {tr.mu, tr.dmu_dT, tr.k, tr.dk_dT, cp_mass(T, X), dcp_dT(T, X) / W};
}

void check_dimensional(const char *fn, double x_start, double x_end, double s,
                       double G_s, double M, double h_duct, double T_s) {
  require(std::isfinite(s) && s > 0.0, fn, "s (slot height) must be > 0");
  require(std::isfinite(G_s) && G_s > 0.0, fn, "G_s must be > 0");
  require(std::isfinite(M) && M >= 0.0, fn, "M must be >= 0");
  require(std::isfinite(h_duct), fn, "h_duct must be finite");
  require(std::isfinite(T_s) && T_s > 0.0, fn, "T_s must be > 0");
  check_span(fn, x_start / s, x_end / s);
}

// h = h_duct + w f_x (h_film - h_duct) and its partials, from the film's own.
void apply_handback(SlotFilmHtc &r, double h_film, double dhf_dG,
                    double dhf_dT, double fraction, double M,
                    double M_handback, double h_duct) {
  double w = 0.0, dw = 0.0;
  handback_weight(M, M_handback, w, dw);
  const double a = w * fraction;
  r.h_film = h_film;
  r.fraction_in_box = fraction;
  r.weight_M = w;
  r.h = h_duct + a * (h_film - h_duct);
  r.dh_dG_s = a * dhf_dG;
  r.dh_dT_s = a * dhf_dT;
  r.dh_dM = dw * fraction * (h_film - h_duct);
  r.dh_dh_duct = 1.0 - a;
}

} // namespace

SlotFilmStanton slot_film_stanton_seban_1960(double x_start_over_s,
                                             double x_end_over_s, double Re_s,
                                             double C_near, double C_far) {
  const char *fn = "slot_film_stanton_seban_1960";
  check_span(fn, x_start_over_s, x_end_over_s);
  require(std::isfinite(Re_s) && Re_s > 0.0, fn, "Re_s must be > 0");
  require(std::isfinite(C_near) && C_near > 0.0 && std::isfinite(C_far) &&
              C_far > 0.0,
          fn, "C_near and C_far must be > 0");
  const double a = x_start_over_s, b = x_end_over_s;
  const SebanCurve curve(C_near, C_far);
  SlotFilmStanton r;
  r.x_in_box = a >= seban1960::x_over_s_min && b <= seban1960::x_over_s_max;
  r.Re_in_box = Re_s >= seban1960::Re_s_min && Re_s <= seban1960::Re_s_max;
  const InBox box = in_box_part(a, b, seban1960::x_over_s_max);
  r.fraction_in_box = box.fraction;
  if (box.fraction <= 0.0) {
    return r;
  }
  const double f_mean =
      (box.hi > a) ? curve.integral(a, box.hi) / (box.hi - a) : curve.f(a);
  const double re_factor = std::pow(Re_s, -seban1960::Re_exponent);
  r.St = re_factor * f_mean;
  r.dSt_dRe = -seban1960::Re_exponent * r.St / Re_s;
  r.branch = curve.branch(a, box.hi);
  return r;
}

SlotFilmNusselt slot_film_nusselt_ballal_lefebvre_1973(double x_start_over_s,
                                                      double x_end_over_s,
                                                      double Re_s) {
  const char *fn = "slot_film_nusselt_ballal_lefebvre_1973";
  check_span(fn, x_start_over_s, x_end_over_s);
  require(std::isfinite(Re_s) && Re_s > 0.0, fn, "Re_s must be > 0");
  const double a = x_start_over_s, b = x_end_over_s;
  SlotFilmNusselt r;
  r.x_in_box = b <= ballal_lefebvre1973::x_over_s_max;
  const InBox box = in_box_part(a, b, ballal_lefebvre1973::x_over_s_max);
  r.fraction_in_box = box.fraction;
  if (box.fraction <= 0.0) {
    return r;
  }
  // h s / k = (Nu_x / (x/s)) = C Re^0.8 (x/s)^(0.44 - 1).
  const double p = ballal_lefebvre1973::x_exponent;
  const double x_mean =
      (box.hi > a)
          ? (std::pow(box.hi, p) - std::pow(a, p)) / (p * (box.hi - a))
          : std::pow(a, p - 1.0);
  r.Nu_s = ballal_lefebvre1973::C_wall_jet *
           std::pow(Re_s, ballal_lefebvre1973::Re_exponent) * x_mean;
  r.dNu_dRe = ballal_lefebvre1973::Re_exponent * r.Nu_s / Re_s;
  r.branch = SLOT_FILM_BRANCH_WALL_JET;
  return r;
}

SlotFilmHtc slot_film_htc_seban_1960(double x_start, double x_end, double s,
                                     double G_s, double M, double h_duct,
                                     double T_s, double P_s,
                                     const std::vector<double> &X_s) {
  check_dimensional("slot_film_htc_seban_1960", x_start, x_end, s, G_s, M,
                    h_duct, T_s);
  const SlotProps pr = slot_props(T_s, P_s, X_s);
  SlotFilmHtc r;
  r.Re_s = G_s * s / pr.mu;
  const SlotFilmStanton st =
      slot_film_stanton_seban_1960(x_start / s, x_end / s, r.Re_s);
  r.branch = st.branch;
  r.x_in_box = st.x_in_box;
  r.Re_in_box = st.Re_in_box;
  r.M_in_box = M >= seban1960::M_min && M <= seban1960::M_max;
  // h_film = St(Re) G cp, Re = G s / mu(T).
  const double h_film = st.St * G_s * pr.cp;
  const double dRe_dG = s / pr.mu;
  const double dRe_dT = -r.Re_s * pr.dmu_dT / pr.mu;
  const double dhf_dG = pr.cp * (st.St + G_s * st.dSt_dRe * dRe_dG);
  const double dhf_dT = G_s * (pr.cp * st.dSt_dRe * dRe_dT + st.St * pr.dcp_dT);
  apply_handback(r, h_film, dhf_dG, dhf_dT, st.fraction_in_box, M,
                 seban1960::M_handback, h_duct);
  return r;
}

SlotFilmHtc slot_film_htc_ballal_lefebvre_1973(double x_start, double x_end,
                                               double s, double G_s, double M,
                                               double h_duct, double T_s,
                                               double P_s,
                                               const std::vector<double> &X_s) {
  check_dimensional("slot_film_htc_ballal_lefebvre_1973", x_start, x_end, s,
                    G_s, M, h_duct, T_s);
  const SlotProps pr = slot_props(T_s, P_s, X_s);
  SlotFilmHtc r;
  r.Re_s = G_s * s / pr.mu;
  const SlotFilmNusselt nu =
      slot_film_nusselt_ballal_lefebvre_1973(x_start / s, x_end / s, r.Re_s);
  r.branch = nu.branch;
  r.x_in_box = nu.x_in_box;
  // No Reynolds range is stated for Eq. 6; m > 1.3, no upper limit.
  r.Re_in_box = true;
  r.M_in_box = M >= ballal_lefebvre1973::M_handback;
  // h_film = Nu_s(Re) k / s.
  const double h_film = nu.Nu_s * pr.k / s;
  const double dRe_dG = s / pr.mu;
  const double dRe_dT = -r.Re_s * pr.dmu_dT / pr.mu;
  const double dhf_dG = pr.k / s * nu.dNu_dRe * dRe_dG;
  const double dhf_dT = (nu.Nu_s * pr.dk_dT + pr.k * nu.dNu_dRe * dRe_dT) / s;
  apply_handback(r, h_film, dhf_dG, dhf_dT, nu.fraction_in_box, M,
                 ballal_lefebvre1973::M_handback, h_duct);
  return r;
}

} // namespace combaero
