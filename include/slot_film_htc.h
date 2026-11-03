#pragma once

// -----------------------------------------------------------------------------
// Film-side heat transfer coefficient downstream of a 2-D slot (#962).
//
// Near a film slot the hot-side boundary layer is a fresh layer (or a wall
// jet) of the slot fluid, and its h is set by the slot flow and the distance
// from the lip, not by the duct's developed flow. These kernels give that h
// on the ADIABATIC-WALL temperature of the same flow (Seban's Eq. 3), so the
// wall heat flux is q = h (T_aw - T_w) with T_aw from the film effectiveness.
//
// SHIPPED DEFAULT: Seban, R. A. (1960). "Heat transfer and effectiveness for
// a turbulent boundary layer with tangential fluid injection." ASME J. Heat
// Transfer 82(4), 303-313. Stanton form on SLOT-fluid properties, Fig. 12:
//
//   St_s = h / (rho_s u_s cp_s) = C_B Re_s^-0.3 (x/s)^-0.3     curve B (Eq. 13)
//   St_s                         = C_C Re_s^-0.3 (x/s)^-0.6     curve C (Eq. 14)
//   Re_s = rho_s u_s s / mu_s
//
// with C_B = 0.135 as DRAWN (the paper prints 0.10) and C_C = 0.40 as drawn
// (printed 0.41). The printed 0.10 contradicts Seban's own curve B, his
// Fig. 9 and his data: the B/C crossing of the drawn curves sits at x/s 37
// (the printed pair would cross at x/s 110, not at the "x/s = 40" the text
// gives), Fig. 9 back-solves to 0.140-0.144, and the digitised Fig. 12
// points fit 0.134 (validation/cooling/extractions/seban_1960_slot_film_h.md;
// signed off by the user 2026-11-03). The two curves are joined by a C1 blend
// in ln(x/s) over +-0.1 decade around their crossing.
//
// What Lefebvre & Ballal's Eq. 8.41, Nu_x = 0.069 [Re_s x/s]^0.7, is: Eq. 13
// with the printed 0.10 times Pr 0.69; it inherits the misprint. The
// Cranfield SME4 report's 0.090 is Seban's curve and data (0.13-0.135) x Pr 0.69.
//
// NAMED OPTION, not the default: Ballal & Lefebvre (1973), Cranfield SME
// Report 4, Eq. 6 = Lefebvre & Ballal (2010) Eq. 8.43/8.44, the m > 1.3
// wall-jet form
//
//   Nu_x = h x / k_s = 0.10 Re_s^0.8 (x/s)^0.44
//
// shipped so that the book's film-cooled worked example can be reproduced.
// It sits 20-67 % ABOVE the Cranfield rig's own data at low Re (SME4 Fig. 2;
// only Pai & Whitelaw's points lie on it). Seban's Fig. 12 is itself the
// high-injection data (m 1.06-8.91), so Seban is the sourced choice for all
// m > 1.
//
// VALIDITY BOX AND HAND-BACK. Seban's heat-transfer runs: x/s about 4-300,
// mass-velocity ratio m 1.06-8.91, Re_s about 1.4e3-8e3, air into air at a
// small temperature ratio. Beyond the far end (x/s 300; the wall jet's
// x/s 50, where Ballal & Lefebvre switch to Colburn) the side hands back to
// the duct's own h, segment by segment. Below m = 1 Seban finds plain
// flat-plate h (his Fig. 8), so the side hands back to the duct h there too:
// a C1 smoothstep in m over [m_min, m_min + 0.1] (a numerical transition,
// not physics). Outside the box otherwise (x/s < 4, Re_s, m above the data)
// the form is applied and flagged.
//
// ELEMENT MEAN. A wall segment's h is the mean of h(x) over its span
// [x_start, x_end] (integrated, not a point value), consistent with a
// ThermalWall's element-mean h; x_start == x_end is the point value.
// -----------------------------------------------------------------------------

#include <vector>

namespace combaero {

namespace seban1960 {

// Fig. 12 curve B (Eq. 13): the constant as DRAWN, and as printed.
constexpr double C_near = 0.135;
constexpr double C_near_printed = 0.10;
// Fig. 12 curve C (Eq. 14): the constant as drawn (trace 0.400), and as
// printed.
constexpr double C_far = 0.40;
constexpr double C_far_printed = 0.41;
// Exponents: Re_s^-0.3 on both curves; (x/s)^-0.3 near, ^-0.6 far.
constexpr double Re_exponent = 0.3;
constexpr double x_exponent_near = 0.3;
constexpr double x_exponent_far = 0.6;
// Half-width of the C1 blend of curves B and C, decades of x/s, centred on
// their crossing (x/s = (C_far / C_near)^(1/0.3) = 37.4).
constexpr double blend_half_width_decades = 0.1;

// The heat-transfer runs (Table 2 and Fig. 12): x/s, mass-velocity ratio,
// slot Reynolds number. Not stated as limits by the paper.
constexpr double x_over_s_min = 4.0;
constexpr double x_over_s_max = 300.0;
constexpr double M_min = 1.06;
constexpr double M_max = 8.91;
constexpr double Re_s_min = 1.4e3;
constexpr double Re_s_max = 8.0e3;
// Below this mass-velocity ratio the side hands back to the duct's h:
// Seban's Fig. 8 finds plain flat-plate h for m < 1.
constexpr double M_handback = 1.0;

} // namespace seban1960

namespace ballal_lefebvre1973 {

// SME4 Eq. 6 (book Eq. 8.43): Nu_x = 0.10 Re_s^0.8 (x/s)^0.44, m > 1.3.
constexpr double C_wall_jet = 0.10;
constexpr double Re_exponent = 0.8;
constexpr double x_exponent = 0.44;
// Its stated range: m > 1.3, up to x/s 50 (SME4 uses Colburn beyond).
constexpr double M_handback = 1.3;
constexpr double x_over_s_max = 50.0;

} // namespace ballal_lefebvre1973

// Width of the smoothstep in m over which a side hands back to the duct h,
// [M_handback, M_handback + SLOT_FILM_HTC_HANDBACK_BAND].
constexpr double SLOT_FILM_HTC_HANDBACK_BAND = 0.1;

// Branch codes reported by the kernels.
constexpr int SLOT_FILM_BRANCH_NEAR = 0;      // Seban curve B only
constexpr int SLOT_FILM_BRANCH_BLEND = 1;     // the span reaches the B/C blend
constexpr int SLOT_FILM_BRANCH_FAR = 2;       // Seban curve C only
constexpr int SLOT_FILM_BRANCH_WALL_JET = 3;  // Ballal & Lefebvre Eq. 6
constexpr int SLOT_FILM_BRANCH_NONE = -1;     // the whole span is past x/s max

// Seban's Stanton number, span mean over the IN-BOX part of [x_start_over_s,
// x_end_over_s] (x/s <= x_over_s_max).
struct SlotFilmStanton {
  double St = 0.0;              // mean St_s over the in-box part [-]
  double dSt_dRe = 0.0;         // d St / d Re_s
  double fraction_in_box = 0.0; // share of the span at x/s <= x_over_s_max
  int branch = SLOT_FILM_BRANCH_NONE;
  bool x_in_box = false;  // the whole span in [x_over_s_min, x_over_s_max]
  bool Re_in_box = false; // Re_s in [Re_s_min, Re_s_max]
};

// St_s = Re_s^-0.3 f(x/s), f curve B, the blend, curve C (as drawn).
// Throws on x_start_over_s < 0, x_end_over_s < x_start_over_s, a point value
// at x/s = 0, or Re_s <= 0. C_near and C_far are the curve constants
// (default: as drawn; seban1960::C_near_printed / C_far_printed give the
// print, whose curves cross at x/s 110); the blend sits at their crossing.
SlotFilmStanton slot_film_stanton_seban_1960(
    double x_start_over_s, double x_end_over_s, double Re_s,
    double C_near = seban1960::C_near, double C_far = seban1960::C_far);

// Ballal & Lefebvre's wall-jet Nusselt number on the slot height, span mean
// over the in-box part: Nu_s = h s / k_s = 0.10 Re_s^0.8 (x/s)^-0.56.
struct SlotFilmNusselt {
  double Nu_s = 0.0;   // mean h s / k_s over the in-box part [-]
  double dNu_dRe = 0.0;
  double fraction_in_box = 0.0;
  int branch = SLOT_FILM_BRANCH_NONE;
  bool x_in_box = false;
};
SlotFilmNusselt slot_film_nusselt_ballal_lefebvre_1973(double x_start_over_s,
                                                      double x_end_over_s,
                                                      double Re_s);

// A wall side's h under the slot film, dimensional: the film h over its
// in-box span, the duct's h past it, and the duct's h below the hand-back m.
//
//   h = h_duct + w(M) f_x (h_film - h_duct)
//
// f_x the share of the span inside x/s max and w(M) the smoothstep from 0 at
// M_handback to 1 at M_handback + SLOT_FILM_HTC_HANDBACK_BAND. Properties
// (mu_s, cp_s, k_s) are the slot fluid's at T_s, P_s, X_s (mole fractions).
// G_s and M are taken as independent inputs (M = G_s / G_gas: chain both).
struct SlotFilmHtc {
  double h = 0.0;          // the side's h [W/(m^2 K)]
  double dh_dG_s = 0.0;    // at fixed M [W/(m^2 K) / (kg/(m^2 s))]
  double dh_dT_s = 0.0;    // through the slot fluid's properties [W/(m^2 K^2)]
  double dh_dM = 0.0;      // through the hand-back weight
  double dh_dh_duct = 0.0; // 1 - w f_x
  double h_film = 0.0;     // the film's span mean over its in-box part
  double Re_s = 0.0;
  double fraction_in_box = 0.0;
  double weight_M = 0.0;   // w(M)
  int branch = SLOT_FILM_BRANCH_NONE;
  bool x_in_box = false;
  bool Re_in_box = false;
  bool M_in_box = false;
};

// Seban 1960 (the default). x_start, x_end [m] downstream of the slot lip,
// s the slot height [m], G_s = m_s / A_slot [kg/(m^2 s)] > 0, M the
// mass-velocity ratio G_s / G_gas, h_duct the side's own h [W/(m^2 K)].
SlotFilmHtc slot_film_htc_seban_1960(double x_start, double x_end, double s,
                                     double G_s, double M, double h_duct,
                                     double T_s, double P_s,
                                     const std::vector<double> &X_s);

// Ballal & Lefebvre 1973 Eq. 6 (book Eq. 8.44), the named option.
SlotFilmHtc slot_film_htc_ballal_lefebvre_1973(double x_start, double x_end,
                                               double s, double G_s, double M,
                                               double h_duct, double T_s,
                                               double P_s,
                                               const std::vector<double> &X_s);

} // namespace combaero
