#pragma once

// -----------------------------------------------------------------------------
// Combustor liner radiation (#793): the Lefebvre-Herbert lumped, gray method
// as printed in Lefebvre & Ballal, Gas Turbine Combustion, 3rd ed. (2010),
// Ch. 8, Eqs. 8.7-8.17.
//
//   R1 = 0.5 sigma (1 + eps_w) eps_g Tg^1.5 (Tg^2.5 - Tw^2.5)          (8.7)
//   eps_g = 1 - exp[-290 P L (q l_b)^0.5 Tg^-1.5],  P in kPa          (8.11)
//   L = 336 / H^2,  H the fuel hydrogen content, mass %               (8.15)
//   l_b = 3.4 V / A                                                   (8.10)
//   R2 = sigma eps_w eps_c / (eps_c + eps_w (1 - eps_c) A_w/A_c)
//        (Tw^4 - Tc^4)                                                (8.17)
//
// What it is: one gas temperature, one beam length, one gray emissivity per
// surface. Eq. 8.7 replaces the gas absorptivity by alpha_g/eps_g =
// (Tg/Tw)^1.5 (Eq. 8.6) and the wall's absorptivity by 0.5 (1 + eps_w)
// (Eq. 8.5). Eq. 8.9 (L = 1) is Reeves' (NGTE Memo M285, 1956) non-luminous
// fit for kerosine products, written in the fuel/air ratio q, so it is a
// kerosine-product law: it is not a composition-based emissivity and does not
// carry over to natural gas or hydrogen products. The book prints Eq. 8.9 as
// "2 - exp"; Eq. 8.11 and the worked example use "1 - exp" (erratum).
//
// The luminosity factor L stands in for soot. Shipped: Eq. 8.15 (Lefebvre's
// fit to engine data). Not shipped, named so a user can supply L directly:
//   Eq. 8.12  L = 7.53 (C/H - 5.5)^0.84   (Lefebvre & Herbert 1960, original)
//   Eq. 8.13  L = 3 (C/H - 5.2)^0.75      (later modification)
//   Eq. 8.14  L = 0.0691 (C/H - 1.82)^2.71 (Kretschmer & Odgers 1978)
// with C/H the fuel's carbon/hydrogen mass ratio. The four forms disagree.
//
// Mean beam length: Eq. 8.10's 3.4 V/A (after Fishenden & Saunders 1950).
// The textbook engineering value is 0.9 x 4 V/A = 3.6 V/A (Howell, Menguc,
// Daun & Siegel, Thermal Radiation Heat Transfer 7th ed., on-line App. B);
// the two differ by 6 %.
// -----------------------------------------------------------------------------

#include <array>
#include <tuple>
#include <vector>

namespace combaero {

// Stefan-Boltzmann constant [W/(m^2 K^4)], CODATA 2018 (exact in SI). The
// book's worked examples use 5.67e-8 (0.007 % lower).
constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;

// Reeves' coefficient in Eqs. 8.9/8.11 [kPa^-1 m^-0.5 K^1.5]: with P in kPa,
// l_b in m and Tg in K the exponent is dimensionless.
constexpr double REEVES_EMISSIVITY_COEFF = 290.0;

// Pa per kPa: the kernels take P in Pa and form Eq. 8.11's P in kPa.
constexpr double LINER_RADIATION_PA_PER_KPA = 1000.0;

// Eq. 8.15: L = 336 / H^2, H in mass percent.
constexpr double LUMINOSITY_HYDROGEN_COEFF = 336.0;

// Mean beam length coefficients, l_b = c V / A: Lefebvre Eq. 8.10 (3.4,
// shipped default) and the textbook 0.9 x 4 = 3.6 (Howell App. B).
constexpr double BEAM_LENGTH_COEFF_LEFEBVRE = 3.4;
constexpr double BEAM_LENGTH_COEFF_HOWELL = 3.6;

// Atomic masses [g/mol] (IUPAC conventional values) for the elemental
// fuel/air ratio; they reproduce the species data's molar masses (CH4
// 16.043).
constexpr double ATOMIC_MASS_C = 12.011;
constexpr double ATOMIC_MASS_H = 1.008;

// Eq. 8.11 gas emissivity and its partial derivatives.
struct GasEmissivityResult {
  double eps = 0.0;       // eps_g [-]
  double deps_dT = 0.0;   // d eps_g / d Tg [1/K]
  double deps_dP = 0.0;   // d eps_g / d P [1/Pa]
  double deps_dq = 0.0;   // d eps_g / d q [-]
  double deps_dl_b = 0.0; // d eps_g / d l_b [1/m]
  double deps_dL = 0.0;   // d eps_g / d L [-]
};

// eps_g = 1 - exp[-290 (P/1000) L (q l_b)^0.5 Tg^-1.5] (Eq. 8.11; L = 1 is
// Reeves' non-luminous Eq. 8.9).
//   T_g : mean gas temperature [K], > 0
//   P   : gas pressure [Pa], >= 0
//   q   : fuel/air mass ratio [-], >= 0
//   l_b : mean beam length [m], >= 0
//   L   : luminosity factor [-], >= 0
// q l_b = 0 gives eps_g = 0 with zero q and l_b derivatives (the one-sided
// derivative there is infinite). Throws on T_g <= 0 or a negative input.
GasEmissivityResult gas_emissivity_luminous(double T_g, double P, double q,
                                            double l_b, double L = 1.0);

// Eq. 8.15: L = 336 / H^2 and dL/dH, H the fuel hydrogen mass percent
// (e.g. 13.9 for Jet A). Throws on H <= 0.
std::tuple<double, double> luminosity_factor_hydrogen(double H_mass_percent);

// Eq. 8.7, gas to liner radiation flux, and its partial derivatives at fixed
// eps_g (chain eps_g's own derivatives for the total).
struct LinerGasRadiationResult {
  double R1 = 0.0;         // [W/m^2], positive gas -> wall
  double dR1_dT_g = 0.0;   // [W/(m^2 K)] at fixed eps_g
  double dR1_dT_w = 0.0;   // [W/(m^2 K)]
  double dR1_deps_g = 0.0; // [W/m^2]
  double dR1_deps_w = 0.0; // [W/m^2]
};

// R1 = 0.5 sigma (1 + eps_w) eps_g Tg^1.5 (Tg^2.5 - Tw^2.5). Throws on
// T_g <= 0 or T_w < 0.
LinerGasRadiationResult liner_gas_radiation(double T_g, double T_w,
                                            double eps_g, double eps_w);

// Eq. 8.17, liner to casing radiation flux [W/m^2] per unit LINER area, gray
// two-surface enclosure with view factor 1:
//   R2 = sigma eps_w eps_c / (eps_c + eps_w (1 - eps_c) A_w/A_c) (Tw^4 - Tc^4)
// Returns (R2, dR2/dT_w, dR2/dT_c). Throws on an emissivity outside (0, 1],
// area_ratio <= 0 or a negative temperature.
std::tuple<double, double, double>
liner_casing_radiation(double T_w, double T_c, double eps_w, double eps_c,
                       double area_ratio);

// sigma eps_w eps_c / (eps_c + eps_w (1 - eps_c) A_w/A_c) [W/(m^2 K^4)], the
// Eq. 8.17 exchange factor (0.4 sigma and 0.6 sigma are the book's Eqs.
// 8.18/8.19 shortcuts for aluminium and steel casings).
double liner_casing_exchange_factor(double eps_w, double eps_c,
                                    double area_ratio);

// Eq. 8.10: l_b = coeff V / A [m]. Throws on V < 0, A <= 0 or coeff <= 0.
double mean_beam_length(double volume, double area,
                        double coeff = BEAM_LENGTH_COEFF_LEFEBVRE);

// Fuel/air mass ratio of a mixture from its elements: the mass of its carbon
// and hydrogen atoms over the mass of everything else, q = w / (1 - w),
// w = sum_k Y_k (C_k M_C + H_k M_H) / M_k. Elemental, so the same q for the
// reactants and for their products. Hydrogen in water vapour the air brought
// counts as fuel hydrogen; a fuel's own oxygen or nitrogen (an alcohol, NH3)
// counts as air. Returns (q, dq/dY) in the mass fractions Y. Throws when the
// mixture is all C and H (w >= 1) or Y has the wrong length.
std::tuple<double, std::vector<double>>
elemental_fuel_air_ratio(const std::vector<double> &Y);

// -----------------------------------------------------------------------------
// One wall with radiation on either face (#793): the local balance solved for
// the two surface temperatures, folded into an effective h and source
// temperature per side so the convective wall kernels apply unchanged.
//
// Per side s the heat into the wall face is
//   q_s = h_s (T_aw_s - T_s) + h_r,s (T_src,s - T_s)
// (side B's sign mirrored), with the radiation written exactly as a
// conductance on its own driving temperature:
//   kind 1 (gas, Eq. 8.7):    h_r = k Tsrc^1.5 (Tsrc^2.5 - T_s^2.5)/(Tsrc - T_s)
//                             k = 0.5 sigma (1 + eps_w) eps_g, Tsrc = Tg
//   kind 2 (surface, Eq. 8.17): h_r = k (T_s + Tsrc)(T_s^2 + Tsrc^2)
//                             k = the exchange factor, Tsrc = Tc
// both divided differences, finite at T_s = Tsrc. Effective side:
//   H_s = h_s + h_r,s,  T*_s = (h_s T_aw_s + h_r,s Tsrc,s) / H_s,
// so H_s (T*_s - T_s) is the face's total flux and the two-face profile
// q = (T*_a - T*_b) / (1/H_a + R_wall + 1/H_b) holds. The surface
// temperatures are its fixed point T_a = T*_a - q/H_a, T_b = T*_b + q/H_b,
// solved by Newton; the Jacobian by implicit differentiation (DualN).
// -----------------------------------------------------------------------------

enum class WallRadiationKind : int { None = 0, Gas = 1, Surface = 2 };

// Parameter order of WallRadiationSides::d columns.
constexpr int WALL_RAD_N_PARAMS = 9;
// [T_aw_a, h_a, T_aw_b, h_b, R_wall, T_src_a, k_a, T_src_b, k_b]

struct WallRadiationSides {
  double H_a = 0.0, T_a = 0.0, H_b = 0.0, T_b = 0.0; // effective sides
  double T_s_a = 0.0, T_s_b = 0.0;                   // surface temperatures [K]
  double h_r_a = 0.0, h_r_b = 0.0;                   // radiative conductances
  double q_rad_a = 0.0; // radiative flux into face A, h_r,a (T_src,a - T_s,a)
  double q_rad_b = 0.0; // radiative flux out of face B, h_r,b (T_s,b - T_src,b)
  double q = 0.0;       // total flux through the wall, A -> B [W/m^2]
  int iterations = 0;
  // d(H_a, T_a, H_b, T_b)/d(params), total (through the surface fixed point)
  std::array<std::array<double, WALL_RAD_N_PARAMS>, 4> d{};
  // d(q_rad_a, q_rad_b)/d(params), total
  std::array<std::array<double, WALL_RAD_N_PARAMS>, 2> dq_rad{};
};

// h_a, h_b > 0 (the caller floors them); R_wall >= 0 the conductive plus
// fouling resistance between the faces [m^2 K/W]; k_s >= 0. kind None
// ignores that side's T_src and k. Throws if the surface temperatures do
// not converge or leave T > 0.
WallRadiationSides wall_radiation_sides(double T_aw_a, double h_a,
                                        double T_aw_b, double h_b,
                                        double R_wall, int kind_a,
                                        double T_src_a, double k_a, int kind_b,
                                        double T_src_b, double k_b);

// -----------------------------------------------------------------------------
// The face of a wall whose OTHER side is adiabatic (#970): no heat crosses the
// wall, so the face's own balance closes on itself,
//   h (T_aw - T_s) + h_r(T_s) (T_src - T_s) = 0,
// with h_r as in wall_radiation_sides (kind 0 none, 1 gas Eq. 8.7, 2 surface
// Eq. 8.17). Without radiation T_s = T_aw. The left side falls strictly with
// T_s, so the root is unique and lies between T_aw and T_src; it is solved by
// a bracketed Newton and differentiated implicitly.
// -----------------------------------------------------------------------------

// Parameter order of WallAdiabaticFace::dT_s.
constexpr int WALL_ADIABATIC_N_PARAMS = 4;
// [T_aw, h, T_src, k]

struct WallAdiabaticFace {
  double T_s = 0.0;   // surface temperature [K]
  double h_r = 0.0;   // radiative conductance at T_s [W/(m^2*K)]
  double q_rad = 0.0; // radiative flux into the face, h_r (T_src - T_s) [W/m^2]
  int iterations = 0;
  // dT_s/d[T_aw, h, T_src, k]
  std::array<double, WALL_ADIABATIC_N_PARAMS> dT_s{};
};

// h > 0 (the caller floors it); k >= 0; kind None ignores T_src and k. Throws
// on a non-finite input, an invalid kind or temperature, or no convergence.
WallAdiabaticFace wall_adiabatic_face(double T_aw, double h, int kind,
                                      double T_src, double k);

} // namespace combaero
