#ifndef COOLING_CORRELATIONS_H
#define COOLING_CORRELATIONS_H

#include "correlation_status.h"

#include <limits>
#include <tuple>
#include <utility>
#include <vector>

namespace combaero::cooling {


// Thermal-hydraulic performance factor (Webb & Eckert 1972)
// eta = (Nu/Nu0) / (f/f0)^(1/3)
// Compares heat transfer gain against pumping power penalty at equal velocity.
// eta > 1: surface is a net improvement; eta < 1: pressure-drop generator.
// Generic: works for ribs, dimples, pin fins, impingement, or any surface.
//
// Parameters:
//   Nu_ratio: Nu_enhanced / Nu_smooth [-]
//   f_ratio : f_enhanced  / f_smooth  [-]
//
// Returns: thermal-hydraulic performance factor eta [-]
//
// Source: Webb, R.L. & Eckert, E.R.G. (1972) Int. J. Heat Mass Transfer 15(8), 1647-1658
double thermal_performance_factor(double Nu_ratio, double f_ratio);


// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

// Adiabatic wall temperature from cooling effectiveness
// Uses definition: eta = (T_hot - T_aw) / (T_hot - T_coolant)
//
// Parameters:
//   T_hot      : mainstream hot-gas temperature [K]
//   T_coolant  : coolant supply temperature [K]
//   eta        : adiabatic cooling effectiveness [-]
//
// Returns: adiabatic wall temperature [K]
double adiabatic_wall_temperature(double T_hot, double T_coolant, double eta);

// Film/effusion-cooled wall heat flux with single wall layer
// Computes q based on adiabatic wall temperature and thermal resistance network:
//   q = U * (T_aw - T_coolant)
// where U uses both convective sides and wall conduction.
//
// Parameters:
//   T_hot      : mainstream hot-gas driving temperature [K]
//                For M < 0.3: pass static temperature T_static.
//                For M > 0.3 (combustor liner, turbine cooling): pass
//                T_adiabatic_wall() from stagnation.h as T_hot.
//                Using T_total overcorrects (recovery factor r < 1, so T_aw < T0).
//   T_coolant  : coolant supply temperature [K]
//   h_hot      : hot-side HTC [W/(m^2*K)]
//   h_coolant  : coolant-side HTC [W/(m^2*K)]
//   eta        : adiabatic cooling effectiveness [-]
//   t_wall     : wall thickness [m]
//   k_wall     : wall thermal conductivity [W/(m*K)]
//
// Returns: wall heat flux [W/m^2]
double cooled_wall_heat_flux(double T_hot, double T_coolant,
                             double h_hot, double h_coolant,
                             double eta,
                             double t_wall, double k_wall);


// -------------------------------------------------------------
// Baldauf et al. (2002) - laterally averaged film-cooling effectiveness
// -------------------------------------------------------------
//
// Baldauf, S., Scheurlen, M., Schulz, A. and Wittig, S. (2002). "Correlation
// of Film-Cooling Effectiveness From Thermographic Measurements at Enginelike
// Conditions." ASME J. Turbomachinery 124(4), 686-698. DOI 10.1115/1.1504443
// (conference version ASME GT-2002-30180).
// docs/heat_transfer/film/686_1_Baldauf_Film_2002.pdf (gitignored).
//
// A row of CYLINDRICAL, streamwise-inclined holes on a flat plate. Unlike the
// older exponential-decay correlations it is valid from the ejection point
// rather than only far downstream, and it carries the ADJACENT JET
// INTERACTION -- the lateral hole-spacing effect that drives jet lift-off --
// as a correlated parameter rather than an excluded case.
//
// eta here is the paper's Eq. (1), (T_G - T_AW)/(T_G - T_C), which is exactly
// combaero's own convention, so it feeds adiabatic_wall_temperature()
// unchanged.
//
// STRUCTURE. The correlation collapses every measured curve onto one base
// curve in scaled coordinates, then back-transforms:
//
//   x/D -> xi'      Eq. (40) inverted in closed form
//   xi' -> eta*'    Eq. (36), the turbulence-dependent base curve
//   eta*' -> eta*   Eq. (37), undoing the adjacent-jet downstream branch
//   eta* -> eta_c   Eq. (38), blowing up to this case's peak
//   eta_c -> eta    Eq. (39), backscaling out the geometry normalisation
//
// ANGLE UNITS. The paper tabulates alpha in DEGREES but every trigonometric
// function consumes RADIANS. This was not assumed: seven of the worked
// example's coefficients, spanning cos(a), cos(1.5a), cos(2.3a), cos(2.5a),
// sin(2a) and cos^0.65(a), reproduce the paper's Table 4 on radians and fail
// on degrees by 10-30%. The API therefore takes degrees and converts.
//
// A PUBLISHED INCONSISTENCY, IMPLEMENTED AS PRINTED. Eq. (31) for b_0 does
// not reproduce the paper's own Table 4:
//
//     Eq. (31) as printed -> 0.83612467
//     Table 4             -> 0.61626073        (36% apart)
//
// Every other one of the 17 checkable Table 4 coefficients reproduces to
// better than 2e-6 relative, so this is the paper's, not ours. Eq. (31) was
// confirmed from four independent channels (publisher MathML, the page at
// 600 dpi, the PDF text layer, and a reader's own transcription); no angle
// convention, unit choice or single-token variant reaches Table 4's value,
// and Table 4 is self-consistent with Eq. (32) to the digit. We implement the
// PRINTED equation, because fidelity means implementing what the paper
// states, and `b_0_override` lets a caller test the Table 4 reading without
// the library inventing a second formula -- the paper contains none. Table 4
// gives one number at the Table 3 conditions; reaching it needs
// sin(...) = -0.167 where the printed equation gives +0.470, and no sign,
// unit or angle convention makes that argument negative at alpha = 30 deg.
//
// THE EFFECT IS GOVERNED BY x/D, NOT BY M. An earlier note here said "under
// 5% below M ~ 0.5 and up to 50% at M = 2.5", which misattributes it: b_0
// feeds b_1 through Eq. (32), and b_1 is the DESCENDING-branch gradient, so
// it does nothing near the hole whatever the blowing. Measured by sweeping
// the envelope with b_0_override:
//
//   x/D <= 20            under 1% at every M
//   x/D = 200, M = 2.0   -19%
//   x/D = 400, M = 2.0   -28%
//   worst over the box   -75%  (s/D 5, alpha 90, Tu 0.0035, M 2.5, x/D 400)
//
// Practical consequence: at effusion row spacings (x/D of 6 to 9) the two
// readings differ by under 1%, so no effusion dataset can arbitrate Eq. (31).
// Only far-downstream single-row data could. See
// validation/cooling/extractions/baldauf_2002_film_effectiveness.md.
namespace baldauf2002 {

// Table 2 -- constants of the base curve fit. Valid for all geometry and
// density ratios, which is the whole point of the base-curve collapse.
constexpr double xi_0 = 9.0;
constexpr double eta_0 = 5.8;      // superseded by eta_0T, Eq. (35)
constexpr double a_star = 4.0;     // ascending-branch gradient
constexpr double b_star = 0.7;     // descending-branch gradient (see Eq. 34)
constexpr double c_star = 0.24;    // apex sharpness, fits the apex to eta* = 1

// Eq. (35) rebases eta_0 when the turbulence-dependent slope b*_T moves.
constexpr double eta_0T_c0 = 2.5;
constexpr double eta_0T_c1 = 5.8;

// The paper's stated envelope. Outside it the correlation is extrapolation:
// Tu in particular appears inside exp[2.6 Tu - 0.0012/Tu^2 - 1.76], which
// diverges as Tu -> 0.
constexpr double M_min = 0.2, M_max = 2.5;
constexpr double P_min = 1.2, P_max = 1.8;
constexpr double s_over_D_min = 2.0, s_over_D_max = 5.0;
constexpr double alpha_deg_min = 30.0, alpha_deg_max = 90.0;
constexpr double Tu_min = 0.0035, Tu_max = 0.075;

// The paper's own accuracy: overall RMS deviation of the measurements about
// the correlation, 5% at the apex and 3% on the descending branch.
constexpr double rms_deviation = 0.055;

// Eq. (31)'s b_0 as the worked example's Table 4 reports it, at the Table 3
// conditions. Kept so the discrepancy above is on the record in code, not
// only in prose. NOT used by the correlation.
constexpr double b0_table4_at_table3 = 0.61626073;

} // namespace baldauf2002

// Laterally averaged adiabatic film-cooling effectiveness downstream of one
// row of cylindrical, streamwise-inclined holes.
//
// Parameters:
//   x_over_D  : streamwise distance from the ejection point, in hole diameters
//   M         : blowing rate (rho u)_C / (rho u)_G [-]
//   P         : density ratio rho_C / rho_G [-]
//   alpha_deg : ejection angle to the SURFACE [deg], 30 to 90
//   s_over_D  : lateral hole spacing / diameter [-], 2 to 5
//   Tu        : mainstream turbulence intensity [-], e.g. 0.015 for 1.5%
//
// Returns: laterally averaged eta [-], Eq. (1) convention.
double film_effectiveness_baldauf_2002(double x_over_D, double M, double P,
                                       double alpha_deg, double s_over_D,
                                       double Tu,
                                       CorrelationStatus *status = nullptr,
                                       double b_0_override =
                                           std::numeric_limits<double>::quiet_NaN());

// Solver-facing (f, J): (eta, d eta/dM, d eta/dP).
//
// M and P are the two inputs a network solve varies -- M through the coolant
// mass flow, P through the temperature ratio. Geometry and Tu are fixed per
// element and carry no partials. Analytic via forward-mode dual numbers, not
// finite differences.
std::tuple<double, double, double> film_effectiveness_baldauf_2002_and_derivatives(
    double x_over_D, double M, double P, double alpha_deg, double s_over_D,
    double Tu, CorrelationStatus *status = nullptr,
    double b_0_override = std::numeric_limits<double>::quiet_NaN());

// -------------------------------------------------------------
// Multi-row film superposition
// -------------------------------------------------------------
//
// Sellers, J.P. (1963). "Gaseous film cooling with multiple slot injection."
// AIAA Journal 1(9), 2154-2156. As restated in:
//
//   Gao, Z., Qiu, T., Liu, P., Ding, S., Li, Z., Cheng, R. and Yuan, Q.
//   (2025). "A Study on the Film Superposition Method for the Multi-Row Film
//   Cooling of the Turbine Outer Ring." Processes 13, 143, Eqs. (1)-(11).
//   docs/heat_transfer/film/processes-13-00143-v2.pdf (open access).
//
// WHY A CORRECTION IS NEEDED. Sellers assumes the rows are independent. Gao
// measures what that costs: "the Sellers method accumulates prediction
// errors as the number of hole rows increases, leading to an OVERESTIMATION
// of the cooling efficiency." That is worst exactly where effusion lives --
// many closely spaced rows -- so a film module built for a few rows and then
// reused for effusion would be wrong in the regime it is needed most.
//
// Gao's fix is a per-row mainstream temperature correction alpha. It is not
// an invented knob: Eqs. (3) and (4) are an energy balance on the mainstream
// entrained into the boundary layer at each injection, giving
//
//   (T_g - T'_aw)/(T_g - T_aw) = C (m_c/m_g) / (C (m_c/m_g) + 1)
//
// so alpha is the fraction of the film's temperature deficit that survives
// mixing on the way to the next row. alpha = 1 recovers Sellers exactly.
//
// WHAT GAO DOES NOT PUBLISH. Eq. (5) gives alpha's functional form,
//
//   alpha_i = a r / (a r + 1) + b,     r = m_coolant / m_mainstream
//
// but the paper never prints the fitted a and b -- there is no coefficient
// table and no inline value. So the SHAPE is sourced and the CONSTANTS are
// not. alpha is therefore exposed as a caller input defaulting to 1
// (i.e. plain Sellers), with Eq. (5) available for anyone fitting their own
// rig. This is the tuner slot the validation policy reserves: matching a
// specific rig is the user's job and the harness never scores it.
namespace film_superposition {

// alpha = 1 is plain Sellers. Anything below it damps the accumulated
// effectiveness, which is the direction Gao's measurements require.
constexpr double alpha_sellers = 1.0;

} // namespace film_superposition

// Sellers superposition, Gao Eq. (1):
//
//   eta = eta_1 + sum_{i=2..n} eta_i prod_{j<i} (1 - eta_j)
//       = 1 - prod_i (1 - eta_i)
//
// The two forms are algebraically identical; the product form is used
// because it is O(n) and cannot accumulate the pairwise rounding the sum
// form does. Known to OVERESTIMATE as the row count grows.
double film_superposition_sellers(const std::vector<double>& eta_rows);

// (eta, d eta / d eta_i) for the above. d eta/d eta_i = prod_{j != i}
// (1 - eta_j), computed without division so a fully effective row
// (eta_j = 1) does not produce a NaN.
std::pair<double, std::vector<double>> film_superposition_sellers_and_gradient(
    const std::vector<double>& eta_rows);

// Gao Eq. (7), Sellers with the per-row mainstream temperature correction:
//
//   eta = sum_i [ eta_i prod_{j=i..n-1} alpha_j prod_{k=i+1..n} (1 - eta_k) ]
//
// alpha_between_rows has n-1 entries for n rows: alpha_j is the correction
// applied between row j and row j+1. Passing all ones reproduces
// film_superposition_sellers exactly.
double film_superposition_corrected(const std::vector<double>& eta_rows,
                                    const std::vector<double>& alpha_between_rows);

// (eta, d eta / d eta_i) for Eq. (7).
std::pair<double, std::vector<double>> film_superposition_corrected_and_gradient(
    const std::vector<double>& eta_rows,
    const std::vector<double>& alpha_between_rows);

// Gao Eq. (5): the published FORM of the mainstream temperature correction.
// The coefficients are NOT published -- see the note above -- so a and b are
// required arguments rather than defaulted, to stop a made-up number
// acquiring the authority of a default.
//
//   mass_flow_ratio : m_coolant / m_mainstream for this row [-]
double mainstream_temperature_correction(double mass_flow_ratio, double a,
                                         double b);

// Gao Eq. (9): equivalent slot width of a row of holes, s = A_hole / pitch.
// This is what lets one row's measured distribution stand in for another
// spacing, via the X/(M s) scaling of Eq. (8).
double equivalent_slot_width(double hole_area, double pitch);

// Gao Eq. (10): equivalent blowing ratio, M_e = M_0 A_0 / A_e. Normalises a
// hole count and spacing onto the baseline single-row configuration.
double equivalent_blowing_ratio(double M_baseline, double area_baseline,
                                double area_equivalent);

} // namespace combaero::cooling

#endif // COOLING_CORRELATIONS_H
