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
// GAO'S COEFFICIENTS, AND WHY THEY ARE STILL NOT A DEFAULT. Eq. (5) gives
//
//   alpha_i = a r / (a r + 1) + b,     r = m_coolant / m_mainstream
//
// and section 4.3.1 DOES print the fitted values: "The empirical
// coefficients in Equation (5), a and b, were determined to be 12 and
// 0.9465, respectively." They are recorded below as gao2025::a_case1 and
// b_case1. (An earlier note here said they were never published; that was
// wrong, and it was wrong because a regex looked for "a = 12" while the
// paper states it in prose.)
//
// They are NOT defaulted, because they cannot reach what a tighter plate
// needs -- and that follows from the constants alone, with no rig
// arithmetic. Since a r / (a r + 1) >= 0 for a > 0 and r >= 0,
//
//     alpha >= b = 0.9465   for EVERY r,
//
// so the published pair can never damp a row by more than 5.35% while
// staying at or below 1 (above 1 the superposition refuses -- alpha > 1
// would create coolant). Murray & Ireland's 5.75 D staggered plate needs a
// per-row alpha of about 0.85 at low blowing and 0.69 near M = 1, both
// BELOW that floor. No choice of r reaches them.
//
// This argument is deliberately independent of how r is scaled: Gao's test
// section dimensions are not stated in the text, so r cannot be put on
// their scale, and any claim that depended on doing so would be
// unverifiable.
//
// The likely reason is structural: Eq. (5) carries NO streamwise-spacing
// term, while Murray shows streamwise spacing to be the dominant driver --
// tripling it to 17.25 D restores plain superposition to under 10% error.
// Gao calibrated on a plate with 10.5 d streamwise spacing; Murray's
// effective spacing is 2.875 D, 3.7x tighter. Two plates with the same
// coolant fraction and different row spacing get the same alpha from this
// form, which cannot be right if spacing dominates.
//
// So alpha stays a caller input defaulting to 1 (plain Sellers). This is
// the tuner slot the validation policy reserves: matching a specific rig is
// the user's job and the harness never scores it. See
// validation/cooling/extractions/gao_alpha_fit_on_murray.md.
namespace film_superposition {

// alpha = 1 is plain Sellers. Anything below it damps the accumulated
// effectiveness, which is the direction Gao's measurements require.
constexpr double alpha_sellers = 1.0;

// Gao's own fitted coefficients for Eq. (5), section 4.3.1. Recorded so the
// published values are in code rather than only in prose, NOT so they can be
// used as defaults -- see the note above for why they do not transfer.
// Calibrated on Case 1 (d = 1.2 mm, spanwise pitch 3.5 d, streamwise spacing
// 10.5 d) at blowing ratios 0.3 and 1.0.
constexpr double gao_a_case1 = 12.0;
constexpr double gao_b_case1 = 0.9465;

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

// Gao Eq. (5), the mainstream temperature correction. Both the form and the
// coefficients are published (gao_a_case1, gao_b_case1 above), but a and b
// stay REQUIRED rather than defaulted: the published pair is calibrated on
// one plate and returns alpha > 1 on a tighter one -- see the note above.
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

// -------------------------------------------------------------
// Effusion plate INTERNAL heat transfer
// -------------------------------------------------------------
//
// Andrews, Alikhanizadeh, Asere, Hussain, Khoshkbar Azari and Mkpadi (1986),
// "Small Diameter Film Cooling Holes: Wall Convective Heat Transfer",
// ASME 86-GT-225. See
// validation/cooling/extractions/andrews_effusion_internal_h.md.
//
// A coolant hole cools its wall in TWO places, and the paper's central point
// is that the first dominates:
//
//   1. the APPROACH flow over the coolant-side plate surface, converging
//      into the hole (Sparrow's multi-hole correlation, Eq. 15);
//   2. the THROAT, a short tube with a sharp-edged entry (Mills' tabulated
//      data, Eq. 12 with the entry-length factor R_Nu of Eqs. 13 and 14).
//
// Andrews sums them -- "the authors have treated the wall heat transfer as
// the summation of Equations 13 or 14 and 15" -- after putting both on the
// same Nusselt definition, which is Eq. 18's job. Both Nu below are
// therefore referenced to the HOLE DIAMETER and the HOLE INTERNAL AREA
// pi D L, so they add directly.
//
// CONVERTING TO A PLATE-AREA COEFFICIENT. Effusion elements want a
// coefficient per unit plate area, which is what the source's own Fig. 8
// plots. That is
//
//   h_plate = Nu k / D * A_h / A,   A_h = pi D L,  A = X^2 - pi D^2 / 4
//
// and the A_h/A factor is not optional: for Andrews' plate C it is 3.46, so
// omitting it overstates the coefficient by that factor. The conversion is
// left to the caller because it needs no correlation -- only geometry.
namespace andrews1986 {

// Eq. 18's leading constant, 0.881, is Sparrow's multi-hole approach
// correlation Nu_l = 0.881 Re^0.476 Pr^(1/3) rebased from the surface
// dimension l (Eq. 16) onto the hole internal area. Independent of
// pitch-to-diameter ratio, per the paper.
constexpr double sparrow_coefficient = 0.881;
constexpr double sparrow_re_exponent = 0.476;

// Mills' throat data, Eq. 12, is Dittus-Boelter in form with a tabulated
// entry-length multiplier.
constexpr double mills_coefficient = 0.023;
constexpr double mills_re_exponent = 0.8;

// R_Nu is curve-fitted in two branches that meet at L/D = 2. They agree
// there to 1.3e-3, which is how the transcription was checked.
constexpr double mills_branch_L_over_D = 2.0;

} // namespace andrews1986

// Mills' entry-length factor R_Nu for a sharp-edged tube entry, Andrews
// Eqs. (13) and (14). A short hole never reaches a developed profile, so
// R_Nu > 1; it decays to 1 as the hole lengthens.
//
//   L/D <= 2 :  R_Nu = 0.13 z^3 - 0.75 z^2 + 1.04 z + 2.24,       z = L/D
//   L/D >  2 :  R_Nu = 1 - 49.3 w^4 + 58.6 w^3 - 26.5 w^2 + 7.48 w, w = D/L
double mills_entry_length_factor(double L_over_D);

// Sparrow's hole-approach contribution, Andrews Eq. (18): the coolant-side
// plate surface converging into the hole, expressed on the hole's own
// Nusselt definition so it can be summed with the throat term.
//
//   Nu = 0.881 Re^0.476 Pr^(1/3) * X / (pi L)
//
//   Re        : on the HOLE diameter [-]
//   X_over_L  : hole pitch / hole length [-]
double effusion_approach_nusselt(double Re, double Pr, double X_over_L);

// Mills' short-hole throat contribution, Andrews Eq. (12).
//
//   Nu = 0.023 Re^0.8 Pr^(1/3) R_Nu(L/D)
double effusion_throat_nusselt(double Re, double Pr, double L_over_D);

// Andrews Eq. (19): the two summed, both on the hole diameter and hole
// internal area.
//
//   Nu = (0.881 (X/pi L) Re^0.476 + 0.023 Re^0.8 R_Nu) Pr^(1/3)
//
// Multiply by k/D for a hole-area coefficient, and by A_h/A again for the
// plate-area coefficient an effusion element needs.
double effusion_internal_nusselt(double Re, double Pr, double X_over_L,
                                 double L_over_D);

} // namespace combaero::cooling

#endif // COOLING_CORRELATIONS_H
