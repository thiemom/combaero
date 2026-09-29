#ifndef COOLING_CORRELATIONS_H
#define COOLING_CORRELATIONS_H

#include <tuple>
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
// states. The effect is under 5% below M ~ 0.5 and up to 50% at M = 2.5; see
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
                                       double Tu);

// Solver-facing (f, J): (eta, d eta/dM, d eta/dP).
//
// M and P are the two inputs a network solve varies -- M through the coolant
// mass flow, P through the temperature ratio. Geometry and Tu are fixed per
// element and carry no partials. Analytic via forward-mode dual numbers, not
// finite differences.
std::tuple<double, double, double> film_effectiveness_baldauf_2002_and_derivatives(
    double x_over_D, double M, double P, double alpha_deg, double s_over_D,
    double Tu);

} // namespace combaero::cooling

#endif // COOLING_CORRELATIONS_H
