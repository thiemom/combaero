#ifndef COOLING_CORRELATIONS_H
#define COOLING_CORRELATIONS_H

#include <vector>

namespace combaero::cooling {

// -------------------------------------------------------------
// Rib-Enhanced Cooling
// -------------------------------------------------------------

// Rib enhancement factor for heat transfer
// Accounts for increased surface area and turbulence from ribs
//
// Parameters:
//   e_D   : rib height to hydraulic diameter ratio [-]
//   P_e   : rib pitch to rib height ratio [-]
//   alpha : rib angle [degrees]
//
// Returns: Enhancement factor [-] (typically 1.5-4.0)
//
// Source: Han et al. (1988), "Heat Transfer and Friction in Channels
//         with Two Opposite Rib-Roughened Walls"
// Valid: e_D = 0.02-0.1, P_e = 5-20, alpha = 30-90 deg
double rib_enhancement_factor(double e_D, double P_e, double alpha);

// Rib friction factor multiplier
// Accounts for increased pressure drop from ribs
//
// Parameters:
//   e_D   : rib height to hydraulic diameter ratio [-]
//   P_e   : rib pitch to rib height ratio [-]
//
// Returns: Friction multiplier [-] (typically 2-10)
//
// Source: Han et al. (1988)
// Valid: e_D = 0.02-0.1, P_e = 5-20
double rib_friction_multiplier(double e_D, double P_e);

// -------------------------------------------------------------
// Impingement Cooling
// -------------------------------------------------------------

// Impingement jet Nusselt number correlation
// For single jet or jet array impinging on flat plate
//
// Parameters:
//   Re_jet : Reynolds number based on jet diameter [-]
//   Pr     : Prandtl number [-]
//   z_D    : jet-to-plate distance / jet diameter [-]
//   x_D    : streamwise spacing / jet diameter [-] (default: 0 for single jet)
//   y_D    : spanwise spacing / jet diameter [-] (default: 0 for single jet)
//
// Returns: Average Nusselt number [-]
//
// Source: Florschuetz et al. (1981), "Streamwise Flow and Heat Transfer
//         Distributions for Jet Array Impingement with Crossflow"
// Valid: Re_jet = 5000-80000, z_D = 1-12, x_D/y_D = 4-16
double impingement_nusselt(double Re_jet, double Pr, double z_D,
                          double x_D = 0.0, double y_D = 0.0);

// -------------------------------------------------------------
// Film Cooling Effectiveness
// -------------------------------------------------------------

// Adiabatic film cooling effectiveness
// Measures how well coolant film protects surface from hot gas
//
// Parameters:
//   x_D       : downstream distance / hole diameter [-]
//   M         : blowing ratio (rho_c*v_c)/(rho_inf*v_inf) [-]
//   DR        : density ratio rho_c/rho_inf [-]
//   alpha_deg : injection angle [degrees]
//
// Returns: Adiabatic effectiveness eta [-] (0 = no cooling, 1 = perfect)
//
// Source: Baldauf et al. (2002), "Correlation of Film-Cooling Effectiveness
//         from Thermographic Measurements at Enginelike Conditions"
// Valid: M = 0.3-2.5, DR = 1.2-2.0, alpha = 20-90 deg, x_D = 0-40
double film_cooling_effectiveness(double x_D, double M, double DR, double alpha_deg);

// Laterally averaged film cooling effectiveness
// Averaged across span for design calculations
//
// Parameters:
//   x_D       : downstream distance / hole diameter [-]
//   M         : blowing ratio [-]
//   DR        : density ratio [-]
//   alpha_deg : injection angle [degrees]
//   s_D       : hole spacing / hole diameter [-] (default: 3.0)
//
// Returns: Laterally averaged effectiveness eta_avg [-]
//
// Source: Baldauf et al. (2002)
// Valid: M = 0.3-2.5, DR = 1.2-2.0, s_D = 2-6
double film_cooling_effectiveness_avg(double x_D, double M, double DR,
                                     double alpha_deg, double s_D = 3.0);

// Multi-row film cooling using Sellers (1963) superposition
// Combines effectiveness from multiple upstream rows
// Parameters:
//   row_positions_xD : streamwise positions of hole rows [x/D]
//   eval_xD          : evaluation location [x/D]
//   M                : blowing ratio [-]
//   DR               : density ratio [-]
//   alpha_deg        : injection angle [deg]
// Returns: total adiabatic effectiveness [-]
// Reference: Sellers (1963), Baldauf et al. (2002)
// Accuracy: +/- 15-20% (flat plate; curvature requires Ito correction)
double film_cooling_multirow_sellers(
    const std::vector<double>& row_positions_xD,
    double eval_xD,
    double M,
    double DR,
    double alpha_deg
);

// -------------------------------------------------------------
// Effusion Cooling
// -------------------------------------------------------------

// Effusion cooling effectiveness for high-density hole arrays
// Used in combustor liners with continuous coolant injection
//
// Parameters:
//   x_D      : streamwise distance from first row / hole diameter [-]
//   M        : blowing ratio (rho_c*v_c)/(rho_inf*v_inf) [-]
//   DR       : density ratio rho_c/rho_inf [-]
//   porosity : hole area / total surface area [-]
//   s_D      : hole pitch / hole diameter [-]
//   alpha_deg: injection angle [degrees]
//
// Returns: Adiabatic effectiveness eta [-] (0 = no cooling, 1 = perfect)
//
// Source: Lefebvre (1984), "Gas Turbine Combustion"
//         L'Ecuyer & Matsuura (1985), crossflow effects
// Valid: M = 1-4, DR = 1.2-2.0, porosity = 0.02-0.10, s_D = 4-8, alpha = 20-45 deg
// Accuracy: +/- 20-25% (high-density arrays, combustor liner geometry)
double effusion_effectiveness(
    double x_D,
    double M,
    double DR,
    double porosity,
    double s_D,
    double alpha_deg
);

// Effusion hole discharge coefficient
// Accounts for compressibility, inclination, and hole geometry
//
// Parameters:
//   Re_d      : hole Reynolds number based on diameter [-]
//   P_ratio   : plenum/mainstream pressure ratio [-]
//   alpha_deg : hole inclination angle [degrees]
//   L_D       : hole length / diameter [-] (default: 4.0)
//
// Returns: Discharge coefficient Cd [-]
//
// Source: Hay & Spencer (1992), compressible flow through inclined holes
//         Gritsch et al. (1998), shaped hole discharge
// Valid: Re_d > 3000, 1.02 < P_ratio < 1.15, 20 < alpha < 45 deg, 2 < L_D < 8
// Accuracy: +/- 10-15%
double effusion_discharge_coefficient(
    double Re_d,
    double P_ratio,
    double alpha_deg,
    double L_D = 4.0
);

// -------------------------------------------------------------
// Pin Fin Arrays
// -------------------------------------------------------------

// Pin fin array Nusselt number for internal cooling
// Staggered or inline arrangements of cylindrical pins
//
// Parameters:
//   Re_d         : Reynolds number based on pin diameter [-]
//   Pr           : Prandtl number [-]
//   L_D          : pin length / diameter [-]
//   S_D          : spanwise spacing / diameter [-]
//   X_D          : streamwise spacing / diameter [-]
//   is_staggered : true for staggered, false for inline (default: true)
//
// Returns: Average Nusselt number Nu_d [-]
//
// Source: Metzger et al. (1982), "Effects of Pin Shape and Array Orientation"
// Valid: Re_d = 3000-90000, L_D = 0.5-4.0, S_D = 1.5-4.0, X_D = 1.5-4.0
// Accuracy: +/-15-20%
double pin_fin_nusselt(
    double Re_d,
    double Pr,
    double L_D,
    double S_D,
    double X_D,
    bool is_staggered = true
);

// Pin fin array pressure drop friction coefficient
// For use in: dP = N_rows * f_pin * (rho * v_max² / 2)
// where v_max is the velocity at the minimum cross-section.
//
// Parameters:
//   Re_d         : Reynolds number based on pin diameter [-]
//   is_staggered : true for staggered (default), false for inline
//
// Returns: friction coefficient f_pin [-]
//
// Source: Metzger et al. (1982), Simoneau & VanFossen (1984)
// Valid: Re_d = 3000-90000
// Accuracy: +/-20%
double pin_fin_friction(double Re_d, bool is_staggered = true);

// -------------------------------------------------------------
// Dimpled Surfaces
// -------------------------------------------------------------

// Dimpled surface heat transfer enhancement factor
// Semi-spherical indentations for internal cooling
//
// Parameters:
//   Re_Dh : Reynolds number based on channel hydraulic diameter [-]
//   d_Dh  : dimple diameter / channel height [-]
//   h_d   : dimple depth / diameter [-]
//   S_d   : dimple spacing / diameter [-]
//
// Returns: Enhancement factor Nu_dimple/Nu_smooth [-]
//
// Source: Chyu et al. (1997), "Heat Transfer of Arrays of Semi-Spherical Indentations"
// Valid: Re_Dh = 10000-80000, d_Dh = 0.1-0.3, h_d = 0.1-0.3, S_d = 1.5-3.0
// Accuracy: +/-20%
double dimple_nusselt_enhancement(
    double Re_Dh,
    double d_Dh,
    double h_d,
    double S_d
);

// Dimpled surface friction multiplier
// Friction penalty for dimpled surfaces
//
// Parameters:
//   Re_Dh : Reynolds number based on channel hydraulic diameter [-]
//   d_Dh  : dimple diameter / channel height [-]
//   h_d   : dimple depth / diameter [-]
//
// Returns: Friction multiplier f_dimple/f_smooth [-]
//
// Source: Chyu et al. (1997)
// Valid: Re_Dh = 10000-80000, d_Dh = 0.1-0.3, h_d = 0.1-0.3
// Accuracy: +/-15%
double dimple_friction_multiplier(
    double Re_Dh,
    double d_Dh,
    double h_d
);

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

void validate_rib_params(double e_D, double P_e);
void validate_rib_params(double e_D, double P_e, double alpha);
void validate_impingement_params(double Re_jet, double z_D, double x_D, double y_D);
void validate_film_cooling_params(double M, double DR, double alpha_deg);
void validate_effusion_params(double M, double DR, double porosity, double s_D, double alpha_deg);
void validate_effusion_discharge_params(double Re_d, double P_ratio, double alpha_deg, double L_D);
void validate_pin_fin_params(double Re_d, double L_D, double S_D, double X_D);
void validate_dimple_params(double Re_Dh, double d_Dh, double h_d);
void validate_dimple_params(double Re_Dh, double d_Dh, double h_d, double S_d);

} // namespace combaero::cooling

#endif // COOLING_CORRELATIONS_H
