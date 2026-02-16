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

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

// Validate parameter ranges for correlations
void validate_rib_params(double e_D, double P_e, double alpha);
void validate_impingement_params(double Re_jet, double z_D, double x_D, double y_D);
void validate_film_cooling_params(double M, double DR, double alpha_deg);

} // namespace combaero::cooling

#endif // COOLING_CORRELATIONS_H
