#ifndef COOLING_CORRELATIONS_H
#define COOLING_CORRELATIONS_H

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


} // namespace combaero::cooling

#endif // COOLING_CORRELATIONS_H
