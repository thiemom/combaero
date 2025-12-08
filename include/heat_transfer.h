#ifndef HEAT_TRANSFER_H
#define HEAT_TRANSFER_H

// -------------------------------------------------------------
// Heat Transfer Correlations for Forced Convection
// -------------------------------------------------------------
// All Nusselt correlations return Nu [-], from which:
//   h = Nu * k / L   [W/(m²·K)]
// where k is thermal conductivity [W/(m·K)] and L is characteristic length [m].
//
// For pipe flow: L = D (diameter)
// For flat plate: L = x (distance from leading edge) or plate length
//
// References:
// - Dittus-Boelter (1930): Univ. California Publ. Eng., 2, 443
// - Gnielinski (1976): Int. Chem. Eng., 16, 359
// - Sieder-Tate (1936): Ind. Eng. Chem., 28, 1429
// - Petukhov (1970): Advances in Heat Transfer, 6, 503

// -------------------------------------------------------------
// Internal Flow (Pipe/Duct)
// -------------------------------------------------------------

// Dittus-Boelter correlation (1930)
// Nu = 0.023 * Re^0.8 * Pr^n
// where n = 0.4 for heating (fluid being heated)
//       n = 0.3 for cooling (fluid being cooled)
//
// Valid for:
//   - Fully developed turbulent flow
//   - Re > 10,000
//   - 0.6 < Pr < 160
//   - L/D > 10 (entrance effects negligible)
//
// Parameters:
//   Re      : Reynolds number [-]
//   Pr      : Prandtl number [-]
//   heating : true if fluid is being heated, false if cooled
double nusselt_dittus_boelter(double Re, double Pr, bool heating = true);

// Gnielinski correlation (1976)
// Nu = (f/8) * (Re - 1000) * Pr / (1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1))
//
// More accurate than Dittus-Boelter, especially in transition region.
// Valid for:
//   - 2300 < Re < 5×10^6
//   - 0.5 < Pr < 2000
//
// Parameters:
//   Re  : Reynolds number [-]
//   Pr  : Prandtl number [-]
//   f   : Darcy friction factor [-] (use friction_colebrook or similar)
double nusselt_gnielinski(double Re, double Pr, double f);

// Gnielinski with automatic friction factor (smooth pipe)
// Uses Petukhov friction correlation: f = (0.790*ln(Re) - 1.64)^(-2)
double nusselt_gnielinski(double Re, double Pr);

// Sieder-Tate correlation (1936)
// Nu = 0.027 * Re^0.8 * Pr^(1/3) * (μ_bulk / μ_wall)^0.14
//
// Accounts for viscosity variation between bulk and wall temperatures.
// Valid for:
//   - Re > 10,000
//   - 0.7 < Pr < 16,700
//   - L/D > 10
//
// Parameters:
//   Re       : Reynolds number at bulk temperature [-]
//   Pr       : Prandtl number at bulk temperature [-]
//   mu_ratio : μ_bulk / μ_wall [-] (typically 0.5-2.0)
double nusselt_sieder_tate(double Re, double Pr, double mu_ratio = 1.0);

// Petukhov correlation (1970)
// Nu = (f/8) * Re * Pr / (1.07 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1))
//
// Basis for Gnielinski; valid for fully turbulent flow.
// Valid for:
//   - 10^4 < Re < 5×10^6
//   - 0.5 < Pr < 2000
//
// Parameters:
//   Re  : Reynolds number [-]
//   Pr  : Prandtl number [-]
//   f   : Darcy friction factor [-]
double nusselt_petukhov(double Re, double Pr, double f);

// Petukhov with automatic friction factor (smooth pipe)
double nusselt_petukhov(double Re, double Pr);

// -------------------------------------------------------------
// Laminar Flow (for completeness)
// -------------------------------------------------------------

// Fully developed laminar flow in circular pipe
// Nu = 3.66 (constant wall temperature)
// Nu = 4.36 (constant heat flux)
constexpr double NU_LAMINAR_CONST_T = 3.66;
constexpr double NU_LAMINAR_CONST_Q = 4.36;

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

// Heat transfer coefficient from Nusselt number
// h = Nu * k / L  [W/(m²·K)]
//
// Parameters:
//   Nu : Nusselt number [-]
//   k  : thermal conductivity [W/(m·K)]
//   L  : characteristic length [m] (diameter for pipe flow)
double htc_from_nusselt(double Nu, double k, double L);

// Petukhov friction factor for smooth pipes
// f = (0.790 * ln(Re) - 1.64)^(-2)
// Valid for 3000 < Re < 5×10^6
double friction_petukhov(double Re);

// -------------------------------------------------------------
// State-based convenience functions
// -------------------------------------------------------------

struct State;  // Forward declaration

// Nusselt number for pipe flow using State
// Automatically computes Re, Pr from state and velocity/diameter
double nusselt_pipe(const State& s, double velocity, double diameter,
                    bool heating = true, double roughness = 0.0);

// Heat transfer coefficient for pipe flow [W/(m²·K)]
double htc_pipe(const State& s, double velocity, double diameter,
                bool heating = true, double roughness = 0.0);

#endif // HEAT_TRANSFER_H
