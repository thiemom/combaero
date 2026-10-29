#ifndef MATERIALS_H
#define MATERIALS_H

#include <string>
#include <vector>
#include <functional>

namespace combaero::materials {

// -------------------------------------------------------------
// Material Property Database
// -------------------------------------------------------------

// Material thermal property data structure
struct ThermalMaterial {
    std::string name;                           // Material identifier
    std::function<double(double)> k_func;       // k(T) function [W/(m*K)]
    double T_min;                               // Minimum valid temperature [K]
    double T_max;                               // Maximum valid temperature [K]
    std::string source;                         // Data source reference
    std::string notes;                          // Additional information
};

// Database access functions
const ThermalMaterial& get_material(const std::string& name);
std::vector<std::string> list_materials();

// -------------------------------------------------------------
// Superalloys - Thermal Conductivity k(T)
// -------------------------------------------------------------
// All functions take temperature in Kelvin and return k in W/(m*K)

// Inconel 718 (Ni-based superalloy)
// Source: Haynes International, Special Metals Corporation
// Valid: 300-1200 K
// Linear fit: k = 11.4 + 0.0146*T_C
double k_inconel718(double T);

// Haynes 230 (Ni-Cr-W-Mo alloy)
// Source: HAYNES 230 alloy brochure H-3000O (Haynes International, 2021),
// Physical Properties table: 8.9 W/m K at RT, 10.4 at 100 C rising 2.0 per
// 100 C to 28.4 at 1000 C.
// Valid: 300-1273 K, where the table ends
// Linear: k = 8.4 + 0.020*T_C, exact at every printed point (#761)
double k_haynes230(double T);

// -------------------------------------------------------------
// Structural Alloys - Thermal Conductivity k(T)
// -------------------------------------------------------------

// Stainless Steel 316 (austenitic)
// Source: NIST, ASM Handbook
// Valid: 300-1200 K
// Quadratic fit: k = 13.5 + 0.015*T_C - 2.1e-6*T_C^2
double k_stainless_steel_316(double T);

// Aluminum 6061-T6
// Source: ASM Handbook, Aluminum Association
// Valid: 200-600 K (limited by T6 temper stability)
// Linear fit: k = 167.0 - 0.012*T_C
double k_aluminum_6061(double T);

// -------------------------------------------------------------
// Thermal Barrier Coatings - k(T) with Aging
// -------------------------------------------------------------

// Yttria-Stabilized Zirconia (YSZ) - 7-8 wt% Y2O3
// Source: Zhu, D., Miller, R. A. (2000). NASA/TM-2000-210237 (#761).
// Valid: 300-1700 K
//
// Plasma-sprayed (default): as sprayed k0 = 1.0 W/m K, flat in T (Zhu &
// Miller: "about 1 W/m-K"). Sintering after `hours` at T follows their Eq. 1:
//   k = k0 + (k_inf - k0) A exp(-Q/RT) [1 - exp(-t/tau)],
//   tau = tau0 exp(Q_tau/RT), t in seconds,
// A = 102.2, Q = 68228 J/mol, tau0 = 572.5 s, Q_tau = 41710 J/mol. k_inf =
// 1.74 W/m K is inferred: with Eq. 1's own prefactor it reproduces the three
// plateaus of their Fig. 4 (1.44 / 1.17 / 1.09 W/m K at 1320 / 1100 / 990 C
// surface; read off the figure, unconfirmed). The fit is the overall
// conductivity of a 254 um coating under a 64 W/cm^2 gradient: an effective
// through-thickness value, best for coatings of that order. A user with a
// characterised coating supplies its own conductivity.
//
// EB-PVD: as-deposited 1.5 + 0.0002 T_C only, uncited (kept until a public
// source is held); no sintering.
//
// Parameters:
//   T        : (surface) temperature [K]
//   hours    : exposure at T [h] (default: 0 = as sprayed)
//   is_ebpvd : true for EB-PVD, false for plasma-sprayed (default: false)
// Returns: thermal conductivity k [W/(m*K)]
constexpr double YSZ_T_MIN = 300.0;
constexpr double YSZ_T_MAX = 1700.0;
constexpr double YSZ_APS_K0 = 1.0;        // W/(m K), as sprayed
constexpr double YSZ_APS_K_INF = 1.74;    // W/(m K), inferred from Fig. 4
constexpr double YSZ_ZHU_A = 102.2;       // -
constexpr double YSZ_ZHU_Q = 68228.0;     // J/mol
constexpr double YSZ_ZHU_TAU0_S = 572.5;  // s
constexpr double YSZ_ZHU_Q_TAU = 41710.0; // J/mol
double k_tbc_ysz(double T, double hours = 0.0, bool is_ebpvd = false);

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

// Validate temperature is within material's valid range
// Throws std::runtime_error if out of range
void validate_temperature(const std::string& material_name, double T, double T_min, double T_max);

} // namespace combaero::materials

#endif // MATERIALS_H
