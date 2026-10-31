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

// k(T) of a database material [W/(m*K)]. Inside [T_min, T_max] this is
// k_func(T). Outside it, extrapolate = false holds the end value (T is
// clamped to the range) and extrapolate = true continues the end slope
// linearly: an opt-in, unsourced continuation the caller asks for. Throws
// for a non-finite T, T <= 0 or a continued k <= 0.
double material_conductivity(const std::string& name, double T, bool extrapolate);

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
// Physical Properties: k = alpha * cp * rho from its thermal-diffusivity and
// specific-heat columns and its density (8.97 g/cm3), 25-1000 C (#881). The
// brochure's printed conductivity column (8.9 + 2.0 per 100 C) is 13% above
// that product at 500 C; Gen3 CSP's measurement (Georgia Tech, 340-780 C)
// agrees with the product within its 95% band and not with the column.
// Valid: 298-1273 K, where the table ends
double k_haynes230(double T);

// Combustor liner / transition-piece alloys. Source: the manufacturers'
// alloy brochures (Haynes International), Physical Properties table, metric
// column, interpolated linearly (exact at every printed point; RT = 25 C).
// Convenience data: a user with qualified properties supplies their own.

// HASTELLOY X (Ni-Cr-Fe-Mo). Valid: 300-1273 K: the brochure table to 900 C
// plus 400 C and 1000 C from the manufacturer's data sheet ("sheet, tested at
// 760 C"), which agrees with the brochure within 0.1 W/m K where both print.
double k_hastelloy_x(double T);

// HAYNES 188 (Co-Ni-Cr-W). Valid: 300-1273 K (table to 1000 C).
double k_haynes188(double T);

// HAYNES 263 (Ni-Co-Cr-Mo, age-hardenable). Valid: 300-1273 K (table to 1000 C).
double k_haynes263(double T);

// Further tabulated materials (#761), interpolated linearly in the source
// table and refused outside it:
// HAYNES 25 (L-605, Co-Cr-W-Ni), brochure table 25-1000 C. Valid: 298-1273 K.
double k_haynes25(double T);
// INCONEL 625 (Ni-Cr-Mo-Nb), Special Metals Table 3, -157 to 982 C. Valid: 116-1255 K.
// Battelle measurements on annealed material. Gen3 CSP (Georgia Tech) measures
// 14% higher at 260-580 C and 25% higher at 620-700 C, from a higher
// diffusivity; unresolved without a third source (#881), so kept.
double k_inconel625(double T);
// NIMONIC 75 (Ni-Cr-Ti), Virgamet N75 table 20-1000 C. Valid: 293-1273 K.
double k_nimonic75(double T);
// High-purity copper, NSRDS-NBS 8 Table 2b. Valid: 100-1200 K.
double k_copper(double T);
// Al2O3, 99.5% pure, 98% dense, NSRDS-NBS 8 Table 14. Valid: 100-2000 K.
double k_alumina(double T);
// X22CrMoV12-1 (1.4923), data sheet values at 20 and 650 C. Valid: 293-923 K.
double k_x22crmov12_1(double T);

// NiCoCrAlY TBC bond coat. Source: Zhu & Miller 1999, NASA/TM-1999-209069,
// Fig. 2 (literature values): k = 6.729572 + 0.02344*T_C. Valid: 848-1373 K,
// the range the line is drawn over (about 575-1100 C).
double k_bondcoat_nicocraly(double T);

// -------------------------------------------------------------
// Structural Alloys - Thermal Conductivity k(T)
// -------------------------------------------------------------

// Stainless Steel 316 (austenitic)
// Source: Georgia Tech Gen3 CSP thermophysical property database, measured
// k (laser-flash alpha * cp * rho), 260-820 C, 95% bands 0.9-2.7 W/m K (#881).
// Valid: 533-1093 K (the measured range)
double k_stainless_steel_316(double T);

// Stainless Steel 310 (austenitic, 25Cr-20Ni)
// Source: Georgia Tech Gen3 CSP thermophysical property database, measured
// k, 260-800 C, 95% bands 0.8-2.5 W/m K (#881).
// Valid: 533-1073 K (the measured range)
double k_stainless_steel_310(double T);

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
// 1.74 W/m K is inferred from the 1320 C plateau of their Fig. 4 (1.44 W/m K).
// Fidelity on Fig. 4 at 32 h (overlay, user-checked 2026-10-29): on the
// 1320 C curve, 2% high at 1100 C (1.19 vs 1.17) and 3% high at 990 C (1.11
// vs 1.08). The printed A and Q give too weak a temperature dependence for
// one k_inf to fit all three, and the paper's own curves are not
// least-squares fits of its data either; Eq. 1 is kept as printed. The fit is the overall
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
