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
    std::function<double(double)> k_func;       // k(T) function [W/(m·K)]
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
// All functions take temperature in Kelvin and return k in W/(m·K)

// Inconel 718 (Ni-based superalloy)
// Source: Haynes International, Special Metals Corporation
// Valid: 300-1200 K
// Linear fit: k = 11.4 + 0.0146*T_C
double k_inconel718(double T);

// Haynes 230 (Ni-Cr-W-Mo alloy)
// Source: Haynes International Technical Data
// Valid: 300-1400 K
// Quadratic fit: k = 11.8 + 0.0158*T_C + 1.1e-6*T_C^2
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
// Source: NASA TM-2010-216765 (Zhu/Miller sintering model)
// Valid: 300-1700 K
//
// Sintering model:
// - Initial k increases linearly with T
// - Over time, pores close (sintering) -> k increases
// - Rate depends on temperature (Arrhenius)
//
// Parameters:
//   T     : temperature [K]
//   hours : operating hours at temperature (default: 0 = as-sprayed)
//
// Returns: k [W/(m·K)]
//   - As-sprayed: ~0.8-1.5 W/(m·K)
//   - Fully sintered: ~1.8-2.0 W/(m·K)
double k_tbc_ysz(double T, double hours = 0.0);

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

// Validate temperature is within material's valid range
// Throws std::runtime_error if out of range
void validate_temperature(const std::string& material_name, double T, double T_min, double T_max);

} // namespace combaero::materials

#endif // MATERIALS_H
