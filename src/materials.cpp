#include "materials.h"
#include <cmath>
#include <stdexcept>
#include <map>

namespace combaero::materials {

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

void validate_temperature(const std::string& material_name, double T, double T_min, double T_max) {
    if (T < T_min || T > T_max) {
        throw std::runtime_error(
            "Temperature " + std::to_string(T) + " K is outside valid range [" +
            std::to_string(T_min) + ", " + std::to_string(T_max) + "] K for " + material_name
        );
    }
}

// -------------------------------------------------------------
// Superalloys
// -------------------------------------------------------------

double k_inconel718(double T) {
    constexpr double T_min = 300.0;
    constexpr double T_max = 1200.0;
    validate_temperature("Inconel 718", T, T_min, T_max);
    
    double T_C = T - 273.15;
    return 11.4 + 0.0146 * T_C;
}

double k_haynes230(double T) {
    constexpr double T_min = 300.0;
    constexpr double T_max = 1400.0;
    validate_temperature("Haynes 230", T, T_min, T_max);
    
    double T_C = T - 273.15;
    return 11.8 + 0.0158 * T_C + 1.1e-6 * T_C * T_C;
}

// -------------------------------------------------------------
// Structural Alloys
// -------------------------------------------------------------

double k_stainless_steel_316(double T) {
    constexpr double T_min = 300.0;
    constexpr double T_max = 1200.0;
    validate_temperature("Stainless Steel 316", T, T_min, T_max);
    
    double T_C = T - 273.15;
    return 13.5 + 0.015 * T_C - 2.1e-6 * T_C * T_C;
}

double k_aluminum_6061(double T) {
    constexpr double T_min = 200.0;
    constexpr double T_max = 600.0;
    validate_temperature("Aluminum 6061", T, T_min, T_max);
    
    double T_C = T - 273.15;
    return 167.0 - 0.012 * T_C;
}

// -------------------------------------------------------------
// Thermal Barrier Coatings
// -------------------------------------------------------------

double k_tbc_ysz(double T, double hours, bool is_ebpvd) {
    // Validate temperature
    if (T < 300.0 || T > 1700.0) {
        throw std::runtime_error("k_tbc_ysz: Temperature must be in range [300, 1700] K, got " + std::to_string(T));
    }
    
    double T_C = T - 273.15;
    
    // As-sprayed thermal conductivity (different for APS vs EB-PVD)
    double k_initial;
    if (is_ebpvd) {
        // EB-PVD: Electron beam physical vapor deposition
        // Columnar microstructure -> higher initial conductivity
        k_initial = 1.5 + 0.0002 * T_C;
    } else {
        // APS: Atmospheric plasma spray (default)
        // Splat boundaries -> lower initial conductivity
        k_initial = 0.8 + 0.00045 * T_C;
    }
    
    // No aging or low temperature (sintering negligible below ~1073 K / 800 C)
    // Sintering is diffusion-driven and requires high temperature
    if (hours <= 0.0 || T < 1073.0) {
        return k_initial;
    }
    
    // Sintering model (NASA Zhu/Miller TM-2010-216765)
    // k increases as pores close over time at high temperature
    constexpr double k_fully_sintered = 1.85;  // W/(m·K) - fully dense YSZ limit
    constexpr double E_activation = 12000.0;   // Activation energy [K]
    
    // Sintering rate (Arrhenius relationship)
    // Rate increases exponentially with temperature
    double rate = 0.01 * std::exp(-E_activation / T);
    
    // Sintering progress (stretched exponential, n=0.4 for APS YSZ)
    // zeta = 1 - exp(-(rate*t)^n) where 0 <= zeta <= 1
    double progress = 1.0 - std::exp(-std::pow(rate * hours, 0.4));
    
    // Interpolate between as-sprayed and fully sintered
    double k_sintered = k_initial + (k_fully_sintered - k_initial) * progress;
    
    // Clamp to physical limit (fully dense YSZ)
    return std::min(k_sintered, k_fully_sintered);
}

// -------------------------------------------------------------
// Material Database
// -------------------------------------------------------------

static const std::map<std::string, ThermalMaterial> material_database = {
    {"inconel718", {
        "inconel718",
        k_inconel718,
        300.0, 1200.0,
        "Haynes International, Special Metals Corporation",
        "Ni-based superalloy, common in turbine applications"
    }},
    {"haynes230", {
        "haynes230",
        k_haynes230,
        300.0, 1400.0,
        "Haynes International Technical Data",
        "Ni-Cr-W-Mo alloy, high-temperature oxidation resistance"
    }},
    {"ss316", {
        "ss316",
        k_stainless_steel_316,
        300.0, 1200.0,
        "NIST, ASM Handbook",
        "Austenitic stainless steel, corrosion resistant"
    }},
    {"al6061", {
        "al6061",
        k_aluminum_6061,
        200.0, 600.0,
        "ASM Handbook, Aluminum Association",
        "Aluminum alloy, limited to <300°C due to T6 temper"
    }},
    {"ysz", {
        "ysz",
        [](double T) { return k_tbc_ysz(T, 0.0); },  // Default: as-sprayed
        300.0, 1700.0,
        "NASA TM-2010-216765 (Zhu/Miller sintering model)",
        "7-8 wt% Y2O3 stabilized zirconia, thermal barrier coating"
    }}
};

const ThermalMaterial& get_material(const std::string& name) {
    auto it = material_database.find(name);
    if (it == material_database.end()) {
        throw std::runtime_error("Material '" + name + "' not found in database");
    }
    return it->second;
}

std::vector<std::string> list_materials() {
    std::vector<std::string> names;
    names.reserve(material_database.size());
    for (const auto& pair : material_database) {
        names.push_back(pair.first);
    }
    return names;
}

} // namespace combaero::materials
