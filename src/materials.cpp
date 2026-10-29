#include "materials.h"
#include "input_validation.h"
#include "thermo.h"
#include <cmath>
#include <stdexcept>
#include <map>
#include <string>
#include <vector>

namespace combaero::materials {

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

void validate_temperature(const std::string& material_name, double T, double T_min, double T_max) {
    if (!validation::is_finite_in_range(T, T_min, T_max)) {
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
    constexpr double T_max = 1000.0 + 273.15;
    validate_temperature("Haynes 230", T, T_min, T_max);

    // HAYNES 230 brochure H-3000O (2021), "Physical Properties": k = 8.9 W/m K
    // at room temperature and 10.4 + 2.0 per 100 C to 28.4 at 1000 C. The
    // line below reproduces every printed value from 25 to 1000 C, where the
    // table (and the valid range) ends. The quadratic used before was 37%
    // high at room temperature (#761).
    double T_C = T - 273.15;
    return 8.4 + 0.020 * T_C;
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
    if (!validation::is_finite_in_range(T, YSZ_T_MIN, YSZ_T_MAX)) {
        throw std::runtime_error("k_tbc_ysz: Temperature must be in range [300, 1700] K, got " + std::to_string(T));
    }
    if (std::isnan(hours)) {
        throw std::runtime_error("k_tbc_ysz: hours must not be NaN");
    }

    if (is_ebpvd) {
        // EB-PVD (columnar): as-deposited only, an uncited legacy fit kept until
        // a public source is held (sources.yaml wanted list). Zhu & Miller's
        // kinetics below are for plasma-sprayed coatings and are not applied.
        return 1.5 + 0.0002 * (T - 273.15);
    }

    // Plasma-sprayed ZrO2-8wt%Y2O3 (Zhu & Miller 2000, NASA/TM-2000-210237).
    // As sprayed about 1 W/m K, flat in T (p. 6); sintering raises it as
    //   (k - k0)/(k_inf - k0) = A exp(-Q / R T) [1 - exp(-t / tau)],
    //   tau = tau0 exp(Q_tau / R T)  [s],
    // their Eq. 1 (p. 7) for the overall conductivity under a 64 W/cm^2
    // gradient at surface temperature T: an effective through-thickness value.
    if (!(hours > 0.0)) {
        return YSZ_APS_K0;
    }
    const double RT = thermo::R_GAS * T;
    const double tau_s = YSZ_ZHU_TAU0_S * std::exp(YSZ_ZHU_Q_TAU / RT);
    const double reach = YSZ_ZHU_A * std::exp(-YSZ_ZHU_Q / RT);
    const double progress = reach * (1.0 - std::exp(-hours * 3600.0 / tau_s));
    return YSZ_APS_K0 + (YSZ_APS_K_INF - YSZ_APS_K0) * progress;
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
        300.0, 1000.0 + 273.15,
        "HAYNES 230 alloy brochure H-3000O (Haynes International, 2021)",
        "Ni-Cr-W-Mo alloy; table to 1000 C"
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
        "Aluminum alloy, limited to <300degC due to T6 temper"
    }},
    {"ysz", {
        "ysz",
        [](double T) { return k_tbc_ysz(T, 0.0); },  // Default: as-sprayed
        300.0, 1700.0,
        "Zhu & Miller 2000, NASA/TM-2000-210237 (plasma-sprayed, as sprayed)",
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
