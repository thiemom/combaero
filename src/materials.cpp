#include "materials.h"
#include "input_validation.h"
#include "thermo.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
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

namespace {

// Linear interpolation in a manufacturer's k(T_C) table: exact at every
// printed point. T_C must lie within the table (the callers validate T).
template <std::size_t N>
double table_k(const std::array<double, N>& T_C, const std::array<double, N>& k, double t) {
    for (std::size_t i = 1; i < N; ++i) {
        if (t <= T_C[i]) {
            const double w = (t - T_C[i - 1]) / (T_C[i] - T_C[i - 1]);
            return k[i - 1] + w * (k[i] - k[i - 1]);
        }
    }
    return k[N - 1];
}

} // namespace

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
    validate_temperature("Haynes 230", T, 25.0 + 273.15, 1000.0 + 273.15);
    // HAYNES 230 brochure H-3000O (2021), Physical Properties: k = alpha cp
    // rho from its diffusivity and specific-heat columns, rho = 8970 kg/m3
    // (#881). The 800 C row is skipped: its specific heat prints 5595 J/kg C
    // (574 at 700 C, 609 at 900 C), a misprint. The step from 600 to 700 C
    // follows the specific heat (486 -> 574); Gen3 CSP measures the same step.
    // The brochure's own conductivity column (8.4 + 0.020 T_C, shipped until
    // #881) is 13% above this product at 500 C and outside Gen3 CSP's band.
    static constexpr std::array<double, 10> T_C = {25.0, 100.0, 200.0, 300.0, 400.0,
                                                   500.0, 600.0, 700.0, 900.0, 1000.0};
    static constexpr std::array<double, 10> k = {8.62, 10.07, 11.67, 13.22, 14.89,
                                                 16.33, 18.27, 22.14, 24.25, 26.68};
    return table_k(T_C, k, T - 273.15);
}

// Combustor liner / transition-piece alloys: the manufacturers' Physical
// Properties tables, metric column, interpolated linearly. "RT" is taken
// as 25 C. Brochures accessed 2026-10-29 (sources.yaml).

double k_hastelloy_x(double T) {
    validate_temperature("Hastelloy X", T, 300.0, 1273.15);
    // HASTELLOY X brochure (25-900 C; its metric rows jump from 200 to 600 C),
    // with 400 C and 1000 C from the manufacturer's data sheet "sheet, tested
    // at 760 C": 16.9 W/m K at 400 C (its "390 F" label repeats the 200 C
    // row above; its 117 Btu-in/ft2-h-F is 16.9 W/m K, and the brochure's
    // British column interpolates to 117 at 752 F) and 28.7 at 1000 C.
    // Where both print a value they agree within 0.1.
    static constexpr std::array<double, 12> T_C = {25.0, 100.0, 200.0, 400.0, 600.0, 650.0,
                                                   700.0, 750.0, 800.0, 850.0, 900.0, 1000.0};
    static constexpr std::array<double, 12> k = {9.2, 11.2, 14.1, 16.9, 20.9, 21.9,
                                                 22.8, 23.8, 24.7, 25.7, 26.7, 28.7};
    return table_k(T_C, k, T - 273.15);
}

double k_haynes188(double T) {
    validate_temperature("Haynes 188", T, 300.0, 1273.15);
    static constexpr std::array<double, 11> T_C = {25.0, 100.0, 200.0, 300.0, 400.0, 500.0,
                                                   600.0, 700.0, 800.0, 900.0, 1000.0};
    static constexpr std::array<double, 11> k = {10.4, 12.2, 14.3, 15.9, 17.5, 19.3,
                                                 21.1, 23.0, 24.8, 25.5, 27.6};
    return table_k(T_C, k, T - 273.15);
}

double k_haynes263(double T) {
    validate_temperature("Haynes 263", T, 300.0, 1273.15);
    static constexpr std::array<double, 11> T_C = {25.0, 100.0, 200.0, 300.0, 400.0, 500.0,
                                                   600.0, 700.0, 800.0, 900.0, 1000.0};
    static constexpr std::array<double, 11> k = {11.7, 13.0, 14.7, 16.3, 18.0, 19.7,
                                                 21.4, 23.0, 24.7, 26.8, 28.5};
    return table_k(T_C, k, T - 273.15);
}

// -------------------------------------------------------------
// Further alloys, copper and alumina (#761): tables interpolated linearly
// -------------------------------------------------------------

double k_haynes25(double T) {
    // HAYNES 25 (L-605) brochure, Haynes International.
    validate_temperature("Haynes 25", T, 25.0 + 273.15, 1000.0 + 273.15);
    static constexpr std::array<double, 11> T_C = {25.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0};
    static constexpr std::array<double, 11> k = {10.5, 12.0, 14.0, 15.9, 17.7, 19.5, 21.2, 22.9, 24.5, 26.0, 27.5};
    return table_k(T_C, k, T - 273.15);
}

double k_inconel625(double T) {
    // INCONEL alloy 625, Special Metals, Table 3 (Battelle measurements,
    // annealed 2100 F / 1 h); the 1093 C row prints no conductivity. Gen3 CSP
    // measures 14% more at 260-580 C (its diffusivity is 22-40% higher, its
    // specific heat lower) and 25% more at 620-700 C, where its specific heat
    // steps up as the alloy ages. Kept until a third source settles it (#881).
    validate_temperature("Inconel 625", T, -157.0 + 273.15, 982.0 + 273.15);
    static constexpr std::array<double, 15> T_C = {-157.0, -129.0, -73.0, -18.0, 21.0, 38.0, 93.0, 204.0, 316.0, 427.0, 538.0, 649.0, 760.0, 871.0, 982.0};
    static constexpr std::array<double, 15> k = {7.2, 7.5, 8.4, 9.2, 9.8, 10.1, 10.8, 12.5, 14.1, 15.7, 17.5, 19.0, 20.8, 22.8, 25.2};
    return table_k(T_C, k, T - 273.15);
}

double k_nimonic75(double T) {
    // Virgamet data page for Alloy N75 / 2.4951; the 20 C value equals
    // the Special Metals NIMONIC 75 sheet (11.7), which prints no other.
    validate_temperature("Nimonic 75", T, 20.0 + 273.15, 1000.0 + 273.15);
    static constexpr std::array<double, 7> T_C = {20.0, 400.0, 600.0, 700.0, 800.0, 900.0, 1000.0};
    static constexpr std::array<double, 7> k = {11.7, 18.6, 22.7, 24.7, 26.5, 28.4, 30.1};
    return table_k(T_C, k, T - 273.15);
}

double k_copper(double T) {
    // High-purity copper, NSRDS-NBS 8 (Powell, Ho & Liley 1966) Table 2b;
    // printed in W/(cm K), x100 here; the extrapolated 1300 K value is not used.
    validate_temperature("Copper", T, 100.0, 1200.0);
    static constexpr std::array<double, 16> T_K = {100.0, 150.0, 200.0, 250.0, 273.0, 300.0, 350.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0};
    static constexpr std::array<double, 16> k = {483.0, 428.0, 413.0, 404.0, 401.0, 398.0, 394.0, 392.0, 388.0, 383.0, 377.0, 371.0, 364.0, 357.0, 350.0, 342.0};
    return table_k(T_K, k, T);
}

double k_alumina(double T) {
    // Al2O3 99.5% pure, 98% dense, polycrystalline: NSRDS-NBS 8 Table 14;
    // printed in W/(cm K), x100 here; the extrapolated 2100 K value is not used.
    validate_temperature("Alumina", T, 100.0, 2000.0);
    static constexpr std::array<double, 24> T_K = {100.0, 150.0, 200.0, 250.0, 273.0, 300.0, 350.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0};
    static constexpr std::array<double, 24> k = {133.0, 77.0, 55.0, 43.4, 39.7, 36.0, 30.7, 26.4, 20.2, 15.8, 12.6, 10.4, 8.9, 7.85, 7.1, 6.55, 6.13, 5.85, 5.66, 5.56, 5.54, 5.59, 5.74, 6.0};
    return table_k(T_K, k, T);
}

double k_x22crmov12_1(double T) {
    // Edelstahlwerke Suedwestfalen 1.4923 data sheet: 24 W/m K at 20 C and
    // 29 at 650 C, the only two values it prints; linear between them.
    validate_temperature("X22CrMoV12-1 (1.4923)", T, 20.0 + 273.15, 650.0 + 273.15);
    static constexpr std::array<double, 2> T_C = {20.0, 650.0};
    static constexpr std::array<double, 2> k = {24.0, 29.0};
    return table_k(T_C, k, T - 273.15);
}

// -------------------------------------------------------------
// Bond coat
// -------------------------------------------------------------

double k_bondcoat_nicocraly(double T) {
    // Zhu & Miller 1999, NASA/TM-1999-209069, Fig. 2 (literature values):
    // k = 6.729572 + 0.02344 T_C, the line drawn over about 575-1100 C.
    validate_temperature("NiCoCrAlY bond coat", T, 848.15, 1373.15);
    return 6.729572 + 0.02344 * (T - 273.15);
}

// -------------------------------------------------------------
// Structural Alloys
// -------------------------------------------------------------

double k_stainless_steel_316(double T) {
    // Georgia Tech Gen3 CSP database, SS316 measured k (laser-flash alpha
    // times cp times a constant rho), 95% band 0.9 W/m K at 260 C to 2.7 at
    // 820 C (#881). The quadratic shipped before (unsourced) lay inside the
    // band at every point but ran to 1200 K, past the data.
    validate_temperature("Stainless Steel 316", T, 260.0 + 273.15, 820.0 + 273.15);
    static constexpr std::array<double, 15> T_C = {260.0, 300.0, 340.0, 380.0, 420.0,
                                                   460.0, 500.0, 540.0, 580.0, 620.0,
                                                   660.0, 700.0, 740.0, 780.0, 820.0};
    static constexpr std::array<double, 15> k = {17.83, 18.32, 18.90, 19.37, 19.83,
                                                 20.33, 20.78, 21.16, 21.78, 22.50,
                                                 22.79, 22.97, 22.87, 22.90, 23.12};
    return table_k(T_C, k, T - 273.15);
}

double k_stainless_steel_310(double T) {
    // Georgia Tech Gen3 CSP database, SS310 measured k, 95% band 0.8 W/m K
    // at 260 C to 2.5 at 800 C (#881). Not joined to NIST's cryogenic fit
    // (11.6 at 300 K): this data's trend runs back to about 13.7 there, so
    // the two are different samples, and the gap is not bridged.
    validate_temperature("Stainless Steel 310", T, 260.0 + 273.15, 800.0 + 273.15);
    static constexpr std::array<double, 15> T_C = {260.0, 300.0, 340.0, 380.0, 420.0,
                                                   460.0, 500.0, 540.0, 580.0, 620.0,
                                                   660.0, 700.0, 740.0, 780.0, 800.0};
    static constexpr std::array<double, 15> k = {17.13, 17.72, 18.31, 18.91, 19.47,
                                                 20.04, 20.55, 21.11, 22.08, 23.08,
                                                 23.58, 24.13, 24.55, 24.83, 25.12};
    return table_k(T_C, k, T - 273.15);
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
    // Fidelity on their Fig. 4 at 32 h: on the 1320 C curve, +2% at 1100 C
    // and +3% at 990 C (the printed A, Q give too weak a T dependence; the
    // paper's own curves are not least-squares fits either). Kept as printed.
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

// No held source prints this material's solidus (#882).
constexpr double NO_SOLIDUS = std::numeric_limits<double>::quiet_NaN();

static const std::map<std::string, ThermalMaterial> material_database = {
    {"inconel718", {
        "inconel718",
        k_inconel718,
        300.0, 1200.0,
        "Haynes International, Special Metals Corporation",
        "Ni-based superalloy, common in turbine applications",
        1260.0 + 273.15,
        "INCONEL alloy 718 (Special Metals) Table 2, melting range 1260-1336 C"
    }},
    {"haynes230", {
        "haynes230",
        k_haynes230,
        25.0 + 273.15, 1000.0 + 273.15,
        "HAYNES 230 brochure H-3000O (Haynes International, 2021): alpha * cp * rho, 25-1000 C",
        "Ni-Cr-W-Mo alloy; table to 1000 C",
        1301.0 + 273.15,
        "HAYNES 230 brochure H-3000O, melting 1301-1371 C"
    }},
    {"hastelloy_x", {
        "hastelloy_x",
        k_hastelloy_x,
        300.0, 1273.15,
        "HASTELLOY X brochure + data sheet (Haynes International), table to 1000 C",
        "Ni-Cr-Fe-Mo combustor liner alloy",
        1260.0 + 273.15,
        "HASTELLOY X brochure, melting range 1260-1355 C"
    }},
    {"haynes188", {
        "haynes188",
        k_haynes188,
        300.0, 1273.15,
        "HAYNES 188 alloy brochure (Haynes International), table to 1000 C",
        "Co-Ni-Cr-W liner / transition-piece alloy",
        1315.0 + 273.15,
        "HAYNES 188 brochure, melting 1315-1410 C"
    }},
    {"haynes263", {
        "haynes263",
        k_haynes263,
        300.0, 1273.15,
        "HAYNES 263 alloy brochure (Haynes International), table to 1000 C",
        "Ni-Co-Cr-Mo age-hardenable liner / transition-piece alloy",
        1300.0 + 273.15,
        "HAYNES 263 brochure, melting range 1300-1355 C"
    }},
    {"haynes25", {
        "haynes25",
        k_haynes25,
        25.0 + 273.15, 1000.0 + 273.15,
        "HAYNES 25 alloy brochure (Haynes International), table to 1000 C",
        "Co-Cr-W-Ni (L-605) combustor / turbine alloy",
        1330.0 + 273.15,
        "HAYNES 25 brochure, melting range 1330-1410 C"
    }},
    {"inconel625", {
        "inconel625",
        k_inconel625,
        -157.0 + 273.15, 982.0 + 273.15,
        "INCONEL alloy 625 (Special Metals) Table 3, -157 to 982 C",
        "Ni-Cr-Mo-Nb liner / duct / bellows alloy",
        1290.0 + 273.15,
        "INCONEL alloy 625 (Special Metals) Table 2, melting range 1290-1350 C"
    }},
    {"nimonic75", {
        "nimonic75",
        k_nimonic75,
        20.0 + 273.15, 1000.0 + 273.15,
        "Virgamet Alloy N75 data page; 20 C value = Special Metals sheet",
        "Ni-Cr-Ti liner / sheet alloy",
        1340.0 + 273.15,
        "NIMONIC alloy 75 (Special Metals SMC-058), melting range 1340-1380 C"
    }},
    {"copper", {
        "copper",
        k_copper,
        100.0, 1200.0,
        "NSRDS-NBS 8 (Powell, Ho & Liley 1966) Table 2b",
        "High-purity copper; alloyed or cold-worked copper conducts less",
        1083.0 + 273.15,
        "ASM Handbook Vol. 2 (1990), C10100/C10200: solidus 1083 C"
    }},
    {"alumina", {
        "alumina",
        k_alumina,
        100.0, 2000.0,
        "NSRDS-NBS 8 (Powell, Ho & Liley 1966) Table 14",
        "Dense polycrystalline Al2O3 (sheaths, insulators); porous grades conduct less",
        NO_SOLIDUS,
        ""
    }},
    {"x22crmov12_1", {
        "x22crmov12_1",
        k_x22crmov12_1,
        20.0 + 273.15, 650.0 + 273.15,
        "Edelstahlwerke Suedwestfalen 1.4923 data sheet (20 and 650 C)",
        "12% Cr martensitic steel (steam turbine casings, bolts)",
        NO_SOLIDUS,
        ""
    }},
    {"nicocraly", {
        "nicocraly",
        k_bondcoat_nicocraly,
        848.15, 1373.15,
        "Zhu & Miller 1999, NASA/TM-1999-209069 Fig. 2 (literature values)",
        "NiCoCrAlY TBC bond coat; fit drawn over 575-1100 C only",
        NO_SOLIDUS,
        ""
    }},
    {"ss316", {
        "ss316",
        k_stainless_steel_316,
        260.0 + 273.15, 820.0 + 273.15,
        "Georgia Tech Gen3 CSP database, measured 260-820 C",
        "Austenitic stainless steel, corrosion resistant",
        NO_SOLIDUS,
        ""
    }},
    {"ss310", {
        "ss310",
        k_stainless_steel_310,
        260.0 + 273.15, 800.0 + 273.15,
        "Georgia Tech Gen3 CSP database, measured 260-800 C",
        "Austenitic 25Cr-20Ni heat-resisting stainless steel (cascade vanes, furnace parts)",
        (2470.0 - 32.0) / 1.8 + 273.15,
        "Rolled Alloys 310 data sheet, melting range 2470-2555 F"
    }},
    {"al6061", {
        "al6061",
        k_aluminum_6061,
        200.0, 600.0,
        "ASM Handbook, Aluminum Association",
        "Aluminum alloy, limited to <300degC due to T6 temper",
        NO_SOLIDUS,
        ""
    }},
    {"ysz", {
        "ysz",
        [](double T) { return k_tbc_ysz(T, 0.0); },  // Default: as-sprayed
        300.0, 1700.0,
        "Zhu & Miller 2000, NASA/TM-2000-210237 (plasma-sprayed, as sprayed)",
        "7-8 wt% Y2O3 stabilized zirconia, thermal barrier coating",
        NO_SOLIDUS,
        ""
    }}
};

const ThermalMaterial& get_material(const std::string& name) {
    auto it = material_database.find(name);
    if (it == material_database.end()) {
        throw std::runtime_error("Material '" + name + "' not found in database");
    }
    return it->second;
}

double material_conductivity(const std::string& name, double T, bool extrapolate) {
    const auto& mat = get_material(name);
    // NaN would fail both range tests below and come out as the clamped end
    // value (#769): refuse a non-finite T on every path.
    validation::require_finite<std::runtime_error>(T, "material_conductivity: T must be finite");
    if (T >= mat.T_min && T <= mat.T_max) {
        return mat.k_func(T);
    }
    if (!extrapolate) {
        return mat.k_func(T < mat.T_min ? mat.T_min : mat.T_max);
    }
    if (!(T > 0.0)) {
        throw std::runtime_error("material_conductivity: T must be > 0 K for " + name +
                                 ", got " + std::to_string(T));
    }
    // End slope from a one-sided difference inside the range: exact on a
    // table's end segment (every segment is longer than h), the tangent to
    // O(h) on a fitted polynomial.
    constexpr double h = 1e-2;
    const bool below = T < mat.T_min;
    const double T_end = below ? mat.T_min : mat.T_max;
    const double T_in = below ? mat.T_min + h : mat.T_max - h;
    const double k_end = mat.k_func(T_end);
    const double slope = (k_end - mat.k_func(T_in)) / (T_end - T_in);
    const double k = k_end + slope * (T - T_end);
    if (!(k > 0.0)) {
        throw std::runtime_error("material_conductivity: extrapolating " + name + " to " +
                                 std::to_string(T) + " K gives k = " + std::to_string(k) +
                                 " W/(m*K) <= 0");
    }
    return k;
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
