#pragma once

#include <string>
#include <unordered_map>
#include <stdexcept>

namespace combaero {

// Mapping from canonical species symbols (as used in thermo_transport_data.h)
// to human-readable common names.
//
// This header is kept separate from the generated thermo_transport_data.h so
// that the human-readable naming can evolve independently of the data tables.
inline const std::unordered_map<std::string, std::string> formula_to_name{
    {"N2",     "Nitrogen"},
    {"O2",     "Oxygen"},
    {"AR",     "Argon"},
    {"CO2",    "Carbon dioxide"},
    {"H2O",    "Water"},
    {"CH4",    "Methane"},
    {"C2H6",   "Ethane"},
    {"C3H8",   "Propane"},
    {"IC4H10", "Isobutane"},
    {"NC5H12", "n-Pentane"},
    {"NC6H14", "n-Hexane"},
    {"NC7H16", "n-Heptane"},
    {"CO",     "Carbon monoxide"},
    {"H2",     "Hydrogen"},
    {"NH3",    "Ammonia"},
};

// Inverse mapping: common name -> formula
inline const std::unordered_map<std::string, std::string> name_to_formula{
    {"Nitrogen",        "N2"},
    {"Oxygen",          "O2"},
    {"Argon",           "AR"},
    {"Carbon dioxide",  "CO2"},
    {"Water",           "H2O"},
    {"Methane",         "CH4"},
    {"Ethane",          "C2H6"},
    {"Propane",         "C3H8"},
    {"Isobutane",       "IC4H10"},
    {"n-Pentane",       "NC5H12"},
    {"n-Hexane",        "NC6H14"},
    {"n-Heptane",       "NC7H16"},
    {"Carbon monoxide", "CO"},
    {"Hydrogen",        "H2"},
    {"Ammonia",         "NH3"},
};

// Lookup functions with error handling

/// Get common name from formula (e.g., "CH4" -> "Methane")
inline std::string common_name(const std::string& formula) {
    auto it = formula_to_name.find(formula);
    if (it == formula_to_name.end()) {
        throw std::out_of_range("common_name: unknown formula '" + formula + "'");
    }
    return it->second;
}

/// Get formula from common name (e.g., "Methane" -> "CH4")
inline std::string formula(const std::string& name) {
    auto it = name_to_formula.find(name);
    if (it == name_to_formula.end()) {
        throw std::out_of_range("formula: unknown name '" + name + "'");
    }
    return it->second;
}

} // namespace combaero
