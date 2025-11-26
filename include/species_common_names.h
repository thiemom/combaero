#pragma once

#include <string>
#include <unordered_map>

namespace combaero {

// Mapping from canonical species symbols (as used in thermo_transport_data.h)
// to human-readable common names.
//
// This header is kept separate from the generated thermo_transport_data.h so
// that the human-readable naming can evolve independently of the data tables.
inline const std::unordered_map<std::string, std::string> species_common_names{
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
};

} // namespace combaero
