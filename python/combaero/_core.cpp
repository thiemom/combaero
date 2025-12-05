#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>

#include "thermo.h"
#include "transport.h"
#include "combustion.h"
#include "equilibrium.h"
#include "humidair.h"
#include "species_common_names.h"
#include "state.h"

namespace py = pybind11;

static std::vector<double> to_vec(
    py::array_t<double, py::array::c_style | py::array::forcecast> arr)
{
    auto buf = arr.unchecked<1>();
    std::vector<double> v(buf.shape(0));
    for (ssize_t i = 0; i < buf.shape(0); ++i)
        v[static_cast<std::size_t>(i)] = buf(i);
    return v;
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "Python bindings for combaero core";

    // Species metadata helpers
    m.def("num_species", &num_species, "Number of thermo species in the internal tables.");

    m.def(
        "species_name",
        &species_name,
        py::arg("index"),
        "Return the canonical species name for a given index."
    );

    m.def(
        "species_index_from_name",
        &species_index_from_name,
        py::arg("name"),
        "Return the species index for a given canonical name."
    );

    m.def(
        "species_molar_mass",
        &species_molar_mass,
        py::arg("index"),
        "Return the molar mass [g/mol] for a given species index."
    );

    m.def(
        "species_molar_mass_from_name",
        &species_molar_mass_from_name,
        py::arg("name"),
        "Return the molar mass [g/mol] for a given species name."
    );

    // Common names for species symbols
    m.def(
        "species_common_names",
        []() {
            return combaero::species_common_names;
        },
        "Mapping from canonical species symbols (e.g. 'CH4') to human-readable common names (e.g. 'Methane')."
    );

    m.def(
        "mixture_h",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return h(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture enthalpy at temperature T and mole fractions X."
    );

    // Thermodynamic mixture properties
    m.def(
        "cp",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return cp(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture isobaric heat capacity Cp(T, X) [J/(mol·K)]."
    );

    m.def(
        "h",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return h(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture enthalpy h(T, X) [J/mol]."
    );

    m.def(
        "s",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double P,
           double P_ref)
        {
            auto X = to_vec(X_arr);
            return s(T, X, P, P_ref);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P"),
        py::arg("P_ref") = 101325.0,
        "Mixture entropy s(T, X, P, P_ref) [J/(mol·K)]."
    );

    m.def(
        "cv",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return cv(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture isochoric heat capacity Cv(T, X) [J/(mol·K)]."
    );

    m.def(
        "density",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return density(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Mixture density rho(T, P, X) [kg/m^3]."
    );

    m.def(
        "specific_gas_constant",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return specific_gas_constant(X);
        },
        py::arg("X"),
        "Mixture specific gas constant R_s(X) [J/(kg·K)]."
    );

    m.def(
        "isentropic_expansion_coefficient",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return isentropic_expansion_coefficient(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Isentropic expansion coefficient gamma(T, X) = Cp/Cv."
    );

    m.def(
        "speed_of_sound",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return speed_of_sound(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Speed of sound c(T, X) [m/s]."
    );

    // Transport properties
    m.def(
        "viscosity",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return viscosity(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Dynamic viscosity mu(T, P, X) [Pa·s]."
    );

    m.def(
        "thermal_conductivity",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return thermal_conductivity(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Thermal conductivity lambda(T, P, X) [W/(m·K)]."
    );

    m.def(
        "kinematic_viscosity",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return kinematic_viscosity(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Kinematic viscosity nu(T, P, X) [m^2/s]."
    );

    m.def(
        "prandtl",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return prandtl(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Prandtl number Pr(T, P, X) [-]."
    );

    // Humid air utilities
    m.def(
        "standard_dry_air_composition",
        &standard_dry_air_composition,
        "Standard dry air mole-fraction composition over the thermo species set."
    );

    m.def(
        "humid_air_composition",
        &humid_air_composition,
        py::arg("T"),
        py::arg("P"),
        py::arg("RH"),
        "Humid air mole-fraction composition for temperature T [K], pressure P [Pa], and relative humidity RH (0-1)."
    );

    m.def(
        "dewpoint",
        &dewpoint,
        py::arg("T"),
        py::arg("P"),
        py::arg("RH"),
        "Dewpoint temperature [K] for ambient T [K], pressure P [Pa], and RH (0-1)."
    );

    m.def(
        "relative_humidity_from_dewpoint",
        &relative_humidity_from_dewpoint,
        py::arg("T"),
        py::arg("Tdp"),
        py::arg("P"),
        "Relative humidity (0-1) from dry-bulb T [K], dewpoint Tdp [K], and pressure P [Pa]."
    );

    // Composition conversion helpers
    m.def(
        "mole_to_mass",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return mole_to_mass(X);
        },
        py::arg("X"),
        "Convert mole fractions X to mass fractions Y.")
    ;

    m.def(
        "mass_to_mole",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_arr)
        {
            auto Y = to_vec(Y_arr);
            return mass_to_mole(Y);
        },
        py::arg("Y"),
        "Convert mass fractions Y to mole fractions X.")
    ;

    // Fraction utilities
    m.def(
        "normalize_fractions",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> arr)
        {
            auto v = to_vec(arr);
            return normalize_fractions(v);
        },
        py::arg("fractions"),
        "Normalize a vector of fractions so it sums to 1.0 (or return all zeros if input is all zeros)."
    );

    m.def(
        "convert_to_dry_fractions",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> arr)
        {
            auto v = to_vec(arr);
            return convert_to_dry_fractions(v);
        },
        py::arg("mole_fractions"),
        "Convert mole fractions to dry fractions by removing water vapor and renormalizing."
    );

    // Equivalence ratio (mole basis)
    m.def(
        "equivalence_ratio_mole",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_mix_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_ox_arr)
        {
            auto X_mix  = to_vec(X_mix_arr);
            auto X_fuel = to_vec(X_fuel_arr);
            auto X_ox   = to_vec(X_ox_arr);
            return equivalence_ratio_mole(X_mix, X_fuel, X_ox);
        },
        py::arg("X_mix"),
        py::arg("X_fuel"),
        py::arg("X_ox"),
        "Equivalence ratio phi (mole basis) for mixture X_mix formed from fuel X_fuel and oxidizer X_ox.")
    ;

    m.def(
        "set_equivalence_ratio_mole",
        [](double phi,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_ox_arr)
        {
            auto X_fuel = to_vec(X_fuel_arr);
            auto X_ox   = to_vec(X_ox_arr);
            return set_equivalence_ratio_mole(phi, X_fuel, X_ox);
        },
        py::arg("phi"),
        py::arg("X_fuel"),
        py::arg("X_ox"),
        "Construct mole-fraction mixture X_mix with target phi from fuel X_fuel and oxidizer X_ox.")
    ;

    // Equivalence ratio (mass basis)
    m.def(
        "equivalence_ratio_mass",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_mix_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_mix  = to_vec(Y_mix_arr);
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return equivalence_ratio_mass(Y_mix, Y_fuel, Y_ox);
        },
        py::arg("Y_mix"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Equivalence ratio phi (mass basis) for mixture Y_mix formed from fuel Y_fuel and oxidizer Y_ox.")
    ;

    m.def(
        "set_equivalence_ratio_mass",
        [](double phi,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return set_equivalence_ratio_mass(phi, Y_fuel, Y_ox);
        },
        py::arg("phi"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Construct mass-fraction mixture Y_mix with target phi from fuel Y_fuel and oxidizer Y_ox.")
    ;

    // Bilger-based helpers (mass basis)
    m.def(
        "bilger_stoich_mixture_fraction_mass",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return bilger_stoich_mixture_fraction_mass(Y_fuel, Y_ox);
        },
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Stoichiometric Bilger mixture fraction Z_st (mass basis) for fuel Y_fuel and oxidizer Y_ox.");

    m.def(
        "bilger_Z_from_equivalence_ratio_mass",
        [](double phi,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return bilger_Z_from_equivalence_ratio_mass(phi, Y_fuel, Y_ox);
        },
        py::arg("phi"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Convert mass-basis equivalence ratio phi to Bilger mixture fraction Z for fuel/oxidizer streams.");

    m.def(
        "equivalence_ratio_from_bilger_Z_mass",
        [](double Z,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return equivalence_ratio_from_bilger_Z_mass(Z, Y_fuel, Y_ox);
        },
        py::arg("Z"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Convert Bilger mixture fraction Z to mass-basis equivalence ratio phi for fuel/oxidizer streams.");

    // Inverse temperature calculations
    m.def(
        "calc_T_from_h",
        [](double h_target,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double T_guess,
           double tol,
           std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return calc_T_from_h(h_target, X, T_guess, tol, max_iter);
        },
        py::arg("h_target"),
        py::arg("X"),
        py::arg("T_guess") = 300.0,
        py::arg("tol") = 1.0e-6,
        py::arg("max_iter") = 50,
        "Solve for temperature T given target enthalpy h_target and composition X."
    );

    m.def(
        "calc_T_from_s",
        [](double s_target,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double T_guess,
           double tol,
           std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return calc_T_from_s(s_target, P, X, T_guess, tol, max_iter);
        },
        py::arg("s_target"),
        py::arg("P"),
        py::arg("X"),
        py::arg("T_guess") = 300.0,
        py::arg("tol") = 1.0e-6,
        py::arg("max_iter") = 50,
        "Solve for temperature T given target entropy s_target, pressure P, and composition X."
    );

    m.def(
        "calc_T_from_cp",
        [](double cp_target,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double T_guess,
           double tol,
           std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return calc_T_from_cp(cp_target, X, T_guess, tol, max_iter);
        },
        py::arg("cp_target"),
        py::arg("X"),
        py::arg("T_guess") = 300.0,
        py::arg("tol") = 1.0e-6,
        py::arg("max_iter") = 50,
        "Solve for temperature T given target Cp cp_target and composition X."
    );

    // Oxygen requirement helpers (scalar fuel index and mixtures)
    m.def(
        "oxygen_required_per_mol_fuel",
        &oxygen_required_per_mol_fuel,
        py::arg("fuel_index"),
        "Moles O2 required per mole of pure fuel species (by index)."
    );

    m.def(
        "oxygen_required_per_kg_fuel",
        &oxygen_required_per_kg_fuel,
        py::arg("fuel_index"),
        "Mass of O2 [kg] required per kg of pure fuel species (by index)."
    );

    m.def(
        "oxygen_required_per_mol_mixture",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return oxygen_required_per_mol_mixture(X);
        },
        py::arg("X"),
        "Moles O2 required per mole of fuel mixture with mole fractions X."
    );

    m.def(
        "oxygen_required_per_kg_mixture",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return oxygen_required_per_kg_mixture(X);
        },
        py::arg("X"),
        "Mass of O2 [kg] required per kg of fuel mixture with mole fractions X."
    );

    // Complete combustion helper
    m.def(
        "complete_combustion_to_CO2_H2O",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return complete_combustion_to_CO2_H2O(X);
        },
        py::arg("X"),
        "Complete combustion of mixture X to CO2 and H2O (returns product mole fractions)."
    );

    m.def(
        "complete_combustion_to_CO2_H2O_with_fraction",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            double f = 0.0;
            auto X_out = complete_combustion_to_CO2_H2O(X, f);
            return py::make_tuple(X_out, f);
        },
        py::arg("X"),
        "Complete combustion of mixture X to CO2 and H2O, returning (product mole fractions, fuel burn fraction f)."
    );

    m.def(
        "adiabatic_T_wgs",
        [](double T_in,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_in_arr,
           double P)
        {
            auto X_in = to_vec(X_in_arr);
            if (X_in.empty())
                throw std::runtime_error("X_in must be non-empty");

            State in;
            in.T = T_in;
            in.P = P;
            in.X = X_in;

            State out = wgs_equilibrium_adiabatic(in);
            return out.T;
        },
        py::arg("T_in"),
        py::arg("X_in"),
        py::arg("P") = 101325.0,
        "Adiabatic flame temperature with WGS limited equilibrium.\n\n"
        "T_in : inlet temperature [K]\n"
        "X_in : inlet mole fractions (1D NumPy array)\n"
        "P    : pressure [Pa] (default 101325)"
    );

    // State struct binding
    py::class_<State>(m, "State")
        .def(py::init<>())
        .def_readwrite("T", &State::T, "Temperature [K]")
        .def_readwrite("P", &State::P, "Pressure [Pa]")
        .def_readwrite("X", &State::X, "Mole fractions [-]")
        // Property getters
        .def("mw", &State::mw, "Molecular weight [g/mol]")
        .def("cp", &State::cp, "Specific heat at constant pressure [J/(mol·K)]")
        .def("h", &State::h, "Specific enthalpy [J/mol]")
        .def("s", &State::s, "Specific entropy [J/(mol·K)]")
        .def("cv", &State::cv, "Specific heat at constant volume [J/(mol·K)]")
        .def("u", &State::u, "Specific internal energy [J/mol]")
        .def("rho", &State::rho, "Density [kg/m³]")
        .def("R", &State::R, "Specific gas constant [J/(mol·K)]")
        .def("gamma", &State::gamma, "Isentropic expansion coefficient [-]")
        .def("a", &State::a, "Speed of sound [m/s]")
        // Setters
        .def("set_T", &State::set_T, py::arg("T"), "Set temperature [K]")
        .def("set_P", &State::set_P, py::arg("P"), "Set pressure [Pa]")
        .def("set_X", [](State& s, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) -> State& {
            s.X = to_vec(X_arr);
            return s;
        }, py::arg("X"), "Set mole fractions [-]");

    // Stream struct binding
    py::class_<Stream>(m, "Stream")
        .def(py::init<>())
        .def_readwrite("state", &Stream::state, "Thermodynamic state")
        .def_readwrite("mdot", &Stream::mdot, "Mass flow rate [kg/s]")
        // Convenience accessors
        .def("T", &Stream::T, "Temperature [K]")
        .def("P", &Stream::P, "Pressure [Pa]")
        .def("X", &Stream::X, "Mole fractions [-]")
        // Property getters
        .def("mw", &Stream::mw, "Molecular weight [g/mol]")
        .def("cp", &Stream::cp, "Specific heat at constant pressure [J/(mol·K)]")
        .def("h", &Stream::h, "Specific enthalpy [J/mol]")
        .def("s", &Stream::s, "Specific entropy [J/(mol·K)]")
        .def("rho", &Stream::rho, "Density [kg/m³]")
        // Setters
        .def("set_T", &Stream::set_T, py::arg("T"), "Set temperature [K]")
        .def("set_P", &Stream::set_P, py::arg("P"), "Set pressure [Pa]")
        .def("set_X", [](Stream& s, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) -> Stream& {
            s.state.X = to_vec(X_arr);
            return s;
        }, py::arg("X"), "Set mole fractions [-]")
        .def("set_mdot", &Stream::set_mdot, py::arg("mdot"), "Set mass flow rate [kg/s]");

    // Stream mixing function
    m.def(
        "mix",
        [](const std::vector<Stream>& streams, double P_out) {
            return mix(streams, P_out);
        },
        py::arg("streams"),
        py::arg("P_out") = -1.0,
        "Mix multiple streams, solving mass and enthalpy balance.\n\n"
        "streams : list of Stream objects\n"
        "P_out   : output pressure [Pa] (default: minimum inlet pressure)\n\n"
        "Returns mixed Stream."
    );

    // State-based combustion functions
    m.def(
        "complete_combustion",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return complete_combustion(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Adiabatic complete combustion to CO2 and H2O.\n\n"
        "Returns State with adiabatic flame temperature and burned composition."
    );

    m.def(
        "complete_combustion_isothermal",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return complete_combustion_isothermal(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Isothermal complete combustion to CO2 and H2O.\n\n"
        "Returns State with same temperature and burned composition."
    );

    // State-based WGS equilibrium functions
    m.def(
        "wgs_equilibrium",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return wgs_equilibrium(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Isothermal WGS equilibrium (CO + H2O <-> CO2 + H2).\n\n"
        "Returns State with equilibrium composition at input temperature."
    );

    m.def(
        "wgs_equilibrium_adiabatic",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return wgs_equilibrium_adiabatic(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Adiabatic WGS equilibrium (CO + H2O <-> CO2 + H2).\n\n"
        "Returns State with equilibrium temperature and composition."
    );
}