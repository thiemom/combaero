#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>

#include "thermo_transport.h"
#include "equilibrium.h"

namespace py = pybind11;

// Build WgsConfig using species indices looked up by name
static combaero::equilibrium::WgsConfig make_wgs_cfg()
{
    combaero::equilibrium::WgsConfig cfg{};
    cfg.i_CO  = static_cast<std::size_t>(species_index_from_name("CO"));
    cfg.i_H2O = static_cast<std::size_t>(species_index_from_name("H2O"));
    cfg.i_CO2 = static_cast<std::size_t>(species_index_from_name("CO2"));
    cfg.i_H2  = static_cast<std::size_t>(species_index_from_name("H2"));
    return cfg;
}

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

    m.def(
        "adiabatic_T_wgs",
        [](double T_in,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_in_arr,
           py::object T_guess_obj)
        {
            auto X_in = to_vec(X_in_arr);
            if (X_in.empty())
                throw std::runtime_error("X_in must be non-empty");

            double T_guess = T_in + 500.0;
            if (!T_guess_obj.is_none())
                T_guess = T_guess_obj.cast<double>();

            // 1 mol basis
            std::vector<double> n_in = X_in;

            double H_target = h(T_in, X_in);
            auto cfg = make_wgs_cfg();

            double T_eq = combaero::equilibrium::solve_adiabatic_T_wgs(
                n_in, H_target, T_guess, cfg
            );
            return T_eq;
        },
        py::arg("T_in"),
        py::arg("X_in"),
        py::arg("T_guess") = py::none(),
        "Adiabatic flame temperature with WGS limited equilibrium.\n\n"
        "T_in   : inlet temperature [K]\n"
        "X_in   : inlet mole fractions (1D NumPy array)\n"
        "T_guess: optional initial guess for T [K]"
    );
}