import inspect

import combaero as cb

# These are pure Python utilities, constants, or modules that
# do not have physical units.
IGNORE_LIST = {
    # Modules
    "compressible",
    "incompressible",
    "heat_transfer",
    "network",
    "species",
    "geometry",
    "orifice_mdot_and_jacobian",
    "pressure_loss_and_jacobian",
    "lossless_pressure_and_jacobian",
    "nusselt_and_jacobian_dittus_boelter",
    "friction_and_jacobian_haaland",
    "friction_and_jacobian_serghides",
    "friction_and_jacobian_colebrook",
    "friction_and_jacobian_petukhov",
    "friction_and_jacobian",
    "density_and_jacobians",
    "enthalpy_and_jacobian",
    "viscosity_and_jacobians",
    "nusselt_and_jacobian_gnielinski",
    "nusselt_and_jacobian_sieder_tate",
    "nusselt_and_jacobian_petukhov",
    "pin_fin_nusselt_and_jacobian",
    "pin_fin_friction_and_jacobian",
    "dimple_nusselt_enhancement_and_jacobian",
    "dimple_friction_multiplier_and_jacobian",
    "rib_enhancement_factor_high_re_and_jacobian",
    "impingement_nusselt_and_jacobian",
    "film_cooling_effectiveness_and_jacobian",
    "effusion_effectiveness_and_jacobian",
    "mach_number_and_jacobian_v",
    "T_adiabatic_wall_and_jacobian_v",
    "T0_from_static_and_jacobian_M",
    "P0_from_static_and_jacobian_M",
    # Utilities
    "get_warning_handler",
    "set_warning_handler",
    "is_well_behaved",
    "suppress_warnings",
    # Unit query API itself
    "has_units",
    "get_units",
    "output_units",
    "input_units",
    "all_units",
    "list_functions_with_units",
    # String manipulation / helpers
    "common_name",
    "formula_to_name",
    "name_to_formula",
    "formula",
    "list_materials",
    "num_species",
    "species_name",
    "species_index_from_name",
    "species_molar_mass",
    "species_molar_mass_from_name",
    # Mix is a special case function combining streams
    "mix",
    # Enums
    "CombustionMethod",
    "BoundaryCondition",
    "CorrelationValidity",
    "CorrelationResult",
}

# Special magic methods created by pybind11 that aren't public calculation points
IGNORE_METHODS = {
    "__init__",
    "__repr__",
    "__str__",
    "__copy__",
    "__deepcopy__",
    "__doc__",
    "__module__",
    "__class__",
}


def test_api_unit_sync():
    """
    Dynamically discover all exposed functions and classes from combaero,
    and assert that they have corresponding unit metadata defined in units_data.h.
    This prevents the documentation and API from drifting out of sync.
    """
    missing = []

    for name in cb.__all__:
        if name in IGNORE_LIST:
            continue

        obj = getattr(cb, name)

        # Skip strings like __version__
        if isinstance(obj, str):
            continue

        if inspect.isclass(obj):
            # For classes, check their public members
            for attr in dir(obj):
                if attr.startswith("_"):
                    continue
                if attr in IGNORE_METHODS:
                    continue

                method_name = f"{name}::{attr}"
                if not cb.has_units(method_name):
                    missing.append(method_name)
        else:
            # Standalone function or property exported directly
            if not cb.has_units(name):
                missing.append(name)

    # Use standard assert for pytest
    assert len(missing) == 0, (
        f"The following {len(missing)} API members are missing from units_data.h:\n"
        + "\n".join(missing)
    )
