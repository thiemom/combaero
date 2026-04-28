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
    "_solver_tools",
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
    "ChamberResult",
    "MixtureState",
    "Stream",
    "AreaChangeResult",
    "AreaChangeElementResult",
    "State::DP",
    "State::HP",
    "State::PV",
    "State::SH",
    "State::SP",
    "State::SV",
    "State::UP",
    "State::UV",
    "State::VH",
    "State::set_P",
    "State::set_T",
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
                if method_name in IGNORE_LIST:
                    continue
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
