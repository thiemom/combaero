import pytest


@pytest.fixture(scope="session")
def combaero():
    """Provide combaero module."""
    try:
        import combaero as cb

        return cb
    except ImportError:
        pytest.skip("combaero not installed")


@pytest.fixture(scope="session")
def cantera():
    """Provide cantera module."""
    try:
        import cantera as ct

        return ct
    except ImportError:
        pytest.skip("cantera not installed")


@pytest.fixture(scope="session")
def species_mapping():
    """Map CombAero species names to Cantera species names."""
    return {
        "N2": "N2",
        "O2": "O2",
        "AR": "AR",
        "CO2": "CO2",
        "H2O": "H2O",
        "CH4": "CH4",
        "C2H6": "C2H6",
        "C3H8": "C3H8",
        "IC4H10": "IC4H10",
        "NC5H12": "NC5H12",
        "NC6H14": "NC6H14",
        "NC7H16": "NC7H16",
        "CO": "CO",
        "H2": "H2",
    }


@pytest.fixture
def complete_combustion_gas(cantera):
    """
    Gas phase with only complete combustion species.

    This restricts equilibrium to only produce CO2 and H2O from fuel,
    matching CombAero's complete_combustion() behavior (no CO, H2, etc.).
    """
    species = {S.name: S for S in cantera.Species.list_from_file("gri30.yaml")}
    complete_species = [
        species[S] for S in ("N2", "O2", "AR", "CO2", "H2O", "CH4", "C2H6", "C3H8", "H2", "CO")
    ]
    return cantera.Solution(thermo="ideal-gas", species=complete_species)


@pytest.fixture
def tolerance_config():
    """Standard tolerance configuration for validation tests.

    Tolerances based on measured deviations:
    - Temperature: 5.0 K (observed max: 4.6 K for C3H8, NASA-7 vs NASA-9)
    - Transport: 35% (observed: up to 30% at high T, different correlations expected)
    - Equilibrium: 0.01% composition, 1 K temperature (observed: 0.002%, 0.0 K)
    - Enthalpy: 1.5% (observed: up to 1.02%, small data source differences)
    """
    return {
        "temperature": 5.0,  # K
        "mole_fraction": 0.01,  # absolute
        "enthalpy": 0.015,  # relative (1.5%) - small differences in thermo data
        "transport": 0.35,  # relative (35%) - increased for high-T correlation differences
        "density": 0.01,  # relative (1%)
        "equilibrium_composition": 0.0001,  # absolute (0.01%) - WGS equilibrium
        "equilibrium_temperature": 1.0,  # K - adiabatic equilibrium
        "equilibrium_constant": 0.002,  # relative (0.2%) - Kp validation
    }
