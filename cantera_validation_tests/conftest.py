import os
import secrets
import zlib

import numpy as np
import pytest

# Random states are drawn fresh each run: that is how test_property_SH caught
# Cantera's own argon weight (#680), on a draw no fixed seed had made. The
# seed is printed in the header and on every failure; --seed N (or
# COMBAERO_TEST_SEED=N) replays the run, and each test's draw depends only on
# the seed and its own id, so one failing test reruns on its own.
SEED_ENV = "COMBAERO_TEST_SEED"


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--seed",
        type=int,
        default=None,
        help=f"seed for the random states (default: fresh, or ${SEED_ENV})",
    )


def pytest_configure(config: pytest.Config) -> None:
    seed = config.getoption("--seed")
    if seed is None:
        env = os.environ.get(SEED_ENV)
        seed = int(env) if env else secrets.randbelow(2**32)
    config._combaero_seed = seed
    # xdist workers start after this and inherit the environment: one seed
    # for the whole run, the one the header shows.
    os.environ[SEED_ENV] = str(seed)


def pytest_report_header(config: pytest.Config) -> str:
    return f"random seed: {config._combaero_seed} (replay: --seed {config._combaero_seed})"


@pytest.fixture(autouse=True)
def _seeded_random(request: pytest.FixtureRequest) -> None:
    seed = request.config._combaero_seed
    np.random.seed((seed ^ zlib.crc32(request.node.nodeid.encode())) % 2**32)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item: pytest.Item, call: pytest.CallInfo):
    outcome = yield
    rep = outcome.get_result()
    if rep.failed:
        seed = item.config._combaero_seed
        rep.sections.append(("random seed", f"replay {item.nodeid} with --seed {seed}"))


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
        "NH3": "NH3",
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


@pytest.fixture(scope="session")
def nasa9_gas(cantera):
    """Cantera gas phase for NASA9 polynomial validation."""
    return cantera.Solution("gri30.yaml")


@pytest.fixture
def tolerance_config():
    """Standard tolerance configuration for validation tests.

    Tolerances based on measured deviations:
    - Temperature: 5.0 K (observed max: 4.6 K for C3H8, NASA-7 vs NASA-9)
    - Transport: 7% (observed after #485: viscosity <= 3.4%, conductivity <= ~6%
      at 1000-1500 K and with polar species -- Eucken vs Cantera's full model)
    - Equilibrium: 0.01% composition, 1 K temperature (observed: 0.002%, 0.0 K)
    - Enthalpy: 1.5% (observed: up to 1.02%, small data source differences)

    NOTE: this was 35%, which let O2's viscosity sit 44.5% low (and air's 11.8%)
    for months: the generated transport data took NUIGMech1.1 parameters its
    own authors flag as unreferenced estimates (#485). The tolerance is set
    just above the largest difference the models themselves explain.
    """
    return {
        "temperature": 5.0,  # K
        "mole_fraction": 0.01,  # absolute
        "enthalpy": 0.015,  # relative (1.5%) - small differences in thermo data
        "transport": 0.07,  # relative (7%) - conductivity model differences, see above
        "density": 0.01,  # relative (1%)
        "equilibrium_composition": 0.0001,  # absolute (0.01%) - WGS equilibrium
        "equilibrium_temperature": 1.0,  # K - adiabatic equilibrium
        "equilibrium_constant": 0.002,  # relative (0.2%) - Kp validation
    }
