from collections.abc import Mapping
from dataclasses import dataclass
from typing import Self

import numpy as np

from . import _core


@dataclass(frozen=True)
class SpeciesLocator:
    """Helper describing the thermo species layout used by combaero.

    This wraps the C++ thermo tables and provides convenient helpers to
    construct composition vectors for use with functions like
    :func:`combaero.mixture_h`.
    """

    names: list[str]
    indices: dict[str, int]
    molar_masses: np.ndarray

    @classmethod
    def from_core(cls) -> Self:
        """Build a layout from the compiled core metadata."""

        n = _core.num_species()
        names = [_core.species_name(i) for i in range(n)]
        indices: dict[str, int] = {name: i for i, name in enumerate(names)}
        molar_masses = np.array([_core.species_molar_mass(i) for i in range(n)], dtype=float)
        return cls(names=names, indices=indices, molar_masses=molar_masses)

    def empty(self) -> np.ndarray:
        """Return an all-zero composition vector of the correct length."""

        return np.zeros(len(self.names), dtype=float)

    def from_mapping(self, composition: Mapping[str, float]) -> np.ndarray:
        """Create a composition vector from a mapping of species name to value.

        Unknown names are resolved via the core's species_index_from_name,
        so aliases supported by the C++ side will also work here.
        """

        X = self.empty()
        for name, value in composition.items():
            idx = self.indices.get(name)
            if idx is None:
                idx = _core.species_index_from_name(name)
            if idx < 0 or idx >= len(X):
                raise ValueError(f"Invalid species index {idx} for name '{name}'")
            X[idx] = float(value)
        return X

    def pure_species(self, name: str) -> np.ndarray:
        """Return a composition vector of a pure species."""
        X = self.empty()
        idx = self.indices.get(name)
        if idx is None:
            idx = _core.species_index_from_name(name)
        if idx < 0 or idx >= len(X):
            raise ValueError(f"Invalid species index {idx} for name '{name}'")
        X[idx] = 1.0
        return X

    def to_mass(self, X: np.ndarray) -> np.ndarray:
        """Convert mole fractions to mass fractions."""
        return _core.mole_to_mass(X)

    def to_mole(self, Y: np.ndarray) -> np.ndarray:
        """Convert mass fractions to mole fractions."""
        return _core.mass_to_mole(Y)

    def dry_air(self) -> np.ndarray:
        """Return the standard dry air composition (mole fractions)."""
        return _core.dry_air()

    def dry_air_mass(self) -> np.ndarray:
        """Return the standard dry air composition (mass fractions)."""
        return self.to_mass(self.dry_air())

    def humid_air(self, T: float, P: float, RH: float) -> np.ndarray:
        """Return the humid air composition for T [K], P [Pa], and RH [0-1] (mole fractions)."""
        return _core.humid_air_composition(T, P, RH)

    def humid_air_mass(self, T: float, P: float, RH: float) -> np.ndarray:
        """Return the humid air composition for T [K], P [Pa], and RH [0-1] (mass fractions)."""
        return self.to_mass(self.humid_air(T, P, RH))


# Create a default instance for the module-level API
_default_locator = SpeciesLocator.from_core()


# Module-level access to species metadata
indices = _default_locator.indices
names = _default_locator.names
num_species = len(names)


def species_index(name: str) -> int:
    """Return the index of a species by name or alias."""
    idx = indices.get(name)
    if idx is None:
        idx = _core.species_index_from_name(name)
    return idx


def empty() -> np.ndarray:
    """Return an all-zero composition vector."""
    return _default_locator.empty()


def from_mapping(composition: Mapping[str, float]) -> np.ndarray:
    """Create a composition vector from a mapping of species name to value."""
    return _default_locator.from_mapping(composition)


def to_mass(X: np.ndarray) -> np.ndarray:
    """Convert mole fractions to mass fractions."""
    return _default_locator.to_mass(X)


def to_mole(Y: np.ndarray) -> np.ndarray:
    """Convert mass fractions to mole fractions."""
    return _default_locator.to_mole(Y)


def pure_species(name: str) -> np.ndarray:
    """Return a composition vector of a pure species."""
    return _default_locator.pure_species(name)


def dry_air() -> np.ndarray:
    """Return the standard dry air composition (mole fractions)."""
    return _default_locator.dry_air()


def dry_air_mass() -> np.ndarray:
    """Return the standard dry air composition (mass fractions)."""
    return _default_locator.dry_air_mass()


def humid_air(T: float, P: float, RH: float) -> np.ndarray:
    """Return the humid air composition (mole fractions)."""
    return _default_locator.humid_air(T, P, RH)


def _parse_comp_str(s: str) -> np.ndarray:
    """Parse a Cantera-style composition string into a mole-fraction array.

    Accepts 'CH4:1, O2:2, N2:7.52'. Values need not sum to 1; from_mapping
    stores them unnormalized as the C++ side expects raw mole fractions.
    """
    mapping: dict[str, float] = {}
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        name, _, val = part.partition(":")
        mapping[name.strip()] = float(val.strip())
    return from_mapping(mapping)


def humid_air_mass(T: float, P: float, RH: float) -> np.ndarray:
    """Return the humid air composition (mass fractions)."""
    return _default_locator.humid_air_mass(T, P, RH)
