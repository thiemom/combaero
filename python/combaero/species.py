from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping

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
    indices: Dict[str, int]
    molar_masses: np.ndarray

    @classmethod
    def from_core(cls) -> "SpeciesLocator":
        """Build a layout from the compiled core metadata."""

        n = _core.num_species()
        names = [_core.species_name(i) for i in range(n)]
        indices: Dict[str, int] = {name: i for i, name in enumerate(names)}
        molar_masses = np.array(
            [_core.species_molar_mass(i) for i in range(n)], dtype=float
        )
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
            X[idx] = float(value)
        return X
