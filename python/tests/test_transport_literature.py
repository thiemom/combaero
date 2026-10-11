"""Viscosity against literature, independent of any mechanism's database (#485).

The generated transport data once took NUIGMech1.1 Lennard-Jones parameters
that its own authors flag as unreferenced estimates ('theoret trans'): O2's
viscosity was 44.5% low, air's 11.8%, and the Cantera comparison (35%
tolerance) let it through. Pinned here to the reference values at 300 K,
low pressure.

References (dilute-gas viscosity, 300 K):
  N2, O2, Ar, air -- Lemmon & Jacobsen (2004), Int. J. Thermophys. 25, 21
  H2, CO2          -- NIST Chemistry WebBook fluid properties
"""

from __future__ import annotations

import pytest

import combaero as cb

REF_300K_UPAS = {"N2": 17.89, "O2": 20.65, "AR": 22.73, "H2": 8.96, "CO2": 15.02}


def _pure(name: str) -> list[float]:
    X = [0.0] * cb.num_species()
    X[cb.species_index_from_name(name)] = 1.0
    return X


@pytest.mark.parametrize("species", sorted(REF_300K_UPAS))
def test_pure_species_viscosity_matches_literature(species: str) -> None:
    """Chapman-Enskog with Lennard-Jones parameters is good to ~2% here."""
    mu = cb.viscosity(300.0, 1.0e5, _pure(species)) * 1e6
    assert mu == pytest.approx(REF_300K_UPAS[species], rel=0.025)


def test_air_viscosity_matches_literature() -> None:
    """18.54 uPa.s at 300 K. It read 16.35 with O2's flagged parameters."""
    mu = cb.viscosity(300.0, 1.0e5, cb.species.dry_air()) * 1e6
    assert mu == pytest.approx(18.54, rel=0.015)
