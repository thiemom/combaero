from __future__ import annotations

import numpy as np

import combaero as ca
from combaero.species import SpeciesLocator


def main() -> None:
    species_locator = SpeciesLocator.from_core()

    print("Species layout:")
    for i, name in enumerate(species_locator.names):
        mm = species_locator.molar_masses[i]
        print(f"  {i:2d}: {name:>6s}  M = {mm:.3f} g/mol")

    # Example: pure CH4 mixture on the combaero species grid
    X = species_locator.empty()
    idx_ch4 = species_locator.indices["CH4"]
    X[idx_ch4] = 1.0

    print("\nMixture enthalpy for pure CH4:")
    for T in (300.0, 600.0, 1200.0):
        h_val = ca.mixture_h(T, X)
        print(f"  T = {T:7.1f} K  h = {h_val: .6e} J/mol")


if __name__ == "__main__":
    main()
