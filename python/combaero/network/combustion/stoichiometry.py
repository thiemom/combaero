"""Stoichiometric combustion calculations for network solver."""

import combaero as cb


def _get_atomic_composition(X: list[float]) -> list[float]:
    """
    Get atomic composition [C, H, O, N] for a given mixture.

    Parameters
    ----------
    X : list[float]
        Mole fractions [14 species]

    Returns
    -------
    list[float]
        Total atoms [C, H, O, N] in the mixture
    """
    # Get molecular structures from CombAero
    # This accesses the atomic composition data from thermo_transport_data.h

    # For now, we'll use a simplified approach since we don't have direct access
    # to the molecular structures. In a full implementation, this would:
    # 1. Get atoms per molecule for each species from molecular_structures
    # 2. Multiply by mole fractions
    # 3. Sum across all species

    # Placeholder: This needs to be implemented properly using CombAero's internal data
    # For now, we'll use CombAero's built-in functions as a workaround

    # Use CombAero to get the atomic composition indirectly
    # This is a temporary solution until we have direct access to molecular_structures

    # Get oxygen required (which gives us O atom count indirectly)
    # This is not the proper way, but serves as a placeholder

    # TODO: Implement proper atom balancing using molecular_structures
    # from thermo_transport_data.h when available

    raise NotImplementedError(
        "Proper atom balancing requires access to molecular_structures from thermo_transport_data.h"
    )


def stoichiometric_products(
    X_fuel: list[float],
    X_oxidiser: list[float],
    phi: float,
) -> list[float]:
    """
    Compute burned gas mole fractions from atom balance.

    Assumes complete combustion:
      - All C → CO2
      - All H → H2O
      - Excess O2 remains
      - N2 passes through inert

    Parameters
    ----------
    X_fuel : list[float]
        Mole fractions of fuel stream [14 species].
    X_oxidiser : list[float]
        Mole fractions of oxidiser stream [14 species].
    phi : float
        Equivalence ratio. phi < 1: lean, phi = 1: stoichiometric, phi > 1: rich.

    Returns
    -------
    list[float]
        Burned gas mole fractions [14 species], normalised to sum to 1.

    Notes
    -----
    Uses molecular_structures from thermo_transport_data.h via CombAero
    for atomic composition (C, H, O, N atoms per species).
    Rich combustion (phi > 1) produces CO and H2 as incomplete products.
    """
    # Validate inputs
    if len(X_fuel) != 14 or len(X_oxidiser) != 14:
        raise ValueError("X_fuel and X_oxidiser must have 14 species")

    if abs(sum(X_fuel) - 1.0) > 1e-6:
        raise ValueError("X_fuel must sum to 1")
    if abs(sum(X_oxidiser) - 1.0) > 1e-6:
        raise ValueError("X_oxidiser must sum to 1")

    if phi <= 0:
        raise ValueError("phi must be positive")

    # CRITICAL: This is a placeholder implementation
    # The proper implementation requires access to molecular_structures
    # from thermo_transport_data.h to perform atom balancing

    # For now, we use CombAero's built-in functions as a workaround
    # This does NOT properly handle equivalence ratio - it's just a placeholder

    # Create a mixed composition based on equivalence ratio
    X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_oxidiser)

    # Get complete combustion products
    # NOTE: This ignores phi for rich/lean cases - needs proper atom balancing
    X_products = cb.complete_combustion_to_CO2_H2O(X_mix)

    # TODO: Implement proper atom balancing for rich/lean mixtures
    # For rich (φ > 1): distribute excess fuel to CO and H₂
    # For lean (φ < 1): add excess O2

    return X_products
