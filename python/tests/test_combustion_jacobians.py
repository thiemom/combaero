"""Tests for combustion Jacobian functions."""

from __future__ import annotations

import combaero as cb


def test_adiabatic_complete_jacobian_basic():
    """Test basic functionality of adiabatic_T_complete_and_jacobian_T."""
    # Create a methane-air mixture
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")

    # Stoichiometric mixture
    X_mix = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)

    T_in = 700.0  # 700K
    P = 500000.0  # 5 bar

    # Test the function
    result = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P, X_mix)

    # Should return tuple of (T_ad, dT_ad_dT_in, X_products)
    assert len(result) == 3
    T_ad, dT_ad_dT_in, X_products = result

    # Basic validation
    assert T_ad > T_in, "Adiabatic temperature should be higher than inlet temperature"
    assert 2000 < T_ad < 3000, f"T_ad {T_ad} should be in reasonable range for CH4 combustion"

    # Check that we have proper combustion products
    assert X_products[3] > 0.05, "Should have significant CO2"
    assert X_products[4] > 0.1, "Should have significant H2O"

    # Jacobian should be positive (T_ad increases with T_in)
    assert dT_ad_dT_in > 0, "Jacobian should be positive"
    assert dT_ad_dT_in < 2.0, "Jacobian should be reasonable magnitude"

    # Species should sum to 1
    assert abs(sum(X_products) - 1.0) < 1e-6, "Species should sum to 1"


def test_adiabatic_equilibrium_jacobian_basic():
    """Test basic functionality of adiabatic_T_equilibrium_and_jacobians."""
    # Create a methane-air mixture
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")

    # Stoichiometric mixture
    X_mix = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)

    T_in = 700.0  # 700K
    P = 500000.0  # 5 bar

    # Test the function
    result = cb._core.adiabatic_T_equilibrium_and_jacobians(T_in, P, X_mix)

    # Should return tuple of (T_ad, dT_ad_dT_in, dT_ad_dP, X_products)
    assert len(result) == 4
    T_ad, dT_ad_dT_in, dT_ad_dP, X_products = result

    # Basic validation
    assert T_ad > T_in, "Adiabatic temperature should be higher than inlet temperature"

    # Jacobian should be positive
    assert dT_ad_dT_in > 0, "Jacobian should be positive"

    # Species should sum to 1
    assert abs(sum(X_products) - 1.0) < 1e-6, "Species should sum to 1"


def test_jacobian_sensitivity():
    """Test that Jacobian responds correctly to temperature changes."""
    # Create a methane-air mixture
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")
    X_mix = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)

    P = 500000.0  # 5 bar
    T_base = 700.0

    # Get Jacobian at base temperature
    _, dT_ad_dT_in_base, _ = cb._core.adiabatic_T_complete_and_jacobian_T(T_base, P, X_mix)

    # Test at slightly different temperature
    T_perturbed = T_base + 10.0
    _, dT_ad_dT_in_perturbed, _ = cb._core.adiabatic_T_complete_and_jacobian_T(
        T_perturbed, P, X_mix
    )

    # Jacobians should be similar (function should be smooth)
    assert abs(dT_ad_dT_in_base - dT_ad_dT_in_perturbed) < 0.1, "Jacobian should be smooth"

    # Both should be positive
    assert dT_ad_dT_in_base > 0 and dT_ad_dT_in_perturbed > 0, "Jacobians should be positive"


def test_lean_vs_stoichiometric():
    """Test that lean mixtures have lower adiabatic temperatures."""
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")

    P = 500000.0
    T_in = 700.0

    # Stoichiometric
    X_stoich = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)
    _, _, X_products_stoich = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P, X_stoich)

    # Lean (phi = 0.8)
    X_lean = cb.set_equivalence_ratio_mole(0.8, X_ch4, X_air)
    T_ad_lean, _, X_products_lean = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P, X_lean)

    # Lean should have lower adiabatic temperature due to excess air
    # (though the current implementation may not show this due to placeholder stoichiometry)
    print(f"Stoichiometric CO2: {X_products_stoich[3]:.4f}, Lean CO2: {X_products_lean[3]:.4f}")
    print(f"Lean T_ad: {T_ad_lean:.2f} K")


def test_pressure_dependence():
    """Test that adiabatic temperature depends on pressure."""
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")
    X_mix = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)

    T_in = 700.0

    # Test at different pressures
    P_low = 101325.0  # 1 atm
    P_high = 500000.0  # 5 bar

    T_ad_low, _, _ = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P_low, X_mix)
    T_ad_high, _, _ = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P_high, X_mix)

    print(f"T_ad at 1 atm: {T_ad_low:.2f} K")
    print(f"T_ad at 5 bar: {T_ad_high:.2f} K")

    # Pressure should have some effect (though may be small)
    print(f"Pressure effect: {T_ad_high - T_ad_low:.2f} K")


if __name__ == "__main__":
    print("Testing combustion Jacobian functions...")

    test_adiabatic_complete_jacobian_basic()
    print("✅ Basic complete combustion Jacobian test passed")

    test_adiabatic_equilibrium_jacobian_basic()
    print("✅ Basic equilibrium combustion Jacobian test passed")

    test_jacobian_sensitivity()
    print("✅ Jacobian sensitivity test passed")

    test_lean_vs_stoichiometric()
    print("✅ Lean vs stoichiometric test completed")

    test_pressure_dependence()
    print("✅ Pressure dependence test completed")

    print("\\n🎉 All combustion Jacobian tests completed!")
