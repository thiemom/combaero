import pytest
import combaero as cb
import numpy as np

def test_phi_stoich():
    # Methane (CH4) + Air
    # CH4 + 2 O2 -> CO2 + 2 H2O
    # Stoichiometric air: O2_req = 2.0 mol/mol_fuel
    
    # Standard Air: ~21% O2, ~79% N2
    # In 1 mol of mixture at phi=1:
    # 1 mol CH4 requires 2 mol O2.
    # Total O2 needed = 2.0.
    # If we have 2 mol O2, we are stoich.
    
    # Composition mapping
    X_air = cb.species.dry_air()
    X_fuel = cb.species.pure_species("CH4")
    
    phi_target = 1.0
    # X_mix = (phi * r_st * X_fuel + X_ox) / (phi * r_st + 1)
    X_mix = cb.set_equivalence_ratio_mole(phi_target, X_fuel, X_air)
    
    phi_calc = cb.equivalence_ratio(X_mix)
    assert pytest.approx(phi_calc, rel=1e-5) == phi_target

def test_phi_rich_lean():
    X_air = cb.species.dry_air()
    X_fuel = cb.species.pure_species("CH4")
    
    for phi_target in [0.5, 0.8, 1.2, 2.0]:
        X_mix = cb.set_equivalence_ratio_mole(phi_target, X_fuel, X_air)
        phi_calc = cb.equivalence_ratio(X_mix)
        assert pytest.approx(phi_calc, rel=1e-5) == phi_target

def test_phi_mass():
    Y_air = cb.species.dry_air_mass()
    Y_fuel = cb.species.to_mass(cb.species.pure_species("CH4"))
    
    for phi_target in [0.7, 1.0, 1.5]:
        Y_mix = cb.set_equivalence_ratio_mass(phi_target, Y_fuel, Y_air)
        phi_calc = cb.equivalence_ratio_mass(Y_mix)
        assert pytest.approx(phi_calc, rel=1e-5) == phi_target

def test_phi_exceptions():
    # Pure air (no fuel)
    X_air = cb.species.dry_air()
    with pytest.raises(RuntimeError, match="No combustible species"):
        cb.equivalence_ratio(X_air)
        
    # Pure fuel (no O2)
    X_fuel = cb.species.pure_species("CH4")
    with pytest.raises(RuntimeError, match="No O2"):
        cb.equivalence_ratio(X_fuel)

if __name__ == "__main__":
    # Manually run a quick check
    print("Methane Air Phi=1 (Molar):", cb.equivalence_ratio(cb.set_equivalence_ratio_mole(1.0, cb.species.pure_species("CH4"), cb.species.dry_air())))
    print("Methane Air Phi=0.5 (Mass):", cb.equivalence_ratio_mass(cb.set_equivalence_ratio_mass(0.5, cb.species.to_mass(cb.species.pure_species("CH4")), cb.species.dry_air_mass())))
