import numpy as np

import combaero as cb


def test_combustion_smoothing_impact_and_jacobian():
    """
    Verify smoothing impact on temperature (< 1K) and Jacobian continuity.
    Tests across equivalence ratio range [0, 2.0] including kinks at 0 and 1.
    """
    n_sp = cb.species.num_species
    X_fuel = np.zeros(n_sp)
    X_fuel[cb.species.species_index("CH4")] = 1.0

    X_ox = np.zeros(n_sp)
    X_ox[cb.species.species_index("O2")] = 0.21
    X_ox[cb.species.species_index("N2")] = 0.79

    T_reactants = 300.0
    P = 101325.0

    # Increase resolution around the kink at Phi=1
    phis = np.linspace(1e-9, 2.0, 10001)

    t_no_smooth = []
    t_smooth = []

    for phi in phis:
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_ox)
        s_raw = cb.complete_combustion(T_reactants, X_mix, P, smooth_phi0=False, smooth_phi1=False)
        t_no_smooth.append(s_raw.T)

        s_sm = cb.complete_combustion(T_reactants, X_mix, P, smooth_phi0=True, smooth_phi1=True)
        t_smooth.append(s_sm.T)

    t_no_smooth = np.array(t_no_smooth)
    t_smooth = np.array(t_smooth)
    diff = np.abs(t_no_smooth - t_smooth)

    # 1. VERIFY IMPACT: Max temperature difference should be < 1K with k=20000
    max_diff = np.max(diff)
    print(f"\nMax temperature difference with smoothing (k=20000): {max_diff:.4f} K")
    assert max_diff < 1.0, f"Smoothing impact too large: {max_diff:.2f} K > 1.0 K"

    # 2. VERIFY SMOOTHNESS: Max curvature (numerical second derivative)
    # We check both raw and smoothed.
    # dT/dphi jump at Phi=1 in complete combustion is inherent.
    # Smoothing should reduce the magnitude of the discrete jump.
    dT_dphi_raw = np.diff(t_no_smooth) / np.diff(phis)
    dT_dphi_sm = np.diff(t_smooth) / np.diff(phis)

    jump_raw = np.max(np.abs(np.diff(dT_dphi_raw)))
    jump_sm = np.max(np.abs(np.diff(dT_dphi_sm)))

    print(f"Max curvature (jump in dT/dphi) raw: {jump_raw:.4f}")
    print(f"Max curvature (jump in dT/dphi) smoothed: {jump_sm:.4f}")

    # At k=20000, transition is sharp relative to grid, but still finite.
    assert jump_sm < 10000.0, f"Smoothed Jacobian seems to have excessive curvature: {jump_sm:.2f}"
    assert jump_sm < jump_raw, "Smoothing didn't reduce slope discontinuity"


if __name__ == "__main__":
    test_combustion_smoothing_impact_and_jacobian()
