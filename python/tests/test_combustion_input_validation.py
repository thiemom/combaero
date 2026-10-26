"""Combustion helpers refuse invalid input instead of returning a wrong answer.

Found in the 0.7.0 combustion audit: a NaN slipped past every
`x <= 0` / `|sum - 1| > tol` guard (NaN compares false), so

  - complete_combustion_to_CO2_H2O / _products_and_jacobian returned NaN
    products, complete_combustion a calc_T_from_h convergence error;
  - set_equivalence_ratio_mole(nan or inf) returned an all-NaN mixture;
  - set_fuel_stream_for_Tad(nan) returned the STOICHIOMETRIC fuel flow
    (Brent's method with a NaN residual falls through to its bracket end),
    and every other set_*_stream_for_* target the same way;
  - fuel_lhv_molar(X, nan) returned 0.0 (std::max(0, NaN) == 0);
  - complete_combustion accepted T <= 0 and P <= 0 (T_ad = 1235 K from a
    -10 K reactant).

and CombustionState.fuel_burn_fraction read 1.0 for a rich complete burn
that leaves 1 - 1/phi of the fuel unburnt.
"""

import math

import numpy as np
import pytest

import combaero as cb

NAN = float("nan")


def _vec(**fractions: float) -> np.ndarray:
    x = np.zeros(cb.num_species())
    for name, value in fractions.items():
        x[cb.species_index_from_name(name)] = value
    return x / x.sum()


AIR = np.array(cb.species.dry_air())
CH4 = _vec(CH4=1.0)


def _stream(T: float, P: float, X: np.ndarray, mdot: float) -> cb.Stream:
    s = cb.Stream()
    s.T = T
    s.P = P
    s.X = X
    s.mdot = mdot
    return s


def _nan_mixture() -> np.ndarray:
    x = np.array(cb.set_equivalence_ratio_mole(0.8, CH4, AIR))
    x[cb.species_index_from_name("CH4")] = NAN
    return x


def test_nan_mole_fractions_are_refused() -> None:
    x = _nan_mixture()
    with pytest.raises(ValueError, match="finite"):
        cb.complete_combustion_to_CO2_H2O(x)
    with pytest.raises(ValueError, match="finite"):
        cb.complete_combustion_products_and_jacobian(x)
    with pytest.raises(ValueError, match="finite"):
        cb.complete_combustion_to_CO2_H2O_with_fraction(x)


def test_nan_mole_fractions_are_refused_by_the_stoichiometry_helpers() -> None:
    x = _nan_mixture()
    with pytest.raises(ValueError, match="finite"):
        cb.equivalence_ratio(x)
    with pytest.raises(ValueError, match="finite"):
        cb.oxygen_required_per_mol_mixture(x)
    with pytest.raises(ValueError, match="finite"):
        cb.oxygen_required_per_kg_mixture(x)
    fuel = CH4.copy()
    fuel[cb.species_index_from_name("N2")] = NAN
    with pytest.raises(ValueError, match="finite"):
        cb.set_equivalence_ratio_mole(1.0, fuel, AIR)


@pytest.mark.parametrize("fn", [cb.complete_combustion, cb.complete_combustion_isothermal])
@pytest.mark.parametrize(
    ("T", "P"), [(-10.0, 1e5), (0.0, 1e5), (NAN, 1e5), (300.0, 0.0), (300.0, -1e5), (300.0, NAN)]
)
def test_complete_combustion_refuses_non_physical_state(fn, T: float, P: float) -> None:
    x = np.array(cb.set_equivalence_ratio_mole(0.5, CH4, AIR))
    with pytest.raises(ValueError):
        fn(T, x, P)


@pytest.mark.parametrize("phi", [NAN, math.inf])
def test_non_finite_phi_is_refused(phi: float) -> None:
    y_f = np.array(cb.mole_to_mass(CH4))
    y_o = np.array(cb.mole_to_mass(AIR))
    with pytest.raises(ValueError):
        cb.set_equivalence_ratio_mole(phi, CH4, AIR)
    with pytest.raises(ValueError):
        cb.set_equivalence_ratio_mass(phi, y_f, y_o)
    with pytest.raises(ValueError):
        cb.bilger_Z_from_equivalence_ratio_mass(phi, y_f, y_o)
    with pytest.raises(ValueError):
        cb.set_fuel_stream_for_phi(
            phi, _stream(300.0, 1e5, CH4, 0.0), _stream(300.0, 1e5, AIR, 1.0)
        )


def test_nan_mixture_fraction_is_refused() -> None:
    y_f = np.array(cb.mole_to_mass(CH4))
    y_o = np.array(cb.mole_to_mass(AIR))
    with pytest.raises(ValueError):
        cb.equivalence_ratio_from_bilger_Z_mass(NAN, y_f, y_o)


@pytest.mark.parametrize(
    "fn",
    [
        cb.set_fuel_stream_for_Tad,
        cb.set_fuel_stream_for_O2,
        cb.set_fuel_stream_for_O2_dry,
        cb.set_fuel_stream_for_CO2,
        cb.set_fuel_stream_for_CO2_dry,
        cb.set_oxidizer_stream_for_Tad,
        cb.set_oxidizer_stream_for_O2,
        cb.set_oxidizer_stream_for_O2_dry,
        cb.set_oxidizer_stream_for_CO2,
        cb.set_oxidizer_stream_for_CO2_dry,
    ],
)
def test_nan_stream_target_is_refused(fn) -> None:
    """A NaN target returned the stoichiometric flow before (a valid-looking answer)."""
    fuel = _stream(300.0, 1e5, CH4, 0.05)
    ox = _stream(600.0, 1e5, AIR, 1.0)
    with pytest.raises(ValueError):
        fn(NAN, fuel, ox)


def test_fuel_lhv_refuses_invalid_input() -> None:
    with pytest.raises(ValueError):
        cb.fuel_lhv_molar(CH4, NAN)
    with pytest.raises(ValueError):
        cb.fuel_lhv_mass(CH4, NAN)
    x = np.zeros(cb.num_species())
    x[cb.species_index_from_name("CH4")] = 1.1
    x[cb.species_index_from_name("N2")] = -0.1
    with pytest.raises(ValueError):
        cb.fuel_lhv_molar(x)
    x = CH4.copy()
    x[cb.species_index_from_name("N2")] = NAN
    with pytest.raises(ValueError):
        cb.fuel_lhv_molar(x)


@pytest.mark.parametrize("phi", [0.6, 1.0, 1.2, 2.0])
def test_fuel_burn_fraction_is_the_burnt_share_of_the_fuel(phi: float) -> None:
    """Complete combustion burns min(1, 1/phi) of the fuel; the rest is left in
    the products (it read 1.0 at every phi before)."""
    st = cb.combustion_state(CH4, AIR, phi, 300.0, 1e5)
    assert st.fuel_burn_fraction == pytest.approx(min(1.0, 1.0 / phi), rel=1e-12)
    i = cb.species_index_from_name("CH4")
    x_in = np.array(cb.set_equivalence_ratio_mole(phi, CH4, AIR))
    x_b = np.array(cb.complete_combustion_to_CO2_H2O(x_in))
    # unburnt CH4 per mole N2 (N2 is inert here)
    n2 = cb.species_index_from_name("N2")
    left = (x_b[i] / x_b[n2]) / (x_in[i] / x_in[n2])
    assert left == pytest.approx(1.0 - st.fuel_burn_fraction, abs=1e-12)


def test_combustion_state_from_streams_refuses_no_oxidizer_flow() -> None:
    with pytest.raises(ValueError):
        cb.combustion_state_from_streams(
            _stream(300.0, 1e5, CH4, 0.05), _stream(300.0, 1e5, AIR, 0.0)
        )
    with pytest.raises(ValueError):
        cb.combustion_state_from_streams(
            _stream(300.0, 1e5, CH4, -0.01), _stream(300.0, 1e5, AIR, 1.0)
        )


def test_combustion_state_from_streams_uses_the_mixed_pressure() -> None:
    """The reactant temperature comes from mix() at the lowest inlet pressure;
    the products were put at the FUEL stream's pressure (3 bar here)."""
    fuel = _stream(300.0, 3e5, CH4, 0.03)
    ox = _stream(600.0, 1e5, AIR, 1.0)
    st = cb.combustion_state_from_streams(fuel, ox)
    p_mix = cb.mix([fuel, ox]).P
    assert pytest.approx(p_mix, rel=1e-15) == st.reactants.thermo.P
    assert pytest.approx(p_mix, rel=1e-15) == st.products.thermo.P
