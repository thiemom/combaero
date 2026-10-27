"""Standalone gas-dynamics helpers against closed forms and their own inverses.

Argon is calorically perfect in the thermo backend (cp = 5/2 R), so the
constant-gamma relations with gamma = 5/3 are exact references for it.
"""

import math

import pytest

import combaero as cb

GAMMA_AR = 5.0 / 3.0


def _argon() -> list[float]:
    X = [0.0] * cb.num_species()
    X[cb.species_index_from_name("AR")] = 1.0
    return X


def _fanno_fL_star(M: float, g: float) -> float:
    """Darcy f L*/D (= 4 f_Fanning L*/D) for a calorically perfect gas."""
    M2 = M * M
    return (1.0 - M2) / (g * M2) + (g + 1.0) / (2.0 * g) * math.log(
        (g + 1.0) * M2 / (2.0 + (g - 1.0) * M2)
    )


@pytest.mark.parametrize("M_in", [0.1, 0.3, 0.5, 0.7, 0.9, 0.98])
def test_fanno_max_length_matches_closed_form(M_in: float) -> None:
    """L* is the length to M = 1, not to the march's 0.99 / 0.999 cutoff."""
    X = _argon()
    T, P, D, f = 300.0, 2.0e5, 0.02, 0.02
    u = M_in * cb.speed_of_sound(T, X)
    L = cb.fanno_max_length(T, P, u, D, f, X)
    L_ref = _fanno_fL_star(M_in, GAMMA_AR) * D / f
    assert pytest.approx(L_ref, rel=1e-7) == L


@pytest.mark.parametrize("u_in", [30.0, 120.0, 250.0])
def test_fanno_max_length_agrees_with_march_choke(u_in: float) -> None:
    """Variable-cp air: L* agrees with where fanno_channel itself chokes."""
    X = cb.species.dry_air()
    T, P, D, f = 400.0, 2.0e5, 0.05, 0.02
    L_star = cb.fanno_max_length(T, P, u_in, D, f, X)
    sol = cb.fanno_channel(T, P, u_in, 1.5 * L_star, D, f, X, n_steps=4000)
    assert sol.choked
    # The march stops at M = 0.999, (1 - M)^2-close to L*; at the default
    # 100 steps its own L_choke is only good to ~0.5% at M_in ~ 0.07.
    assert sol.L_choke == pytest.approx(L_star, rel=1e-5)
    assert not cb.fanno_channel(T, P, u_in, 0.99 * L_star, D, f, X, n_steps=400).choked


@pytest.mark.parametrize(
    "args",
    [
        (300.0, 2e5, 0.0, 0.05, 0.02),  # no flow: L* infinite
        (300.0, 2e5, -50.0, 0.05, 0.02),
        (300.0, 2e5, 50.0, 0.05, 0.0),  # frictionless: L* infinite
        (300.0, 2e5, 50.0, 0.0, 0.02),
        (-300.0, 2e5, 50.0, 0.05, 0.02),
        (300.0, -2e5, 50.0, 0.05, 0.02),
    ],
)
def test_fanno_max_length_rejects_invalid(args: tuple[float, ...]) -> None:
    with pytest.raises(ValueError):
        cb.fanno_max_length(*args, cb.species.dry_air())


@pytest.mark.parametrize(
    ("mdot", "rel"),
    # 1e-6 kg/s is a 2e-5 Pa drive on 1 bar: the forward h0 - h difference is
    # itself only good to ~1e-6 there.
    [
        (1e-6, 1e-5),
        (1e-4, 1e-6),
        (1e-3, 1e-6),
        (1e-2, 1e-6),
        (0.1, 1e-6),
        (5.0, 1e-6),
        (500.0, 1e-6),
    ],
)
def test_solve_P0_from_mdot_round_trip(mdot: float, rel: float) -> None:
    """Small flows need P0 within 0.1% of P_back; the bracket must reach them."""
    X = cb.species.dry_air()
    T0, P_back, A = 600.0, 1.0e5, 2.0e-4
    P0 = cb.solve_Pt_from_mdot(T0, P_back, A, mdot, X)
    assert P_back < P0
    assert cb.nozzle_flow(T0, P0, P_back, A, X).mdot == pytest.approx(mdot, rel=rel)


@pytest.mark.parametrize("frac", [1e-4, 1e-2, 0.3, 0.7, 0.95, 0.999, 0.99999])
def test_solve_P_back_from_mdot_round_trip(frac: float) -> None:
    """The back pressure reproduces the target flow, including tiny drops."""
    X = cb.species.dry_air()
    T0, P0, A = 600.0, 4.0e5, 2.0e-4
    m_choked = cb.nozzle_flow(T0, P0, 1.0e3, A, X).mdot
    Pb = cb.solve_P_back_from_mdot(T0, P0, A, frac * m_choked, X)
    assert Pb < P0
    assert cb.nozzle_flow(T0, P0, Pb, A, X).mdot == pytest.approx(frac * m_choked, rel=1e-6)


def test_solve_P_back_from_mdot_refuses_above_choked() -> None:
    X = cb.species.dry_air()
    T0, P0, A = 600.0, 4.0e5, 2.0e-4
    m_choked = cb.nozzle_flow(T0, P0, 1.0e3, A, X).mdot
    with pytest.raises(RuntimeError, match="exceeds choked"):
        cb.solve_P_back_from_mdot(T0, P0, A, 1.001 * m_choked, X)


def test_solve_inverse_functions_refuse_invalid_input() -> None:
    X = cb.species.dry_air()
    with pytest.raises(ValueError):
        cb.solve_A_eff_from_mdot(600.0, 4e5, 2e5, -0.05, X)
    with pytest.raises(ValueError):
        cb.solve_A_eff_from_mdot(600.0, 4e5, 2e5, 0.0, X)
    with pytest.raises(ValueError):
        cb.solve_A_eff_from_mdot(-600.0, 4e5, 2e5, 0.05, X)
    with pytest.raises(ValueError):
        cb.solve_Pt_from_mdot(-600.0, 1e5, 2e-4, 0.05, X)
    with pytest.raises(ValueError):
        cb.solve_Pt_from_mdot(600.0, -1e5, 2e-4, 0.05, X)
    with pytest.raises(ValueError):
        cb.solve_Pt_from_mdot(600.0, 1e5, 0.0, 0.05, X)
    with pytest.raises(ValueError):
        cb.solve_P_back_from_mdot(-600.0, 4e5, 2e-4, 0.05, X)
    with pytest.raises(ValueError):
        cb.solve_P_back_from_mdot(600.0, -4e5, 2e-4, 0.05, X)


@pytest.mark.parametrize("M", [0.1, 0.5, 1.0, 2.0])
def test_isentropic_relations_constant_gamma_limit(M: float) -> None:
    X = _argon()
    g = GAMMA_AR
    T, P = 300.0, 1.0e5
    T0_ref = T * (1.0 + 0.5 * (g - 1.0) * M * M)
    P0_ref = P * (T0_ref / T) ** (g / (g - 1.0))
    assert cb.Tt_from_static(T, M, X) == pytest.approx(T0_ref, rel=1e-12)
    assert cb.Pt_from_static(P, T, M, X) == pytest.approx(P0_ref, rel=1e-12)
    assert cb.T_from_stagnation(T0_ref, M, X) == pytest.approx(T, rel=1e-12)
    assert cb.P_from_stagnation(P0_ref, T0_ref, M, X) == pytest.approx(P, rel=1e-12)
    assert cb.mach_from_pressure_ratio(T0_ref, P0_ref, P, X) == pytest.approx(M, rel=1e-10)


def test_choked_flux_constant_gamma_limit() -> None:
    X = _argon()
    g = GAMMA_AR
    R = cb.specific_gas_constant(X)
    T0, P0 = 500.0, 3.0e5
    pr_ref = (2.0 / (g + 1.0)) ** (g / (g - 1.0))
    G_ref = P0 * math.sqrt(g / (R * T0)) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    assert cb.critical_pressure_ratio(T0, P0, X) == pytest.approx(pr_ref, rel=1e-7)
    assert cb.nozzle_flow(T0, P0, 1.0e4, 1.0, X).mdot == pytest.approx(G_ref, rel=1e-12)


def test_stagnation_helpers_documented_in_api_are_exported() -> None:
    """docs/API_PYTHON.md calls cb.bulk_velocity and cb.stagnation_from_static."""
    air = cb.species.dry_air()
    v = cb.bulk_velocity(m_dot=0.5, rho=1.2, area=1e-3)
    assert v == pytest.approx(0.5 / (1.2 * 1e-3))
    Tt, Pt = cb.stagnation_from_static(T=300.0, P=101325.0, v=v, X=air)
    M = cb.mach_number(v, 300.0, air)
    assert Tt == pytest.approx(cb.Tt_from_static(300.0, M, air), rel=1e-12)
    assert Pt == pytest.approx(cb.Pt_from_static(101325.0, 300.0, M, air), rel=1e-12)


def test_incompressible_channel_profile_enthalpy_is_mass_basis() -> None:
    """IncompressibleStation.h is documented [J/kg], like FannoStation.h."""
    X = cb.species.dry_air()
    T, P = 300.0, 1.2e5
    sol = cb.channel_flow(T, P, X, 20.0, 2.0, 0.05, 0.02, 5, True)
    assert sol.profile[0].h == pytest.approx(cb.h_mass(T, X), rel=1e-12)
