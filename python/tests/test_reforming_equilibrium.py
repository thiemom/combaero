"""Tests for reforming equilibrium functions."""

import pytest

import combaero as ca
from combaero.species import SpeciesLocator


@pytest.fixture
def sp() -> SpeciesLocator:
    """Species locator fixture."""
    return SpeciesLocator.from_core()


class TestReformingEquilibrium:
    """Tests for reforming_equilibrium and reforming_equilibrium_adiabatic."""

    def test_reforming_equilibrium_isothermal(self, sp: SpeciesLocator) -> None:
        """Test isothermal reforming equilibrium with CH4."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        T_in = 2000.0
        result = ca.reforming_equilibrium(T_in, X)

        # Temperature should be unchanged (isothermal)
        assert pytest.approx(T_in, rel=1e-6) == result.state.T

        # CH4 should be mostly reformed at high T
        assert result.state.X[sp.indices["CH4"]] < X[sp.indices["CH4"]] * 0.5

        # CO and H2 should be produced
        assert result.state.X[sp.indices["CO"]] > 0.01
        assert result.state.X[sp.indices["H2"]] > 0.01

    def test_reforming_equilibrium_adiabatic(self, sp: SpeciesLocator) -> None:
        """Test adiabatic reforming equilibrium - temperature should drop."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        T_in = 2200.0
        result = ca.reforming_equilibrium_adiabatic(T_in, X)

        # Temperature should decrease (endothermic reforming)
        assert T_in > result.state.T
        assert result.state.T > 1500.0  # But still reasonable

        # CH4 should be reformed
        assert result.state.X[sp.indices["CH4"]] < X[sp.indices["CH4"]] * 0.5

    def test_reforming_multiple_hydrocarbons(self, sp: SpeciesLocator) -> None:
        """Test reforming with multiple hydrocarbons (natural gas)."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.015
        X[sp.indices["C2H6"]] = 0.003
        X[sp.indices["C3H8"]] = 0.002
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        result = ca.reforming_equilibrium(2000.0, X)

        # All hydrocarbons should be reformed
        assert result.state.X[sp.indices["CH4"]] < X[sp.indices["CH4"]] * 0.5
        assert result.state.X[sp.indices["C2H6"]] < X[sp.indices["C2H6"]] * 0.5
        assert result.state.X[sp.indices["C3H8"]] < X[sp.indices["C3H8"]] * 0.5

        # CO and H2 should be produced
        assert result.state.X[sp.indices["CO"]] > 0.01
        assert result.state.X[sp.indices["H2"]] > 0.01

    def test_reforming_c2plus_only(self, sp: SpeciesLocator) -> None:
        """Test reforming with only C2+ hydrocarbons (no CH4)."""
        X = sp.empty()
        X[sp.indices["C2H6"]] = 0.010
        X[sp.indices["C3H8"]] = 0.005
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.085
        X[sp.indices["N2"]] = 0.70

        result = ca.reforming_equilibrium(2000.0, X)

        # C2H6 and C3H8 should be reformed even without CH4
        assert result.state.X[sp.indices["C2H6"]] < X[sp.indices["C2H6"]] * 0.5
        assert result.state.X[sp.indices["C3H8"]] < X[sp.indices["C3H8"]] * 0.5

        # CO and H2 should be produced
        assert result.state.X[sp.indices["CO"]] > 0.01
        assert result.state.X[sp.indices["H2"]] > 0.01

    def test_element_conservation(self, sp: SpeciesLocator) -> None:
        """Test that element ratios are conserved (using N2 as reference)."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.015
        X[sp.indices["C2H6"]] = 0.003
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70
        X[sp.indices["AR"]] = 0.002

        result = ca.reforming_equilibrium(2000.0, X)

        # N2 is inert, so element ratios per N2 should be conserved
        n2_in = X[sp.indices["N2"]]
        n2_out = result.state.X[sp.indices["N2"]]

        # C atoms per N2: CH4 + 2*C2H6 + CO + CO2
        C_per_N2_in = (
            X[sp.indices["CH4"]]
            + 2 * X[sp.indices["C2H6"]]
            + X[sp.indices["CO"]]
            + X[sp.indices["CO2"]]
        ) / n2_in
        C_per_N2_out = (
            result.state.X[sp.indices["CH4"]]
            + 2 * result.state.X[sp.indices["C2H6"]]
            + result.state.X[sp.indices["CO"]]
            + result.state.X[sp.indices["CO2"]]
        ) / n2_out
        assert C_per_N2_in == pytest.approx(C_per_N2_out, rel=1e-5)

        # Ar per N2 should be unchanged
        ar_per_n2_in = X[sp.indices["AR"]] / n2_in
        ar_per_n2_out = result.state.X[sp.indices["AR"]] / n2_out
        assert ar_per_n2_in == pytest.approx(ar_per_n2_out, rel=1e-10)


class TestCh4OnlyReformingEquilibrium:
    """CH4-only reforming (the cases of the removed smr_wgs_equilibrium)."""

    def test_ch4_isothermal(self, sp: SpeciesLocator) -> None:
        """Isothermal CH4 reforming + WGS."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        result = ca.reforming_equilibrium(2000.0, X)

        # CH4 should be reformed
        assert result.state.X[sp.indices["CH4"]] < X[sp.indices["CH4"]] * 0.5

        # CO and H2 should be produced
        assert result.state.X[sp.indices["CO"]] > 0.01
        assert result.state.X[sp.indices["H2"]] > 0.01

    def test_ch4_adiabatic(self, sp: SpeciesLocator) -> None:
        """Adiabatic CH4 reforming + WGS."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        T_in = 2200.0
        result = ca.reforming_equilibrium_adiabatic(T_in, X)

        # Temperature should decrease
        assert T_in > result.state.T

    def test_fallback_to_wgs(self, sp: SpeciesLocator) -> None:
        """With no hydrocarbon present, reforming reduces to WGS."""
        X = sp.empty()
        X[sp.indices["CO"]] = 0.05
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.05
        X[sp.indices["N2"]] = 0.70

        result = ca.reforming_equilibrium(2000.0, X)

        # Should still work (WGS only)
        assert pytest.approx(2000.0, rel=1e-6) == result.state.T
        # WGS shifts CO + H2O -> CO2 + H2
        assert result.state.X[sp.indices["H2"]] > 0.0


# ---------------------------------------------------------------------------
# Simultaneous equilibrium: the reactions are re-equilibrated in rounds until
# nothing moves. A single round (reform, then WGS once) left H2 7-49% low at
# 650-850 K, the pre-reformer range: WGS moved CO, H2 and H2O after reforming
# had equilibrated. Cantera comparisons: cantera_validation_tests.
# ---------------------------------------------------------------------------


def _mix(**kw: float) -> list[float]:
    import numpy as np

    X = np.zeros(ca.num_species())
    for name, v in kw.items():
        X[ca.species_index_from_name(name)] = v
    return list(X / X.sum())


def _elements(X: list[float]) -> dict[str, float]:
    import numpy as np

    counts = {
        "C": {"CH4": 1, "C2H6": 2, "C3H8": 3, "CO": 1, "CO2": 1},
        "H": {"CH4": 4, "C2H6": 6, "C3H8": 8, "H2": 2, "H2O": 2},
        "O": {"CO": 1, "CO2": 2, "H2O": 1},
    }
    n = np.asarray(X)
    return {
        e: sum(k * n[ca.species_index_from_name(s)] for s, k in sp_.items())
        / n[ca.species_index_from_name("N2")]
        for e, sp_ in counts.items()
    }


class TestSimultaneousReformingEquilibrium:
    # Golden values from the two-extent SMR+WGS Newton solver
    # (smr_wgs_equilibrium, removed in 0.8.0), evaluated at v0.7.1 just before
    # its removal: it solved both extents at once, so the simultaneous rounds
    # must reproduce it. Mole fractions of CH4, H2O, CO, CO2, H2.
    _SMR_WGS_X = {
        600.0: (0.084596974, 0.266713558, 0.000086780, 0.012749074, 0.051256639),
        650.0: (0.074570404, 0.245360322, 0.000457781, 0.020733550, 0.084307540),
        712.0: (0.058178650, 0.211869097, 0.002482021, 0.032369103, 0.136922477),
        800.0: (0.029704814, 0.161972378, 0.014278614, 0.044300708, 0.220038672),
        900.0: (0.004793335, 0.131784917, 0.038066025, 0.041272863, 0.279289526),
    }
    _SMR_WGS_T_AD = {700.0: 617.762887, 900.0: 711.998293, 1200.0: 824.688325}

    @pytest.mark.parametrize("T", [600.0, 650.0, 712.0, 800.0, 900.0])
    def test_matches_the_smr_wgs_solver_where_it_holds(self, T: float) -> None:
        """CH4 + steam (forward reforming): the removed smr_wgs solver solved
        both extents at once, so the two must agree. A single round was 0.04
        off at 712 K."""
        X = _mix(CH4=0.1, H2O=0.3, N2=0.6)
        got = ca.reforming_equilibrium(T, X)
        assert got.converged
        names = ("CH4", "H2O", "CO", "CO2", "H2")
        got_X = [got.state.X[ca.species_index_from_name(n)] for n in names]
        assert pytest.approx(self._SMR_WGS_X[T], abs=1e-6) == got_X

    @pytest.mark.parametrize("T_in", [700.0, 900.0, 1200.0])
    def test_adiabatic_matches_the_smr_wgs_solver(self, T_in: float) -> None:
        """At 900 K inlet a single round gave 734.1 K; Cantera gives 711.9."""
        X = _mix(CH4=0.1, H2O=0.3, N2=0.6)
        got = ca.reforming_equilibrium_adiabatic(T_in, X)
        assert got.converged
        assert pytest.approx(self._SMR_WGS_T_AD[T_in], abs=1e-3) == got.state.T

    @pytest.mark.parametrize(
        "mix",
        [
            {"CH4": 0.05, "C2H6": 0.03, "C3H8": 0.02, "H2O": 0.4, "N2": 0.5},
            {"CH4": 0.02, "CO": 0.2, "H2": 0.3, "H2O": 0.1, "CO2": 0.05, "N2": 0.33},
        ],
        ids=["higher-hydrocarbons", "syngas-methanation"],
    )
    @pytest.mark.parametrize("T", [650.0, 800.0])
    def test_conserves_the_elements(self, mix: dict[str, float], T: float) -> None:
        X = _mix(**mix)
        got = ca.reforming_equilibrium(T, X)
        assert got.converged
        e0, e1 = _elements(X), _elements(list(got.state.X))
        for el in e0:
            assert e1[el] == pytest.approx(e0[el], rel=1e-10)

    def test_higher_hydrocarbons_become_methane_at_low_temperature(self) -> None:
        """At 650 K equilibrium holds almost no ethane or propane: they reform,
        and CH4 re-forms from CO + H2 (methanation, a backward reforming
        extent). Forward-only extents left them 0.2 off in mole fraction."""
        X = _mix(CH4=0.05, C2H6=0.03, C3H8=0.02, H2O=0.4, N2=0.5)
        got = ca.reforming_equilibrium(650.0, X)
        for name in ("C2H6", "C3H8"):
            assert got.state.X[ca.species_index_from_name(name)] < 1e-6
        assert got.state.X[ca.species_index_from_name("CH4")] > 0.05


# ---------------------------------------------------------------------------
# Steam-starved rich burn (#759). Complete combustion of rich propane leaves
# steam that the propane reforming then uses up; the rounds (#755) got it back
# only through reverse WGS, a sliver per round, and stopped 2000 rounds short:
# phi 2.5-3.0 ended 124-320 K above Cantera, converged=False. The extents are
# now solved jointly (Newton, element-potential restart). Cantera gri30 HP
# temperatures recorded 2026-10-28 (cantera_validation_tests compares live).
# ---------------------------------------------------------------------------


def _rich(fuel: str, phi: float) -> list[float]:
    import numpy as np

    X_fuel = np.zeros(ca.num_species())
    X_fuel[ca.species_index_from_name(fuel)] = 1.0
    return list(ca.set_equivalence_ratio_mole(phi, X_fuel, ca.species.dry_air()))


_C3H8_RICH_CANTERA_T = {2.5: 1336.72, 2.75: 1200.50, 3.0: 1080.77}


@pytest.mark.parametrize("phi", [2.5, 2.75, 3.0])
def test_steam_starved_rich_propane_converges_to_cantera(phi: float) -> None:
    """Before #759: converged=False and T 124 / 229 / 320 K hot; now within
    1.5 K (the rest is the NASA-7 vs gri30 thermo and gri30's minor species)."""
    got = ca.combustion_equilibrium(300.0, _rich("C3H8", phi), 101325.0)
    assert got.converged
    assert pytest.approx(_C3H8_RICH_CANTERA_T[phi], abs=3.0) == got.state.T
    # The steam is back: Cantera has 0.048 / 0.030 / 0.015, the rounds left 0.
    assert got.state.X[ca.species_index_from_name("H2O")] > 0.01


def test_the_steam_starved_issue_mixture_converges() -> None:
    """The #755 report's case: C3H8 0.09, O2 0.18, N2 0.73 (phi 2.5) at 500 K.
    Cantera gri30 HP: 1446.38 K; the rounds gave 1559.9 K, unconverged."""
    X = _mix(C3H8=0.09, O2=0.18, N2=0.73)
    got = ca.combustion_equilibrium(500.0, X, 1e5)
    assert got.converged
    assert pytest.approx(1446.38, abs=3.0) == got.state.T


@pytest.mark.parametrize("phi", [2.5, 3.0])
def test_the_steam_starved_jacobian_matches_central_differences(phi: float) -> None:
    """combustion_equilibrium_and_jacobian is one implicit-function solve on
    the converged joint system. Before #759 it chained the unconverged rounds:
    dT/dT_in 1.5% off, dX_CO 10x. Two step sizes each; measured <= 3e-7 (P) and <= 6e-8 (T_in, X)."""
    import numpy as np

    X = np.array(_rich("C3H8", phi))
    T, P = 300.0, 101325.0
    j = ca.combustion_equilibrium_and_jacobian(T, X, P)
    idx = [ca.species_index_from_name(n) for n in ("CO", "H2O", "H2", "CO2", "CH4")]
    v = np.zeros(len(X))
    v[ca.species_index_from_name("C3H8")], v[ca.species_index_from_name("N2")] = 1.0, -1.0

    def val(T_: float = T, P_: float = P, d: float = 0.0) -> np.ndarray:
        r = ca.combustion_equilibrium(T_, X + d * v, P_)
        return np.array([r.state.T] + [r.state.X[i] for i in idx])

    dX_dX = np.array(j.dXout_dX)
    analytic = {
        "T_in": np.array([j.dT_dT_in] + [j.dXout_dT_in[i] for i in idx]),
        "P": np.array([j.dT_dP] + [j.dXout_dP[i] for i in idx]),
        "X": np.array([np.dot(j.dT_dX, v)] + [np.dot(dX_dX[i], v) for i in idx]),
    }
    checks = [
        ("T_in", (1e-2, 1e-3), lambda h: (val(T_=T + h) - val(T_=T - h)) / (2 * h), 1e-6),
        ("P", (100.0, 10.0), lambda h: (val(P_=P + h) - val(P_=P - h)) / (2 * h), 1e-5),
        ("X", (1e-6, 1e-7), lambda h: (val(d=h) - val(d=-h)) / (2 * h), 1e-6),
    ]
    for q, steps, fd_of, rtol in checks:
        for h in steps:
            fd = fd_of(h)
            err = np.abs(analytic[q] - fd) / np.abs(fd)
            assert err.max() <= rtol, (q, h, err, analytic[q], fd)
