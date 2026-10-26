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


class TestSmrWgsEquilibrium:
    """Tests for SMR+WGS equilibrium (CH4 only)."""

    def test_smr_wgs_isothermal(self, sp: SpeciesLocator) -> None:
        """Test isothermal SMR+WGS equilibrium."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        result = ca.smr_wgs_equilibrium(2000.0, X)

        # CH4 should be reformed
        assert result.state.X[sp.indices["CH4"]] < X[sp.indices["CH4"]] * 0.5

        # CO and H2 should be produced
        assert result.state.X[sp.indices["CO"]] > 0.01
        assert result.state.X[sp.indices["H2"]] > 0.01

    def test_smr_wgs_adiabatic(self, sp: SpeciesLocator) -> None:
        """Test adiabatic SMR+WGS equilibrium."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        T_in = 2200.0
        result = ca.smr_wgs_equilibrium_adiabatic(T_in, X)

        # Temperature should decrease
        assert T_in > result.state.T

    def test_smr_wgs_fallback_to_wgs(self, sp: SpeciesLocator) -> None:
        """Test that SMR+WGS falls back to WGS when no CH4 present."""
        X = sp.empty()
        X[sp.indices["CO"]] = 0.05
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.05
        X[sp.indices["N2"]] = 0.70

        result = ca.smr_wgs_equilibrium(2000.0, X)

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
    @pytest.mark.parametrize("T", [600.0, 650.0, 712.0, 800.0, 900.0])
    def test_matches_the_smr_wgs_solver_where_it_holds(self, T: float) -> None:
        """CH4 + steam (forward reforming): smr_wgs solves both extents at
        once, so the two must agree. A single round was 0.04 off at 712 K."""
        import warnings

        X = _mix(CH4=0.1, H2O=0.3, N2=0.6)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            ref = ca.smr_wgs_equilibrium(T, X)
        got = ca.reforming_equilibrium(T, X)
        assert got.converged
        assert pytest.approx(ref.state.X, abs=1e-6) == got.state.X

    @pytest.mark.parametrize("T_in", [700.0, 900.0, 1200.0])
    def test_adiabatic_matches_the_smr_wgs_solver(self, T_in: float) -> None:
        """At 900 K inlet a single round gave 734.1 K; Cantera gives 711.9."""
        import warnings

        X = _mix(CH4=0.1, H2O=0.3, N2=0.6)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            ref = ca.smr_wgs_equilibrium_adiabatic(T_in, X)
        got = ca.reforming_equilibrium_adiabatic(T_in, X)
        assert got.converged
        assert pytest.approx(ref.state.T, abs=1e-3) == got.state.T

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

    def test_a_steam_starved_rich_burn_reports_that_it_did_not_converge(self) -> None:
        """Known limit: rich propane at phi 2.5 leaves no steam after complete
        combustion, and round-by-round equilibration crawls through reverse
        WGS. The result is flagged instead of reported as converged."""
        X = _mix(C3H8=0.09, O2=0.18, N2=0.73)
        assert not ca.combustion_equilibrium(500.0, X, 1e5).converged
