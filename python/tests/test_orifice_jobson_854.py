"""Opt-in compressible orifice models (#854): Jobson (1955) and ISO 5167-2.

``jobson1955`` is Jobson's theory for the compressible contraction of a
sharp-edged orifice (Proc. IMechE 169, 767), extended with the approach
velocity of a pipe; ``iso5167`` is ISO 5167-2:2003 Eq. (1) with Eq. (5)'s
empirical expansibility, restricted to p2/p1 >= 0.75. Since 0.8.0 (#861)
jobson1955 is the build-time default for an incompressible-referenced Cd on
a sharp-edged bore at beta <= 0.75 (test_orifice_default_model_861.py); the
nozzle path itself is unchanged.

Two kinds of evidence, kept apart (docs/VALIDATION_POLICY.md):
  * FIDELITY, on the author's own numbers: Jobson's Table 2.
  * CROSS-SOURCE accuracy: the beta extension against ISO 5167-2 Eq. (5),
    which Jobson never saw and nothing here is fitted to.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
X_AIR = list(cb.species.dry_air())
CI_SLIT = math.pi / (math.pi + 2.0)

# Jobson (1955) Table 2: the slit, Ci = pi/(pi + 2), n = 1.4, as printed.
TABLE2 = [
    (0.932, 0.622),
    (0.865, 0.635),
    (0.805, 0.648),
    (0.745, 0.663),
    (0.690, 0.679),
    (0.640, 0.696),
    (0.590, 0.714),
    (0.545, 0.734),
    (0.528, 0.741),
]


def _jobson_literal(r: float, Ci: float, n: float) -> float:
    """Eqs. (5), (8), (11), (12), (15) transcribed from the paper, by hand."""
    f = 1.0 / Ci - 1.0 / (2.0 * Ci**2)
    rc = (2.0 / (n + 1.0)) ** (n / (n - 1.0))
    if r >= rc:
        K2 = 2.0 * n / (n - 1.0) * r ** (2.0 / n) * (1.0 - r ** ((n - 1.0) / n))
        return (1.0 / (2.0 * f * r ** (1.0 / n))) * (
            1.0 - math.sqrt(1.0 - (2.0 * r ** (1.0 / n)) ** 2 * (1.0 - r) * f / K2)
        )
    K2 = n * (2.0 / (n + 1.0)) ** ((n + 1.0) / (n - 1.0))
    a = 1.0 + (rc - r) * rc ** (1.0 / n) / K2
    return (1.0 / (2.0 * f * rc ** (1.0 / n))) * (
        a - math.sqrt(a * a - (2.0 * rc ** (1.0 / n)) ** 2 * (1.0 - r) * f / K2)
    )


# ---------------------------------------------------------------------------
# Fidelity: Jobson's own Table 2
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("r", "C_paper"), TABLE2)
def test_table2_is_reproduced_within_its_printed_digits(r: float, C_paper: float) -> None:
    """FIDELITY: the paper's Table 2 (3 decimals) at every printed r."""
    assert cb.jobson_contraction_coefficient(r, CI_SLIT, 1.4) == pytest.approx(C_paper, abs=1e-3)


@pytest.mark.parametrize("r", [0.95, 0.8, 0.6, 0.528, 0.5, 0.3, 0.1])
@pytest.mark.parametrize("Ci", [0.55, 0.611, 0.7])
def test_the_library_is_the_papers_equations(r: float, Ci: float) -> None:
    """The templated kernel against a by-hand transcription of Eqs. 11/15,
    both branches."""
    assert cb.jobson_contraction_coefficient(r, Ci, 1.4) == pytest.approx(
        _jobson_literal(r, Ci, 1.4), rel=1e-12
    )


@pytest.mark.parametrize("r", [0.95, 0.7, 0.53, 0.5, 0.3, 0.1])
def test_the_borda_mouthpiece_is_eqs_16_and_17(r: float) -> None:
    """f = 0 (Ci = 0.5): the quadratic turns linear, and the paper gives the
    closed forms Eq. (16) and Eq. (17) (simplified: (1-r)/(rc(1+n) - r))."""
    n = 1.4
    rc = (2.0 / (n + 1.0)) ** (n / (n - 1.0))
    if r >= rc:
        K2 = 2.0 * n / (n - 1.0) * r ** (2.0 / n) * (1.0 - r ** ((n - 1.0) / n))
        C_B = r ** (1.0 / n) * (1.0 - r) / K2
    else:
        C_B = (1.0 - r) / (rc * (1.0 + n) - r)
    assert cb.jobson_contraction_coefficient(r, 0.5, n) == pytest.approx(C_B, rel=1e-10)


def test_force_defect_is_eq5_and_the_beta_extension() -> None:
    for Ci in (0.5, 0.611, 0.7):
        assert cb.jobson_force_defect(Ci) == pytest.approx(1.0 / Ci - 0.5 / Ci**2)
        assert cb.jobson_force_defect(Ci, 0.6) == pytest.approx(
            1.0 / Ci - 0.5 / Ci**2 - 0.36 + 0.5 * 0.6**4
        )
    assert cb.jobson_force_defect(0.5) == pytest.approx(0.0, abs=1e-15)  # Borda


def test_the_network_kernel_reproduces_table2_on_real_air() -> None:
    """The kernel the solver runs, beta 0, cold air (gamma ~1.400): Table 2
    still holds to its printed digits."""
    Pt = 1.0e5
    for r, C_paper in TABLE2:
        res = _solver_tools.orifice_jobson_residuals_and_jacobian(
            0.0, Pt, Pt, 288.15, Y_AIR, r * Pt, CI_SLIT, 1e-4, 0.0, 1.0
        )
        assert res.contraction == pytest.approx(C_paper, abs=1e-3), r


# ---------------------------------------------------------------------------
# Cross-source accuracy: the beta extension against ISO 5167-2 Eq. (5)
# ---------------------------------------------------------------------------


def _rhg_C(beta: float, Re_D: float = 1e6) -> float:
    g = cb.OrificeGeometry()
    g.d, g.D, g.tap = 0.1 * beta, 0.1, cb.OrificeTap.Flange
    return cb.metering_cd_and_dRe(cb.MeteringCdCorrelation.ReaderHarrisGallagher, g, Re_D)[0]


@pytest.mark.parametrize("beta", [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75])
def test_cross_source_beta_extension_is_inside_iso_eq5_uncertainty(beta: float) -> None:
    """CROSS-SOURCE: the expansibility the beta-extended theory implies
    against ISO 5167-2:2003 Eq. (5), kappa 1.4, p2/p1 0.75-0.98, at the
    Reader-Harris/Gallagher C the standard pairs it with (Re_D 1e6). Nothing
    is fitted to ISO. The band is ISO's OWN stated uncertainty of eps,
    3.5 dp/(kappa p1) % (5.3.3.2), a measurement band. Worst case on this
    grid is 0.93 of the band (beta 0.5, p2/p1 0.75)."""
    Ci = cb.jobson_ci_from_cd(_rhg_C(beta), beta, 1.0)
    for r in (0.98, 0.95, 0.9, 0.85, 0.8, 0.75):
        eps_j = cb.jobson_expansibility(r, beta, Ci, 1.4)
        eps_iso = cb.expansibility_factor(beta, 1.0 - r, 1.0, 1.4)
        band = 3.5 * (1.0 - r) / 1.4 / 100.0
        assert abs(eps_j / eps_iso - 1.0) <= band, (beta, r)


def test_the_nozzle_misses_iso_eq5_by_far_more() -> None:
    """Why the model exists: today's nozzle expansion is 4-12% below ISO
    over the same grid, i.e. 7-19 times ISO's band at p2/p1 0.75."""
    for beta in (0.3, 0.5, 0.7):
        r, n = 0.75, 1.4
        b4 = beta**4
        eps_noz = math.sqrt(
            n
            / (n - 1)
            * r ** (2 / n)
            * (1 - r ** ((n - 1) / n))
            / (1 - r)
            * (1 - b4)
            / (1 - b4 * r ** (2 / n))
        )
        eps_iso = cb.expansibility_factor(beta, 1.0 - r, 1.0, n)
        assert eps_noz / eps_iso - 1.0 < -3.0 * 3.5 * (1.0 - r) / n / 100.0


# ---------------------------------------------------------------------------
# The kernels: choke continuity and the analytic Jacobian
# ---------------------------------------------------------------------------


def _kernel(model: str, Pt: float, P1: float, Tt: float, Pb: float, Cd: float, beta: float):
    if model == "jobson1955":
        return _solver_tools.orifice_jobson_residuals_and_jacobian(
            0.0, Pt, P1, Tt, Y_AIR, Pb, Cd, 1e-4, beta, 0.98
        )
    return _solver_tools.orifice_iso5167_residuals_and_jacobian(
        0.0, Pt, P1, Tt, Y_AIR, Pb, Cd, 1e-4, beta
    )


@pytest.mark.parametrize("beta", [0.0, 0.5])
def test_jobson_flow_and_slope_are_continuous_through_choke(beta: float) -> None:
    """mdot and d(mdot)/dPb are C1 through r_c; past it the flow keeps
    rising as the jet widens (Jobson's Fig. 3), never flat like a nozzle."""
    Pt, Tt = 3.0e5, 400.0
    P1 = 0.995 * Pt if beta else Pt
    n = cb.isentropic_expansion_coefficient(Tt, X_AIR)
    rc = (2.0 / (n + 1.0)) ** (n / (n - 1.0))
    pbs = np.linspace(rc - 0.02, rc + 0.02, 81) * Pt
    res = [_kernel("jobson1955", Pt, P1, Tt, pb, 0.61, beta) for pb in pbs]
    m = np.array([r.m_dot_calc for r in res])
    d = np.array([r.d_mdot_dP_static_down for r in res])
    assert any(r.choked for r in res) and not all(r.choked for r in res)
    assert np.all(np.diff(m) < 0.0)  # more flow at every lower back pressure
    assert np.all(d < 0.0)
    # No kink: the slope varies smoothly across the switch.
    assert np.max(np.abs(np.diff(d))) <= 0.05 * np.max(np.abs(d))


@pytest.mark.parametrize("model", ["jobson1955", "iso5167"])
@pytest.mark.parametrize("beta", [0.0, 0.6])
@pytest.mark.parametrize("Pb", [1.8e5, 1.07e5, 1.05e5, 0.6e5])
def test_kernel_jacobian_matches_central_differences(model: str, beta: float, Pb: float) -> None:
    """Analytic (DualN) columns against a test-only central difference,
    subcritical, either side of choke (r_c ~ 0.528) and deep choke."""
    x = {"Pt": 2.0e5, "P1": 1.98e5 if beta else 2.0e5, "Tt": 600.0, "Pb": Pb, "Cd": 0.61}
    r = _kernel(model, beta=beta, **x)
    an = {
        "Pt": r.d_mdot_dP_total_up,
        "P1": r.d_mdot_dP_static_up,
        "Tt": r.d_mdot_dTt_up,
        "Pb": r.d_mdot_dP_static_down,
        "Cd": r.d_mdot_dCd,
    }
    for k in x:
        if k == "P1" and beta == 0.0:
            assert an[k] == 0.0  # a reservoir does not read the static
            continue
        h = 1e-6 * abs(x[k])
        xp, xm = dict(x), dict(x)
        xp[k] += h
        xm[k] -= h
        fd = (
            _kernel(model, beta=beta, **xp).m_dot_calc - _kernel(model, beta=beta, **xm).m_dot_calc
        ) / (2 * h)
        assert an[k] == pytest.approx(fd, rel=1e-6, abs=1e-9 * r.m_dot_calc / x[k]), k


def test_iso5167_kernel_flags_below_its_range() -> None:
    for r, ok in ((0.9, True), (0.75, True), (0.74, False), (0.5, False)):
        res = _kernel("iso5167", 2.0e5, 2.0e5, 300.0, r * 2.0e5, 0.6, 0.5)
        assert res.in_range is ok, r


# ---------------------------------------------------------------------------
# The network element
# ---------------------------------------------------------------------------


def _plenum_net(model: str, Pout: float, **kw: object) -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0))
    net.add_node(PressureBoundary("out", Pt=Pout, Tt=300.0))
    kw.setdefault("correlation", "ReaderHarrisGallagher")
    net.add_element(OrificeElement("x", "in", "out", diameter=0.02, compressible_model=model, **kw))
    return net


def _pipe_net(model: str, Pout: float = 1.6e5) -> FlowNetwork:
    """pipe -> orifice (beta 0.5) -> pipe: approach static and recovery."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0))
    net.add_node(PressureBoundary("out", Pt=Pout, Tt=300.0))
    from combaero.network import PlenumNode

    net.add_node(PlenumNode("a"))
    net.add_node(PlenumNode("b"))
    net.add_element(ChannelElement("c1", "in", "a", length=0.5, diameter=0.1))
    net.add_element(OrificeElement("x", "a", "b", diameter=0.05, compressible_model=model))
    net.add_element(ChannelElement("c2", "b", "out", length=0.5, diameter=0.1))
    return net


def _solve(net: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    return s, r


def test_the_nozzle_reports_nothing_new() -> None:
    """The nozzle (explicit here; the default for 'fixed' and Stolz, #861)
    reports no model diagnostics."""
    o = OrificeElement("x", "a", "b", Cd=0.6, diameter=0.01, correlation="fixed")
    assert o.compressible_model == "nozzle" and o.Cv == 1.0
    _, r = _solve(_plenum_net("nozzle", 1.6e5))
    assert "x.dP_tap" in r  # diagnostics are reported...
    for key in ("compressible_model", "expansibility", "Ci", "model_in_range"):
        assert f"x.{key}" not in r  # ...but none of the opt-in models'


def test_the_models_report_their_state_and_iso_flags_its_range() -> None:
    _, r = _solve(_plenum_net("jobson1955", 1.6e5))
    assert r["x.compressible_model"] == "jobson1955"
    assert r["x.Ci"] == pytest.approx(r["x.Cd"], rel=1e-12)  # Cv = 1, beta 0
    assert r["x.contraction_C"] > r["x.Ci"] and r["x.model_in_range"] == 1.0
    assert 0.9 < r["x.expansibility"] < 1.0
    _, r = _solve(_plenum_net("iso5167", 1.8e5))
    assert r["x.model_in_range"] == 1.0
    _, r = _solve(_plenum_net("iso5167", 1.2e5))
    assert r["x.model_in_range"] == 0.0


def test_refusals() -> None:
    with pytest.raises(ValueError, match="compressible_model"):
        OrificeElement("x", "a", "b", diameter=0.01, compressible_model="iso")
    with pytest.raises(ValueError, match="ReaderHarrisGallagher"):
        OrificeElement(
            "x", "a", "b", diameter=0.01, compressible_model="iso5167", correlation="Stolz"
        )
    with pytest.raises(ValueError, match="Cv"):
        OrificeElement("x", "a", "b", diameter=0.01, Cv=0.98)
    with pytest.raises(ValueError, match="Cv"):
        OrificeElement("x", "a", "b", diameter=0.01, compressible_model="jobson1955", Cv=1.2)


def test_incompressible_regime_ignores_the_model() -> None:
    flows = []
    for model in ("nozzle", "jobson1955", "iso5167"):
        net = _plenum_net(model, 1.8e5, regime="incompressible")
        flows.append(_solve(net)[1]["x.m_dot"])
    assert flows[1] == flows[0] and flows[2] == flows[0]


@pytest.mark.parametrize("Pout", [1.9e5, 1.6e5, 0.6e5])
def test_jobson_passes_more_flow_than_the_nozzle_and_iso_tracks_it_in_range(Pout: float) -> None:
    """A sharp plate's jet widens as the pressure ratio falls: Jobson sits
    above the nozzle (which is flat past choke), by more the lower p2/p1;
    inside ISO's range the ISO type and Jobson agree to ~1%."""
    m = {k: _solve(_plenum_net(k, Pout))[1]["x.m_dot"] for k in ("nozzle", "jobson1955", "iso5167")}
    assert m["jobson1955"] > m["nozzle"]
    if Pout / 2.0e5 >= 0.75:
        assert m["iso5167"] == pytest.approx(m["jobson1955"], rel=0.01)
    if Pout == 0.6e5:
        assert m["jobson1955"] > 1.2 * m["nozzle"]


def _row_fd_check(s: NetworkSolver, x: np.ndarray, owner: str) -> None:
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    rows = s._row_kinds(x)
    i = next(k for k, (o, _) in enumerate(rows) if f"'{owner}'" in o)
    scale = (np.abs(J[i]) * np.abs(x)).max()
    for j, name in enumerate(s.unknown_names):
        h = 1e-6 * max(abs(x[j]), 1e-3)
        e = np.zeros_like(x)
        e[j] = h
        fd = (F(x + e)[i] - F(x - e)[i]) / (2 * h)
        assert J[i, j] * abs(x[j]) == pytest.approx(fd * abs(x[j]), rel=1e-5, abs=1e-7 * scale), (
            name
        )


@pytest.mark.parametrize("model", ["jobson1955", "iso5167"])
def test_the_pipe_orifice_row_is_exact(model: str) -> None:
    """Approach static, own-flow consistency and pressure recovery: the
    row's columns against central differences, at and off the root."""
    s, r = _solve(_pipe_net(model))
    x = np.array(r["__x_solution__"])
    _row_fd_check(s, x, "x")
    _row_fd_check(s, x * (1.0 + 0.003 * np.sin(np.arange(x.size))), "x")


def test_a_reynolds_dependent_hole_uses_the_nonlinear_consistency_solve() -> None:
    """McGreehan-Schotsch reads its own flow through Re; Jobson's flow is
    not linear in Cd, so the element solves m* = m_calc(m*) on the kernel's
    own d(m)/d(Cd). The row is exact."""
    net = _plenum_net("jobson1955", 0.9e5, correlation="McGreehanSchotsch", plate_thickness=0.04)
    s, r = _solve(net)
    x = np.array(r["__x_solution__"])
    _row_fd_check(s, x, "x")
    _row_fd_check(s, x * (1.0 + 0.003 * np.cos(np.arange(x.size))), "x")


# ---------------------------------------------------------------------------
# The gap #854 was opened for
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("model", "tol"),
    [("nozzle", None), ("jobson1955", 0.002), ("iso5167", 0.0005), (None, 0.002)],
)
def test_the_iso_metering_example_reads_back_its_own_flow(
    model: str | None, tol: float | None
) -> None:
    """#852's ISO 5167-2 example (air 300 K, D 100 mm, flange taps, ~3% dp)
    reads the network's flow back through ISO's Eqs. (1), (4), (5). The
    default nozzle passes 0.8-0.9% less (the issue); the ISO type closes it
    to the regularisation and Re iteration, Jobson to <0.1% with nothing
    taken from ISO. None is the build-time default since 0.8.0 (#861):
    jobson1955 here (RHG, sharp, beta <= 0.75), so the gap is closed by
    default."""
    import functools
    import importlib.util
    from pathlib import Path

    path = (
        Path(__file__).resolve().parents[1] / "examples" / "papers" / "iso5167_orifice_metering.py"
    )
    spec = importlib.util.spec_from_file_location("iso5167_example_854", path)
    ex = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(ex)
    ex.OrificeElement = functools.partial(OrificeElement, compressible_model=model)
    errs = []
    for beta in ex.BETAS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = ex.solve_run(beta, "compressible")
        q_iso, _ = ex.iso_flow(beta, r["dP_tap"], r["P_approach"], r["T_in"], eps_on=True)
        errs.append(q_iso / r["m_dot"] - 1.0)
    if tol is None:
        assert all(0.006 < e < 0.011 for e in errs), errs  # the issue's numbers
    else:
        assert max(abs(e) for e in errs) < tol, errs
