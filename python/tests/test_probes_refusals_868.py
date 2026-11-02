"""Inner-solve follow-ups of the #868 audit: R4 (#937), R12 (#944), R13 (#945).

R4: chamber_static_state refuses a non-positive or non-finite Tt (a solver
probe's garbage upstream state) instead of returning a negative static T.
R12: dewpoint's RH > 1 behaviour, pinned as documented.
R13: the two ejector bisections flag a target outside their bracket and
clamp; a clamped exit-static result carries the clamped function's Jacobian.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero import _solver_tools as st

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _cd(f, x0: list[float], i: int, eps_rel: float = 1e-6) -> float:
    eps = eps_rel * abs(x0[i])
    xp, xm = list(x0), list(x0)
    xp[i] += eps
    xm[i] -= eps
    return (f(xp) - f(xm)) / (2.0 * eps)


# --- R4 (#937) ---------------------------------------------------------------


@pytest.mark.parametrize("Tt", [-1.0, 0.0, -1.393e19, math.nan, math.inf])
def test_chamber_static_state_refuses_a_garbage_stagnation_T(Tt: float) -> None:
    # Before: Tt = -1 K returned a static T of -1.00001 K, Tt = 0 returned 0.
    with pytest.raises(RuntimeError, match="chamber_static_state: Tt must be finite and > 0"):
        st.momentum_chamber_residual_and_jacobian(1.0e5, 1.01e5, 0.5, Tt, Y_AIR, 0.01)


def test_chamber_static_state_unchanged_at_a_physical_Tt() -> None:
    r = st.momentum_chamber_residual_and_jacobian(1.0e5, 1.01e5, 0.5, 300.0, Y_AIR, 0.01)
    assert 0.0 < r.T_static < 300.0


def test_a_probe_with_a_garbage_chamber_Tt_is_rejected_not_a_crash() -> None:
    """The cascaded branching manifold (test_tee_network, xfail: no root at
    its near-choke ratio) is where the audit found the probes: reversed
    iterates mix a chamber Tt of -4e18 K. Those evaluations now raise the
    refusal, which the solver records as rejected probes; the solve returns
    its best iterate (main: the same |F|, the probes ending at the Newton cap)."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import test_tee_network as tee

    from combaero.network import NetworkSolver

    with pytest.warns(UserWarning):
        sol = NetworkSolver(tee._cascaded_branching_manifold()).solve()
    refused = [
        e
        for e in sol["__caught_exceptions__"]
        if "chamber_static_state: Tt must be finite and > 0" in e["message"]
    ]
    assert refused
    assert all(e["phase"] == "probe" and e["type"] == "RuntimeError" for e in refused)


# --- R12 (#944) --------------------------------------------------------------


def test_dewpoint_rh_range_is_as_documented() -> None:
    assert cb.dewpoint(300.0, 101325.0, 1.0) == 300.0
    assert cb.dewpoint(300.0, 101325.0, 0.5) < 300.0
    with pytest.raises(RuntimeError):
        cb.dewpoint(300.0, 101325.0, 1.5)
    for rh in (0.0, -0.1, math.nan):
        with pytest.raises(ValueError):
            cb.dewpoint(300.0, 101325.0, rh)


# --- R13 (#945): J2, the C-D nozzle exit-static inverse ----------------------

G, R, ETA, EPS = 1.4, 287.05, 0.95, 1e-3
P0, T0, A_T, A_E = 5.0e5, 300.0, 1.0e-4, 1.5e-4


def _exit_static(m_dot: float, p0: float = P0, t0: float = T0):
    return st.ejector_cd_nozzle_exit_static_and_jacobian(m_dot, p0, t0, A_T, A_E, G, R, ETA, EPS)


def _cap() -> float:
    return st.ejector_choked_mass_flow_and_jacobian(P0, T0, A_T, G, R, ETA)[0]


def test_exit_static_in_range_is_unflagged() -> None:
    j = _exit_static(0.5 * _cap())
    assert j.out_of_range == 0
    r_crit = (2.0 / (G + 1.0)) ** (G / (G - 1.0))
    assert r_crit * P0 < j.p_py < P0
    assert j.dp_py_dm_dot < 0.0


@pytest.mark.parametrize(
    ("m_frac", "side", "k"),
    [(1.5, 1, (2.0 / (G + 1.0)) ** (G / (G - 1.0))), (0.0, -1, 1.0 - 1e-12)],
)
def test_exit_static_out_of_range_flags_and_clamps(m_frac: float, side: int, k: float) -> None:
    m = m_frac * _cap()
    j = _exit_static(m)
    assert j.out_of_range == side
    assert j.p_py == pytest.approx(k * P0, rel=1e-12)
    # The clamped function k Pt: flat in m_dot and Tt, slope k in Pt.
    assert (j.dp_py_dm_dot, j.dp_py_dt0) == (0.0, 0.0)
    assert j.dp_py_dp0 == k
    # And it is the derivative of the value (central differences).
    x0 = [m if m > 0.0 else 1e-9 * _cap(), P0, T0]

    def f(x: list[float]) -> float:
        return _exit_static(x[0], x[1], x[2]).p_py

    assert _cd(f, x0, 0) == pytest.approx(0.0, abs=1e-9 * P0 / x0[0])
    assert _cd(f, x0, 1) == pytest.approx(k, rel=1e-6)
    assert _cd(f, x0, 2) == pytest.approx(0.0, abs=1e-9 * P0 / T0)


# --- R13 (#945): J1, the nozzle-exit Mach -----------------------------------


def test_nozzle_exit_mach_flags_an_area_ratio_beyond_the_bracket() -> None:
    inside = st.ejector_nozzle_exit_mach(2.0, 1.4)
    assert inside.out_of_range == 0
    assert inside.mach == pytest.approx(2.1971981216521863, rel=1e-14)
    # A/A*(M = 50) is ~7.8e3 at gamma 5/3 and ~1.5e6 at 1.4.
    assert st.ejector_nozzle_exit_mach(1.0e4, 1.4).out_of_range == 0
    above = st.ejector_nozzle_exit_mach(1.0e4, 5.0 / 3.0)
    assert (above.out_of_range, above.mach) == (1, pytest.approx(50.0, rel=1e-12))


def _element_geo(arn: float, arm: float, gamma: float) -> tuple:
    return (arn, arm, A_T, arn * A_T, (arm - arn) * A_T, gamma, R, 0.95, 0.85, 1.0, 1e-3)


# A clamped M_p1 needs A_p1/A_t above A/A*(M = 50): 1e4 at gamma 5/3. (At
# gamma 1.4 that is 1.5e6, where row 0's smooth-min of an exit flux 1e6 times
# the cap loses ~1e-10 to cancellation and a 1e-6 FD step reads it; the same
# on main, and not the clamp.)
@pytest.mark.parametrize(("arn", "gamma", "flag"), [(1.86, 1.4, 0), (1.0e4, 5.0 / 3.0, 1)])
def test_element_reports_the_clamp_and_its_rows_match_fd(
    arn: float, gamma: float, flag: int
) -> None:
    geo = (*_element_geo(arn, 3.0 * arn, gamma), 0.95, 0.999, 0.98, 1.02)
    cap = st.ejector_choked_mass_flow_and_jacobian(P0, T0, A_T, gamma, R, 0.95)[0]
    x0 = [cap * 0.995, 0.6 * cap, 1.6 * cap, P0, T0, 1.01e5, T0, 1.05e5, 0.5e5]

    def rows(x: list[float]):
        return st.ejector_element_residuals_and_jacobian(*x, *geo, anchor_closure="huang1999")

    res = rows(x0)
    assert res.nozzle_exit_mach_out_of_range == flag
    for i in range(4):
        scale = max(abs(res.residuals[i]), max(abs(res.jacobian[i][j] * x0[j]) for j in range(9)))
        for j in range(9):
            fd = _cd(lambda x, i=i: rows(x).residuals[i], x0, j)
            assert abs(res.jacobian[i][j] - fd) * abs(x0[j]) < 1e-5 * scale, (i, j)


def test_ejector_diagnostics_report_the_nozzle_exit_mach_flag() -> None:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import test_ejector_element as te

    from combaero.network import NetworkSolver

    sol = NetworkSolver(te._build_network()).solve(timeout=60.0)
    assert sol["__success__"]
    flags = [v for k, v in sol.items() if k.endswith(".nozzle_exit_mach_out_of_range")]
    assert flags == [0.0]
    # An A_p1/A_t above A/A*(M = 50) of air (~1.5e6): M_p1 is the clamped end.
    e = te._element(nozzle_exit_area=2.0e6 * te._THROAT_AREA, mixing_area=4.0e6 * te._THROAT_AREA)
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    diag = e.diagnostics(te._states(), 1.0e5, [-0.5, -0.3, 0.8])
    assert diag["nozzle_exit_mach_out_of_range"] == 1.0
    assert (
        te._element().diagnostics(te._states(), 1.0e5, [-0.5, -0.3, 0.8])[
            "nozzle_exit_mach_out_of_range"
        ]
        == 0.0
    )
