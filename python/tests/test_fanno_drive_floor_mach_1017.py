"""Fanno kernels below the 1 Pa drive floor report the flow's Mach (#1017).

Below the floor the flux is blended to h(t) G0, but M_in and M_exit were the
floor's: 23x the flow's at 0.035 Pa on 20.8 bar, 80x at 0.01 Pa on 1 bar.
The second-law heat band squares that Mach, and refused the converged root of
stress random 99 (a heated channel at a 0.035 Pa drive).
"""

import math
import sys
import warnings
from pathlib import Path

import pytest

import combaero as cb
from combaero.network import NetworkSolver

X = list(cb.species.dry_air())


def _mach_from_flux(G: float, Pt0: float, Tt0: float) -> float:
    R = cb.specific_gas_constant(X)
    g = cb.isentropic_expansion_coefficient(Tt0, X)
    return G * math.sqrt(R * Tt0) / (Pt0 * math.sqrt(g))


@pytest.mark.parametrize("heated", [False, True])
@pytest.mark.parametrize(("Pt0", "drive"), [(20.8e5, 0.035), (20.8e5, 0.35), (1e5, 0.01)])
def test_kernel_mach_follows_the_blended_flux(heated: bool, Pt0: float, drive: float) -> None:
    if heated:
        f = cb.fanno_rayleigh_channel_flow(
            Pt0, 300.0, 400.0, X, Pt0 - drive, False, 0.3, 0.1, 0.0, "haaland"
        )
        at = cb.fanno_rayleigh_channel_flow(
            Pt0, 300.0, 400.0, X, Pt0 - 1.0, False, 0.3, 0.1, 0.0, "haaland"
        )
    else:
        f = cb.fanno_channel_flow(Pt0, 300.0, X, Pt0 - drive, False, 0.3, 0.1, 0.0, "haaland")
        at = cb.fanno_channel_flow(Pt0, 300.0, X, Pt0 - 1.0, False, 0.3, 0.1, 0.0, "haaland")
    assert f.M_in == pytest.approx(_mach_from_flux(f.G, Pt0, 300.0), rel=0.01)
    assert f.M_exit / f.M_in == pytest.approx(at.M_exit / at.M_in, rel=1e-9)
    assert f.M_exit_root == at.M_exit


def test_random_99_heated_channel_passes_the_second_law() -> None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
    from networks import FAMILIES

    net, _ = FAMILIES["random"](99)
    e10 = net.elements["e10"]
    s = NetworkSolver(net)
    s.homotopy_retry = True
    s.bc_ramp_retry = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=20)
    assert r["__converged__"]
    # The drive is below the floor: the case this pins.
    assert 0.0 < r["e10.dPt"] < cb.FANNO_FLOW_DRIVE_FLOOR
    # Inlet Mach from the flux and the inlet node's own state, rho a.
    G = r["e10.m_dot"] / (math.pi * e10.diameter**2 / 4.0)
    src = e10.from_node
    M_true = G / (r[f"{src}.rho"] * r[f"{src}.a"])
    assert r["e10.mach_in"] == pytest.approx(M_true, rel=0.01)
    assert r["e10.M_in_duct"] == pytest.approx(M_true, rel=0.01)
    assert r["e10.rho_ref"] == pytest.approx(r[f"{src}.rho"], rel=1e-3)
    assert r["__second_law__"]["e10"]["verdict"] == "pass"
    assert r["__success__"], r["__message__"]
