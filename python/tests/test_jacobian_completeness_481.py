"""Jacobian completeness (#481 C1, C4, C5): the analytic network Jacobian
equals differences where it used to drop a dependence.

C1  A theta-sourced PressureLossElement read the combustor's UNBURNED
    temperature (theta = T_b/T_u - 1, and the correlation's reference T) and,
    for a head loss, the burned composition (its density) -- neither relayed:
    2.6% off in d/d(air flow), 5% in d/d(fuel flow).
C4  On the [200, 5000] K clamp a node's T is flat, but its sensitivities were
    relayed in full (5% off). The clamp is gone (#626): a node heated past
    5000 K relays its own sensitivities.
C5  The wall relay took dT/dQ = 1/(cp sum m); the mixer applies Q phi(m)
    (cb.mixer_heat_per_unit_flow: 1/|m| above MIXER_HEAT_MDOT_FLOOR, fading
    below it): 2x off below the floor and of the wrong sign for a negative
    total.
"""

from __future__ import annotations

import numpy as np
import pytest
import test_combustor_regression as tcr
from scipy.optimize._numdiff import approx_derivative

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    EnergyBoundary,
    FlowNetwork,
    LinearThetaFractionLoss,
    LinearThetaHeadLoss,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j| (see _build_residual_scales)."""
    _, J = s._residuals_and_jacobian(x)
    fd = approx_derivative(
        lambda v: s._residuals_and_jacobian(v, compute_jacobian=False)[0],
        x,
        method="3-point",
        abs_step=np.maximum(np.abs(x) * 1e-6, 1e-9),
    )
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    return float(np.max(np.abs(Js - fds) / scale))


@pytest.mark.parametrize(
    "corr",
    [LinearThetaFractionLoss(k=0.01, xi0=0.02), LinearThetaHeadLoss(k=1.0, zeta0=3.0, area=0.05)],
    ids=["fraction", "head"],
)
def test_a_theta_sourced_loss_carries_the_unburned_state(corr: object) -> None:
    """Errors were 4.4e-4 (fraction) and 6.2e-5 (head) in these scaled units."""
    net, _ = tcr._make_effective_area_inlet_network(loss_correlation=corr)
    s = NetworkSolver(net)
    r = s.solve()
    assert r["__success__"]
    assert net.elements["loss"]._theta_source_resolved == "comb"
    assert _jac_err(s, np.array(r["__x_solution__"])) < 2e-5


def test_a_node_heated_past_the_old_clamp_relays_its_temperature_sensitivity() -> None:
    """20 kW into ~1 g/s: 14 000 K. A [200, 5000] K clamp held such a plenum
    flat (its sensitivities zeroed) and discarded the energy; now it is simply
    hot, energy conserved, and its Jacobian is complete."""
    Y = list(cb.mole_to_mass(cb.species.dry_air()))
    g = FlowNetwork()
    for name, Pt in (("a", 1.1e5), ("b", 1.0e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 300.0, Y
        g.add_node(bnd)
    p = PlenumNode("p")
    g.add_node(p)
    g.add_node(PlenumNode("q"))
    for eid, a, b in (("o1", "a", "p"), ("o2", "p", "q"), ("o3", "q", "b")):
        g.add_element(OrificeElement(eid, a, b, Cd=0.8, diameter=0.01, correlation="fixed"))
    p.add_energy_boundary(EnergyBoundary(id="heat", Q=2e4))
    s = NetworkSolver(g)
    r = s.solve(auto_retry=False)
    assert r["__success__"]
    assert r["__T_clamped__"] == [] and r["__heat_faded__"] == []
    assert r["p.T"] > 1e4
    assert r["__ledger__"]["energy_closure_rel"] < 1e-12
    assert _jac_err(s, np.array(r["__x_solution__"])) < 1e-6


@pytest.mark.parametrize("m", [2.5e-7, 1e-6, 1e-3, -1e-4])
def test_the_mixer_reports_its_own_heat_sensitivity(m: float) -> None:
    """dT_mix/dQ at the heat per unit flow the mixer actually applies."""
    Y = list(cb.mole_to_mass(cb.species.dry_air()))

    def mix(Q: float) -> object:
        return _solver_tools.mixer_from_streams_and_jacobians(
            [cb.MassStream(m, 400.0, 1e5, Y)], Q=Q, fraction=0.0
        )

    h = 1e-7
    fd = (mix(1e-4 + h).T_mix - mix(1e-4 - h).T_mix) / (2 * h)
    assert mix(1e-4).dT_mix_dQ == pytest.approx(fd, rel=1e-6)
    phi, _ = cb.mixer_heat_per_unit_flow(m)
    assert mix(1e-4).dT_mix_dQ == pytest.approx(mix(1e-4).dT_mix_d_delta_h * phi, rel=1e-12)
