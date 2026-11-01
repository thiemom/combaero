"""An in-pipe orifice's consistent flow follows its drive (#875).

With a velocity of approach the orifice row solves m* = m_calc(m*, x), and
m_calc reads |m| through the approach pipe's static pressure. An outer
iterate above the flux at which that static falls to the outlet pressure
gives m_calc < 0; the inner solve then followed the mirrored branch down to
the approach pipe's choke and stopped on a spurious reversed root tens of
times the true flow, where the row has no pressure column at all. A fixed Cd
0.98 at beta 0.95 (the liner-cooling example's exits before #874) hit it and
the coupled solve stopped making progress.

On the mirrored side the slope of m_calc at |a| = m* is +x^2 with
x = C E beta^2 (the jet's effective area over the pipe's): it is the
fixed-point step's growth factor there, so the spurious root is only reached
for x > 1. Measured on the hot-exit network below: none at x = 0.999 from
iterates up to 1000x the root; from 27x at x = 1.001, 9.7x at x = 1.01, 3.2x
at x = 1.10 and 1.28x at x = 2.05 (Cd 0.98, beta 0.95).
"""

import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.components import NetworkMixtureState

T_HOT = 1758.0
P_HOT = 20.0 * 101325.0 * 0.98
D_PIPE = 0.30
DPT = 1000.0
Y_AIR = cb.mole_to_mass(cb.species.dry_air())


def _x(beta: float, cd: float) -> float:
    """C E beta^2: the jet's effective area over the approach pipe's."""
    return cd * beta**2 / math.sqrt(1.0 - beta**4)


def _hot_exit(cd: float, beta: float, regime: str = "incompressible") -> FlowNetwork:
    """The liner example's hot side, without the wall: a 0.3 m annulus at
    25 m/s discharging through a fixed-Cd orifice in its pipe."""
    rho = cb.density(T_HOT, P_HOT, cb.species.dry_air())
    m_dot = rho * 25.0 * math.pi * D_PIPE**2 / 4.0
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("hot_inlet", m_dot=m_dot, Tt=T_HOT, Y=Y_AIR))
    net.add_node(PlenumNode("hot_exit"))
    net.add_node(PressureBoundary("hot_outlet", Pt=P_HOT))
    net.add_element(
        ChannelElement(
            id="hot_channel",
            from_node="hot_inlet",
            to_node="hot_exit",
            length=0.35,
            diameter=D_PIPE,
            roughness=0.0,
        )
    )
    net.add_element(
        OrificeElement(
            id="orf",
            from_node="hot_exit",
            to_node="hot_outlet",
            Cd=cd,
            correlation="fixed",
            in_pipe=True,
            diameter=D_PIPE * beta,
            regime=regime,
        )
    )
    return net


def _orifice(cd: float, beta: float, regime: str = "incompressible") -> OrificeElement:
    net = _hot_exit(cd, beta, regime)
    net.resolve_all_topology()
    return net.elements["orf"]


def _row(e: OrificeElement, m_it: float, dpt: float = DPT) -> tuple[float, float]:
    """(m*, d(row)/d(upstream Pt)) at the outer iterate m_it."""
    up = NetworkMixtureState(Tt=T_HOT, P=P_HOT + dpt, Pt=P_HOT + dpt, Y=Y_AIR, m_dot=m_it)
    down = NetworkMixtureState(Tt=T_HOT, P=P_HOT, Pt=P_HOT, Y=Y_AIR, m_dot=0.0)
    r, j = e.residuals(up, down)
    return m_it - r[0], j[0].get("hot_exit.Pt", 0.0) + j[0].get("hot_exit.P", 0.0)


def test_reproducer_cd_098_beta_095_in_pipe_converges() -> None:
    """The #875 case: it stopped with 'not making good progress'."""
    net = _hot_exit(0.98, 0.95)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sol = NetworkSolver(net).solve()
    assert sol["__success__"]
    assert not [w for w in caught if "did not converge" in str(w.message)]
    assert sol["orf.m_dot"] == pytest.approx(net.nodes["hot_inlet"].m_dot, rel=1e-9)


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize("ratio", [1.0, 1.2, 1.5, 2.0, 3.0])
def test_inner_root_has_the_sign_of_the_drive(regime: str, ratio: float) -> None:
    """From any iterate the row's m* is the forward root, and its pressure
    column is the true one (FD), not the choked branch's zero."""
    e = _orifice(0.98, 0.95, regime)
    root, _ = _row(e, 1e-3)
    assert root > 0.0
    m_star, d_an = _row(e, ratio * root)
    assert m_star == pytest.approx(root, rel=1e-9)
    h = 1e-2
    r0 = ratio * root - _row(e, ratio * root)[0]
    r1 = ratio * root - _row(e, ratio * root, DPT + h)[0]
    assert d_an == pytest.approx((r1 - r0) / h, rel=1e-4)
    assert d_an < -1e-3


def test_liner_example_exits_in_pipe_converge(monkeypatch: pytest.MonkeyPatch) -> None:
    """The liner-cooling example with its two exits back in their pipes
    (in_pipe=True, as before #874's workaround) converges."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
    import combustor_liner_cooling_network as liner
    import test_combustor_liner_cooling_coupled as coupled

    def in_pipe_orifice(*args: object, **kwargs: object) -> OrificeElement:
        kwargs["in_pipe"] = True
        return OrificeElement(*args, **kwargs)

    monkeypatch.setattr(liner, "OrificeElement", in_pipe_orifice)
    coupled.test_coupled_combustor_network_converges_without_warning()


def _cd_at(beta: float, x: float) -> float:
    """The Cd that puts C E beta^2 at x."""
    return x * math.sqrt(1.0 - beta**4) / beta**2


def _unguarded_wrong_sign(cd: float, beta: float) -> bool:
    """Does the inner solve WITHOUT the #875 guard land on a reversed root
    from an iterate between 1x and 100x the root?"""
    e = _orifice(cd, beta)
    root, _ = _row(e, 1e-3)
    return any(_row(e, q * root)[0] < 0.0 for q in np.geomspace(1.0, 100.0, 240))


@pytest.mark.parametrize(
    ("beta", "x", "spurious"),
    [
        (0.85, 0.99, False),
        (0.85, 0.999, False),
        (0.85, 1.001, True),  # from 27x the root
        (0.85, 1.01, True),  # from 9.7x
        (0.95, 0.99, False),
        (0.95, 1.01, True),
        (0.95, 2.05, True),  # Cd 0.98, the #875 exits: from 1.28x
    ],
)
def test_spurious_root_boundary_is_jet_area_equal_pipe_area(
    monkeypatch: pytest.MonkeyPatch, beta: float, x: float, spurious: bool
) -> None:
    """Where the guard matters: the reversed root is reachable exactly for
    x = C E beta^2 above 1, as the mirrored slope +x^2 says."""
    monkeypatch.setattr(
        OrificeElement, "_drive_signed_root", staticmethod(lambda evaluate, a, dpt: a)
    )
    assert _unguarded_wrong_sign(_cd_at(beta, x), beta) == spurious


@pytest.mark.parametrize(("beta", "cd"), [(0.80, 1.00), (0.85, 0.90), (0.90, 0.60)])
def test_guard_is_bit_identical_below_the_bound(
    monkeypatch: pytest.MonkeyPatch, beta: float, cd: float
) -> None:
    """For x < 1 the guard never moves a root: the row is bitwise the
    unguarded one from every iterate."""
    e = _orifice(cd, beta)
    root, _ = _row(e, 1e-3)
    iterates = [q * root for q in (0.5, 1.0, 1.5, 2.5, 3.9)]
    guarded = [_row(e, m) for m in iterates]
    monkeypatch.setattr(
        OrificeElement, "_drive_signed_root", staticmethod(lambda evaluate, a, dpt: a)
    )
    assert [_row(e, m) for m in iterates] == guarded
