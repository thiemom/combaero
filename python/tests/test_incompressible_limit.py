"""The incompressible regime is the low-Mach limit of the compressible one (#588).

For each element with both closures, the same network is solved in both
regimes at the same flow, and the total-pressure drops compared. The gap must
vanish as Ma^2: halving the Mach number quarters it. A gap that stalls at a
finite value means the two regimes model different physics -- the channel's
exit head at a pressure-boundary sink was exactly that (an O(1) gap, the
incompressible drop a third of the compressible one at any Mach), until the
incompressible form honoured the same sink coupling (#360).
"""

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)

P_OUT = 101325.0
T_IN = 300.0
D = 0.05
AREA = math.pi / 4 * D * D
CD = 0.65


def _dpt(kind: str, regime: str, m_dot: float, sink: str = "boundary") -> float:
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("inlet", m_dot=m_dot, Tt=T_IN))
    net.add_node(PressureBoundary("outlet", Pt=P_OUT, Tt=T_IN))
    to = "outlet"
    if sink == "chamber":
        # A momentum chamber between the duct and the sink takes the duct's
        # stagnation pressure; the head is then not lost at the duct exit.
        net.add_node(MomentumChamberNode("mcn", area=AREA))
        to = "mcn"
        net.add_element(
            OrificeElement(
                "exit", "mcn", "outlet", Cd=1.0, diameter=D, correlation="fixed", in_pipe=False
            )
        )
    if kind == "channel":
        el = ChannelElement("el", "inlet", to, length=1.0, diameter=D, regime=regime)
    else:
        el = OrificeElement(
            "el", "inlet", to, Cd=CD, diameter=D, regime=regime, correlation="fixed"
        )
    net.add_element(el)
    sol = NetworkSolver(net).solve()
    assert sol["__success__"], (kind, regime, m_dot, sol["__message__"])
    return float(sol["inlet.Pt"]) - float(sol[f"{to}.Pt"])


def _m_for_mach(kind: str, mach: float) -> float:
    rho = float(cb.density(T_IN, P_OUT, cb.species.dry_air()))
    a = float(cb.speed_of_sound(T_IN, cb.species.dry_air()))
    # Duct Mach for a channel, vena-contracta Mach for an orifice.
    return mach * a * rho * (AREA if kind == "channel" else CD * AREA)


def _gaps(kind: str, machs: tuple[float, ...], sink: str = "boundary") -> list[float]:
    out = []
    for ma in machs:
        m = _m_for_mach(kind, ma)
        dc = _dpt(kind, "compressible", m, sink)
        out.append((_dpt(kind, "incompressible", m, sink) - dc) / dc)
    return out


# Mach ladders end where the gap is still well above the solver's convergence
# tolerance: an orifice at vena-contracta Ma 0.025 drops 44 Pa and its gap
# there, 0.01 Pa, is solver noise rather than physics.
@pytest.mark.parametrize(
    ("kind", "machs"),
    [("channel", (0.1, 0.05, 0.025)), ("orifice", (0.2, 0.1, 0.05))],
)
def test_incompressible_is_the_low_mach_limit(kind: str, machs: tuple[float, ...]) -> None:
    gaps = _gaps(kind, machs)
    # The gap closes: under 0.2% at the bottom of the ladder.
    assert abs(gaps[-1]) < 2e-3, gaps
    # ... and closes as Ma^2: each halving divides it by ~4 (3.4 to 4.6 allows
    # for the Re dependence of the friction factor, which drifts with flow).
    for g_hi, g_lo in zip(gaps, gaps[1:], strict=False):
        assert 3.4 < g_hi / g_lo < 4.6, gaps


def test_channel_exit_head_lost_only_at_a_static_sink() -> None:
    """Into a pressure boundary the incompressible duct loses its exit head;
    into a momentum chamber it does not -- the chamber carries it as its
    dynamic head -- as for the compressible form."""
    m = _m_for_mach("channel", 0.1)
    rho = float(cb.density(T_IN, P_OUT, cb.species.dry_air()))
    q = 0.5 * m * m / (rho * AREA * AREA)
    to_sink = _dpt("channel", "incompressible", m, "boundary")
    to_chamber = _dpt("channel", "incompressible", m, "chamber")
    # The difference is the exit head at the inlet-static density the
    # incompressible form uses, to within that density's few-per-mille spread.
    assert to_sink - to_chamber == pytest.approx(q, rel=2e-2)
    # Both regimes agree on the chamber case too.
    gaps = _gaps("channel", (0.05,), "chamber")
    assert abs(gaps[0]) < 5e-3, gaps
