"""A reversed element reads its upstream side ALONG THE FLOW (#554).

An orifice or an incompressible channel declared against its flow used to keep
reading its declared inlet: the wrong side's density and drive (and, for the
orifice, the wrong side's velocity of approach). The same network declared the
two ways gave different answers -- an orifice passed 0.76x the flow -- and a
compressible orifice declared against its flow did not converge at all. Now
the declaration is only a sign convention.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.components import NetworkMixtureState


def _net(declared: str, regime: str) -> FlowNetwork:
    """a (1 bar) -[channel]- p -[orifice]- b (2 bar): the flow runs b -> a."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=1.0e5, Tt=300.0))
    g.add_node(PressureBoundary("b", Pt=2.0e5, Tt=300.0))
    g.add_node(PlenumNode("p"))
    if declared == "against":
        ends_c, ends_o = ("a", "p"), ("p", "b")
    else:
        ends_c, ends_o = ("p", "a"), ("b", "p")
    g.add_element(ChannelElement("loss", *ends_c, length=2.0, diameter=0.03, regime=regime))
    g.add_element(
        OrificeElement("o", *ends_o, Cd=0.6, diameter=0.02, regime=regime, correlation="fixed")
    )
    return g


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_the_declaration_is_only_a_sign_convention(regime: str) -> None:
    along = NetworkSolver(_net("along", regime)).solve()
    against = NetworkSolver(_net("against", regime)).solve()
    assert along["__success__"], along["__message__"]
    assert against["__success__"], against["__message__"]
    assert against["o.m_dot"] < 0.0 < along["o.m_dot"]
    assert -against["o.m_dot"] == pytest.approx(along["o.m_dot"], rel=1e-8)
    assert -against["loss.m_dot"] == pytest.approx(along["loss.m_dot"], rel=1e-8)
    assert against["p.P"] == pytest.approx(along["p.P"], rel=1e-9)


_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _states(m: float) -> tuple[NetworkMixtureState, NetworkMixtureState]:
    # Declared inlet at 1.2 bar / 300 K, outlet at 1.5 bar / 450 K: reversed,
    # and the two sides' densities differ by ~13%.
    s_in = NetworkMixtureState(P=1.2e5, Pt=1.2e5, Tt=300.0, m_dot=m, Y=_AIR)
    s_out = NetworkMixtureState(P=1.5e5, Pt=1.5e5, Tt=450.0, m_dot=m, Y=_AIR)
    return s_in, s_out


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_a_reversed_orifice_reads_the_outlet_side(regime: str) -> None:
    """Its flow is the outlet-side kernel's, negated; the Jacobian matches
    central differences in the outlet node's columns."""
    o = OrificeElement("o", "A", "B", Cd=0.6, diameter=0.02, regime=regime, correlation="fixed")
    s_in, s_out = _states(-0.1)
    m_calc, d = o._flow_and_partials(s_in, s_out, 0.6)
    m_ref, _ = o._kernel_from(s_out, s_in, "B", "A", 0.6)
    assert m_calc == pytest.approx(-m_ref, rel=1e-12)
    assert "B.Tt" in d and "A.Tt" not in d  # the outlet side's temperature drives it

    def m_at(**kw: float) -> float:
        a, b = _states(-0.1)
        for k, v in kw.items():
            node, attr = k.split("_")
            st = a if node == "A" else b
            setattr(st, attr, v)
            if attr == "Pt":
                st.P = v
        return o._flow_and_partials(a, b, 0.6)[0]

    for col, base, h in (("B_Pt", 1.5e5, 1.0), ("A_Pt", 1.2e5, 1.0), ("B_Tt", 450.0, 1e-3)):
        fd = (m_at(**{col: base + h}) - m_at(**{col: base - h})) / (2 * h)
        name = col.replace("_", ".")
        ana = d.get(name, 0.0) + (d.get(name.replace(".Pt", ".P"), 0.0) if ".Pt" in name else 0.0)
        assert ana == pytest.approx(fd, rel=1e-5), col


def test_a_reversed_incompressible_channel_reads_the_outlet_side() -> None:
    """The drop of a channel declared A -> B with flow B -> A is the drop the
    same channel declared B -> A gives, at the outlet side's density."""
    rev = ChannelElement("c", "A", "B", length=2.0, diameter=0.03, regime="incompressible")
    fwd = ChannelElement("c", "B", "A", length=2.0, diameter=0.03, regime="incompressible")
    s_in, s_out = _states(-0.1)
    r_rev, j_rev = rev.residuals(s_in, s_out)
    f_in, f_out = _states(0.1)
    f_in, f_out = f_out, f_in  # B is the declared inlet now
    r_fwd, _ = fwd.residuals(f_in, f_out)
    # Pt_from - Pt_to - dP: the same equation with every sign flipped.
    assert r_rev[0] == pytest.approx(-r_fwd[0], rel=1e-12)
    assert "B.Tt" in j_rev[0] and "A.Tt" not in j_rev[0]


def test_the_cold_start_stays_between_the_boundaries() -> None:
    """Driven by pressure boundaries alone, a passive node's seed lies within
    their range: wired 1 bar -> p -> 2 bar it was guessed at 0.5 bar."""
    s = NetworkSolver(_net("against", "compressible"))
    values = s._build_x0()
    x0 = dict(zip(s.unknown_names, values, strict=True))
    assert 1.0e5 <= x0["p.P"] <= 2.0e5
    assert math.isfinite(x0["p.P"])
