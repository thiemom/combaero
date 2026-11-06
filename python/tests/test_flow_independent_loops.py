"""A closed loop of elements whose total-pressure relation ignores their own
flow is refused (follow-up to #988).

#988 refused a loop of lossless connections: Pt is equal all round it, so any
circulation is a root. The same holds for every element whose Pt relation is
fixed at any flow: a fraction loss without m_ref (Pt_out = (1 - xi) Pt_in at
every flow, also at zero), a free vortex, a fixed-ratio compressor or turbine.
Round a loop made only of them the product of the Pt ratios is fixed. If it is
not 1 there is no root: a fraction loss in parallel with a lossless link
"solved" only by sitting in its zero-flow blend (stress random 5, 14, 47 after
#996, the fraction path at ~1e-20 kg/s, its second-law verdict not_checked),
or failed with a huge circulation. If it is 1 the split is not determined.
A fraction loss WITH m_ref fades to zero at zero flow, so it is determined in
parallel with a lossless link (it carries nothing), as an orifice is.
"""

from __future__ import annotations

import warnings

import pytest
from test_lossless_loop_988 import Y_AIR, _between

from combaero.network import (
    CompressorElement,
    LosslessConnectionElement,
    MassFlowBoundary,
    NetworkSolver,
    PressureLossElement,
    VortexElement,
)
from combaero.network.components import SplitterOutletElement
from combaero.network.graph import _flow_independent_dpt
from combaero.network.pressure_loss import (
    ConstantFractionLoss,
    ConstantHeadLoss,
    LinearThetaFractionLoss,
)


def _loss(eid: str, a: str, b: str, corr: object) -> PressureLossElement:
    return PressureLossElement(eid, a, b, correlation=corr)


def _refused(g, members: list[str]) -> str:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="closed loop") as exc:
            g.validate()
    msg = str(exc.value)
    named = sorted(eid for eid in g.elements if f"{eid!r}" in msg)
    assert named == sorted(members), msg
    return msg


@pytest.mark.parametrize(
    "make, why",
    [
        (lambda: _loss("x", "p1", "p2", ConstantFractionLoss(xi=0.01)), "without m_ref"),
        (
            lambda: _loss("x", "p1", "p2", LinearThetaFractionLoss(k=0.5, xi0=0.02)),
            "without m_ref",
        ),
        # with m_ref but no loss at any flow: lossless in effect
        (lambda: _loss("x", "p1", "p2", ConstantFractionLoss(xi=0.0, m_ref=0.05)), "zero loss"),
        (lambda: VortexElement("x", "p1", "p2", r_in=0.05, r_out=0.1), "free vortex"),
        (lambda: CompressorElement("x", "p1", "p2", pressure_ratio=1.2, eta_is=0.85), "ratio"),
    ],
)
def test_a_flow_independent_element_in_parallel_with_a_lossless_link_is_refused(
    make, why: str
) -> None:
    g = _between(("l", "p1", "p2"))
    g.add_element(make())
    msg = _refused(g, ["l", "x"])
    assert why in msg
    assert "independently of the flow" in msg
    assert "m_ref" in msg and "flow-dependent loss" in msg


@pytest.mark.parametrize("xi_b", [0.05, 0.01])  # unequal: no root; equal: split undetermined
def test_two_fraction_losses_without_m_ref_in_parallel_are_refused(xi_b: float) -> None:
    g = _between()
    g.add_element(_loss("fa", "p1", "p2", ConstantFractionLoss(xi=0.01)))
    g.add_element(_loss("fb", "p2", "p1", ConstantFractionLoss(xi=xi_b)))
    _refused(g, ["fa", "fb"])


def test_a_mixed_ring_is_refused_and_named() -> None:
    """lossless p1->p2, fraction loss p2->p3, free vortex p3->p1: a ring."""
    g = _between(("la", "p1", "p2"), ("o3", "p1", "p3"))
    g.add_element(_loss("fb", "p2", "p3", ConstantFractionLoss(xi=0.02)))
    g.add_element(VortexElement("vc", "p3", "p1", r_in=0.05, r_out=0.1))
    # o3 (an orifice) is in a loop too, but not in the flow-independent one
    _refused(g, ["la", "fb", "vc"])


def test_a_lossless_only_loop_keeps_the_988_message() -> None:
    g = _between(("la", "p1", "p2"), ("lb", "p2", "p1"))
    msg = _refused(g, ["la", "lb"])
    assert "the lossless connections" in msg and "Remove one of the connections" in msg


def _fraction_m_ref_parallel(Pt_in: float | None = None, m_in: float | None = None):
    g = _between(("l", "p1", "p2"))
    g.add_element(_loss("pl", "p1", "p2", ConstantFractionLoss(xi=0.01, m_ref=0.05)))
    if Pt_in is not None:
        g.nodes["in"].Pt = Pt_in
    if m_in is not None:
        b = MassFlowBoundary("in")
        b.m_dot, b.Tt, b.Y = m_in, 300.0, Y_AIR
        g.nodes["in"] = b
    return g


def _solve(g) -> dict:
    g.validate()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__outcome__")
    return r


def test_a_fraction_loss_with_m_ref_in_parallel_with_a_lossless_link_solves() -> None:
    """It fades to zero at zero flow: the lossless link takes the flow."""
    for g in (_fraction_m_ref_parallel(), _fraction_m_ref_parallel(m_in=0.1)):
        r = _solve(g)
        through = r["o_in.m_dot"]
        assert through > 1e-3
        assert abs(r["l.m_dot"] - through) < 1e-5 * through
        assert abs(r["pl.m_dot"]) < 1e-5 * through


@pytest.mark.parametrize("drive", ["pressure", "mass"])
def test_zero_drive_gives_zero_flow_with_the_m_ref_fraction_loss(drive: str) -> None:
    g = _fraction_m_ref_parallel(Pt_in=1e5) if drive == "pressure" else None
    if g is None:
        g = _fraction_m_ref_parallel(m_in=0.0)
    r = _solve(g)
    # 1e-7 kg/s is ~1e-6 of the 0.08 kg/s this network carries at 2:1 drive;
    # a residual-tolerance circulation l <-> pl of ~1e-8 kg/s remains.
    for eid in ("o_in", "l", "pl", "o_out"):
        assert abs(r[f"{eid}.m_dot"]) < 1e-7, (eid, r[f"{eid}.m_dot"])


def test_a_head_loss_in_parallel_with_a_fraction_loss_solves() -> None:
    """A head loss scales with the flow: it sets the split against a fixed
    fraction, and the loop is not refused. With zeta high enough that the
    head loss alone cannot pass the flow, the fraction loss carries the rest
    and the drop is its xi Pt_in (with a low zeta it would carry ~nothing,
    sitting in its zero-flow blend: a check valve with a cracking pressure)."""
    g = _between()
    g.add_element(_loss("f", "p1", "p2", ConstantFractionLoss(xi=0.01)))
    g.add_element(
        PressureLossElement(
            "h", "p1", "p2", correlation=ConstantHeadLoss(zeta=1000.0, area=4e-3), area=4e-3
        )
    )
    r = _solve(g)
    assert r["f.m_dot"] > 0.1 * r["o_in.m_dot"] and r["h.m_dot"] > 0.0
    dPt = r["p1.Pt"] - r["p2.Pt"]
    assert abs(dPt - 0.01 * r["p1.Pt"]) < 1e-6 * r["p1.Pt"]


def test_orifices_in_parallel_with_a_lossless_link_pass() -> None:
    g = _between(("l", "p1", "p2"), ("oa", "p1", "p2"), ("ob", "p2", "p1"))
    r = _solve(g)
    assert abs(r["oa.m_dot"]) + abs(r["ob.m_dot"]) < 1e-9 * abs(r["o_in.m_dot"])


def test_what_counts_as_flow_independent() -> None:
    def pl(corr: object) -> PressureLossElement:
        return _loss("x", "a", "b", corr)

    assert _flow_independent_dpt(LosslessConnectionElement("x", "a", "b"))
    assert _flow_independent_dpt(pl(ConstantFractionLoss(xi=0.03)))
    assert _flow_independent_dpt(pl(LinearThetaFractionLoss()))
    assert _flow_independent_dpt(pl(ConstantFractionLoss(xi=0.03, m_ref=0.1))) is None
    assert _flow_independent_dpt(pl(LinearThetaFractionLoss(m_ref=0.1))) is None
    assert _flow_independent_dpt(pl(ConstantHeadLoss())) is None
    assert _flow_independent_dpt(pl(lambda ctx: (0.01, 0.0))) is None  # unknown: not refused
    outlets = ["s_a", "s_b"]
    ref = SplitterOutletElement("s_a", "a", "b", 0.5, outlets, reference=True)
    metered = SplitterOutletElement("s_b", "a", "c", 0.5, outlets)
    assert _flow_independent_dpt(ref)
    assert _flow_independent_dpt(metered) is None  # its row fixes the flow, Pt is free
