"""Flow-area inference from topology (issue #262).

Ground truth for these tests is the network's own geometry, not the solver's
previous output: when a component sits next to something whose area is known,
the area it resolves MUST be that area. The bug being pinned here is that
inference searched only for a neighbouring ChannelElement and, failing that,
substituted a hard-coded constant -- so a chamber, combustor or area-change
neighbour produced a fabricated area, a fabricated pressure change, and a
solver stall whose message blamed the solver.

Two behaviours are separated deliberately:

* Where the area IS the geometry being modelled (AreaChangeElement F0/F1,
  ChannelElement diameter, TeeJunctionElement F_C, BorderCarnotLossElement),
  an unresolvable area raises. A default there invents a contraction or
  expansion that does not exist.
* Where the area is only a dynamic-head reference (MomentumChamberNode,
  CombustorNode, PressureLossElement), a nominal default stands, because a
  network can legitimately never read it. It warns only when a correlation
  will actually consume it.

Both must clear the moment the user supplies the value explicitly.
"""

from __future__ import annotations

import contextlib
import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    AreaChangeElement,
    BorderCarnotLossElement,
    ChannelElement,
    CombustorNode,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    TeeJunctionElement,
)
from combaero.network.pressure_loss import ConstantFractionLoss, ConstantHeadLoss
from combaero.network.solver import NetworkSolver

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_CHAMBER_AREA = 0.15
_CHANNEL_D = 0.9
_CHANNEL_AREA = math.pi * (_CHANNEL_D / 2.0) ** 2


def _mfb(net: FlowNetwork, node_id: str = "mfb", m_dot: float = 1.0) -> None:
    net.add_node(MassFlowBoundary(node_id, m_dot=m_dot, Tt=300.0, Y=_Y))


def _pb(net: FlowNetwork, node_id: str = "pb", Pt: float = 1.0e5) -> None:
    net.add_node(PressureBoundary(node_id, Pt=Pt, Tt=300.0, Y=_Y))


# ---------------------------------------------------------------------------
# The reported repro: chamber -> area_change -> channel -> boundary
# ---------------------------------------------------------------------------


def _build_repro(F0: float | None = None) -> tuple[FlowNetwork, AreaChangeElement]:
    """The #262 network: an area change fed by a momentum chamber, not a channel.

    Rebuilt programmatically rather than loaded from the reported export --
    gui/tmp/ is gitignored, so the fixture has to carry the topology itself.
    """
    net = FlowNetwork()
    _mfb(net)
    net.add_node(MomentumChamberNode("chamber", area=_CHAMBER_AREA))
    net.add_node(PlenumNode("mid"))
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "chamber"))
    area_change = AreaChangeElement("expansion", "chamber", "mid", F0=F0)
    net.add_element(area_change)
    net.add_element(ChannelElement("duct", "mid", "pb", diameter=_CHANNEL_D, length=1.0))
    return net, area_change


def test_area_change_takes_upstream_chamber_area_not_a_default():
    """F0 must be the chamber's real area, not the 0.01 m^2 placeholder.

    0.01 against the channel's 0.636 m^2 is a 63.6x sudden expansion that
    exists nowhere in the network; its ~135 Pa dynamic head at the fictitious
    throat is what the solver used to stall on.
    """
    net, area_change = _build_repro()
    net.resolve_all_topology()

    assert pytest.approx(_CHAMBER_AREA) == area_change.F0
    assert pytest.approx(_CHANNEL_AREA) == area_change.F1


def test_fabricated_area_would_corrupt_the_chamber_pressure():
    """The consequence the fabricated area had on the answer, pinned by physics.

    The reported symptom was a solver stall, but that needed the full ejector
    network; in a plain chain the fabricated throat is worse than a stall --
    it CONVERGES, confidently, to a badly wrong pressure. Measured against the
    real 0.15 m^2 upstream area, the 0.01 m^2 placeholder inflated the chamber
    stagnation pressure by 4% at 1 kg/s, 64% at 5 kg/s and 342% at 20 kg/s.

    Ground truth here is the physics, not the previous output: expanding
    0.15 -> 0.636 m^2 at 5 kg/s carries a chamber dynamic head of only ~480 Pa,
    so the total-pressure rise above the outlet boundary cannot be more than a
    few kPa. The fabricated throat produced 64 kPa.
    """
    net, _ = _build_repro()
    net.nodes["mfb"].m_dot = 5.0

    sol = NetworkSolver(net).solve(timeout=60.0)

    assert sol["__success__"], f"solve failed: |F|={sol['__final_norm__']:.3e}"
    rise = sol["chamber.Pt"] - 1.0e5
    assert rise < 5.0e3, (
        f"chamber Pt sits {rise:.0f} Pa above the outlet -- far beyond the "
        f"~480 Pa dynamic head of a real 0.15 -> 0.636 m^2 expansion, which "
        f"means a fabricated upstream area is back."
    )


def test_explicit_F0_is_never_overridden():
    """A user-supplied area wins over anything inference would have found."""
    net, area_change = _build_repro(F0=0.05)
    net.resolve_all_topology()

    assert pytest.approx(0.05) == area_change.F0


# ---------------------------------------------------------------------------
# Unresolvable areas raise, and manual input clears the error
# ---------------------------------------------------------------------------


def _build_unresolvable(F0: float | None = None, F1: float | None = None) -> FlowNetwork:
    """An area change whose neighbours carry no area at all.

    Plenums have no flow area and lossless connections have no geometry, so
    there is genuinely nothing to infer from.
    """
    net = FlowNetwork()
    _mfb(net)
    net.add_node(PlenumNode("up"))
    net.add_node(PlenumNode("down"))
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "up"))
    net.add_element(AreaChangeElement("expansion", "up", "down", F0=F0, F1=F1))
    net.add_element(LosslessConnectionElement("exit", "down", "pb"))
    return net


def test_unresolvable_area_raises_instead_of_fabricating_one():
    net = _build_unresolvable()

    with pytest.raises(ValueError) as excinfo:
        net.resolve_all_topology()

    message = str(excinfo.value)
    assert "AreaChangeElement 'expansion'" in message, "the error must name the component"
    assert "F0" in message, "the error must name the parameter that clears it"
    assert "'up'" in message, "the error must say where it looked"


def test_manual_areas_clear_the_error():
    """The user's fix path: supply the areas, and the network resolves."""
    net = _build_unresolvable(F0=0.05, F1=0.08)

    net.resolve_all_topology()  # must not raise

    element = net.elements["expansion"]
    assert pytest.approx(0.05) == element.F0
    assert pytest.approx(0.08) == element.F1


def test_manual_area_on_the_neighbour_also_clears_the_error():
    """Giving a NEIGHBOUR a known area resolves it too -- inference, not just override."""
    net = FlowNetwork()
    _mfb(net)
    net.add_node(MomentumChamberNode("up", area=0.2))
    net.add_node(MomentumChamberNode("down", area=0.4))
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "up"))
    net.add_element(AreaChangeElement("expansion", "up", "down"))
    net.add_element(LosslessConnectionElement("exit", "down", "pb"))

    net.resolve_all_topology()

    element = net.elements["expansion"]
    assert pytest.approx(0.2) == element.F0
    assert pytest.approx(0.4) == element.F1


# ---------------------------------------------------------------------------
# The same inference gap at the sibling sites named in the issue
# ---------------------------------------------------------------------------


def test_channel_inherits_diameter_from_a_chamber():
    net = FlowNetwork()
    _mfb(net)
    net.add_node(MomentumChamberNode("chamber", area=_CHAMBER_AREA))
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "chamber"))
    channel = ChannelElement("duct", "chamber", "pb", length=1.0)
    net.add_element(channel)

    net.resolve_all_topology()

    assert math.pi * (channel.diameter / 2.0) ** 2 == pytest.approx(_CHAMBER_AREA)


def test_tee_junction_inherits_F_C_from_a_chamber():
    net = FlowNetwork()
    _mfb(net)
    net.add_node(MomentumChamberNode("common", area=_CHAMBER_AREA))
    net.add_node(PlenumNode("straight"))
    net.add_node(PlenumNode("branch"))
    _pb(net, "pb_straight")
    _pb(net, "pb_branch")
    net.add_element(LosslessConnectionElement("feed", "mfb", "common"))
    tee = TeeJunctionElement(
        "tee", common_node="common", straight_node="straight", branch_node="branch", theta=90.0
    )
    net.add_element(tee)
    net.add_element(LosslessConnectionElement("out_s", "straight", "pb_straight"))
    net.add_element(LosslessConnectionElement("out_b", "branch", "pb_branch"))

    net.resolve_all_topology()

    assert pytest.approx(_CHAMBER_AREA) == tee.F_C


def test_border_carnot_loss_inherits_area_from_a_chamber():
    net = FlowNetwork()
    _mfb(net)
    net.add_node(MomentumChamberNode("chamber", area=_CHAMBER_AREA))
    net.add_node(PlenumNode("mid"))
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "chamber"))
    loss = BorderCarnotLossElement("turn", "chamber", "mid", delta_geom_deg=90.0)
    net.add_element(loss)
    net.add_element(LosslessConnectionElement("exit", "mid", "pb"))

    net.resolve_all_topology()

    assert loss.area == pytest.approx(_CHAMBER_AREA)


def test_chamber_node_inherits_area_from_a_combustor():
    net = FlowNetwork()
    _mfb(net)
    net.add_node(CombustorNode("burner", area=_CHAMBER_AREA))
    chamber = MomentumChamberNode("chamber")
    net.add_node(chamber)
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "burner"))
    net.add_element(LosslessConnectionElement("mid", "burner", "chamber"))
    net.add_element(LosslessConnectionElement("exit", "chamber", "pb"))

    net.resolve_all_topology()

    assert chamber.area == pytest.approx(_CHAMBER_AREA)
    assert chamber._area_source != "default"


def test_pressure_loss_inherits_area_from_a_chamber():
    net = FlowNetwork()
    _mfb(net)
    net.add_node(MomentumChamberNode("chamber", area=_CHAMBER_AREA))
    net.add_node(PlenumNode("mid"))
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "chamber"))
    loss = PressureLossElement("loss", "chamber", "mid", correlation=ConstantHeadLoss(zeta=5.0))
    net.add_element(loss)
    net.add_element(LosslessConnectionElement("exit", "mid", "pb"))

    net.resolve_all_topology()

    assert loss.area == pytest.approx(_CHAMBER_AREA)
    assert loss.correlation.area == pytest.approx(_CHAMBER_AREA), "correlation must see it too"


# ---------------------------------------------------------------------------
# Guards on the inference itself
# ---------------------------------------------------------------------------


def test_channel_wins_over_a_chamber():
    """Channels keep priority so networks that resolved before resolve the same.

    A chamber and a channel both border the area change; the channel's area is
    the one the pre-#262 code would have picked, and it must stay that way.
    """
    net = FlowNetwork()
    _mfb(net)
    net.add_node(MomentumChamberNode("chamber", area=_CHAMBER_AREA))
    net.add_node(PlenumNode("mid"))
    _pb(net)
    net.add_element(ChannelElement("duct", "mfb", "chamber", diameter=_CHANNEL_D, length=1.0))
    area_change = AreaChangeElement("expansion", "chamber", "mid")
    net.add_element(area_change)
    # Give the downstream side a resolvable area of its own, so this test
    # isolates the upstream (channel-vs-chamber) priority question.
    net.add_element(ChannelElement("exit", "mid", "pb", diameter=0.4, length=1.0))

    net.resolve_all_topology()

    assert pytest.approx(_CHANNEL_AREA) == area_change.F0


def test_auto_sized_neighbour_is_not_used_as_a_source():
    """An auto-sized area is a placeholder; inheriting it launders a guess.

    The chamber here never received an area, so the area change must raise
    rather than adopt the chamber's own nominal default.
    """
    net = FlowNetwork()
    _mfb(net)
    net.add_node(MomentumChamberNode("chamber"))
    net.add_node(PlenumNode("mid"))
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "chamber"))
    net.add_element(AreaChangeElement("expansion", "chamber", "mid"))
    net.add_element(LosslessConnectionElement("exit", "mid", "pb"))

    with pytest.raises(ValueError, match="cannot determine"):
        net.resolve_all_topology()


# ---------------------------------------------------------------------------
# The dynamic-head tier: default allowed, but never silent when it is read
# ---------------------------------------------------------------------------


def _build_defaulting_loss(correlation: object) -> FlowNetwork:
    net = FlowNetwork()
    _mfb(net)
    net.add_node(PlenumNode("up"))
    net.add_node(PlenumNode("down"))
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "up"))
    net.add_element(PressureLossElement("loss", "up", "down", correlation=correlation))
    net.add_element(LosslessConnectionElement("exit", "down", "pb"))
    return net


def test_defaulted_area_warns_when_a_head_loss_will_consume_it():
    """A head loss scales as 1/area^2, so a nominal area silently rescales it."""
    net = _build_defaulting_loss(ConstantHeadLoss(zeta=5.0))

    with pytest.warns(UserWarning, match="no flow area could be inferred"):
        net.resolve_all_topology()

    assert net.elements["loss"]._area_source == "default"


@contextlib.contextmanager
def _no_area_warning():
    """Assert that no flow-area warning is emitted inside the block."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        yield
    offending = [
        str(w.message) for w in caught if "no flow area could be inferred" in str(w.message)
    ]
    assert not offending, f"unexpected area warning: {offending}"


def test_defaulted_area_is_silent_when_nothing_reads_it():
    """A fraction-based loss never touches the area -- warning there is noise."""
    net = _build_defaulting_loss(ConstantFractionLoss(xi=0.02))

    with _no_area_warning():
        net.resolve_all_topology()


def test_explicit_area_suppresses_the_warning():
    """Manual input clears the warning as well as the error."""
    net = FlowNetwork()
    _mfb(net)
    net.add_node(PlenumNode("up"))
    net.add_node(PlenumNode("down"))
    _pb(net)
    net.add_element(LosslessConnectionElement("feed", "mfb", "up"))
    net.add_element(
        PressureLossElement("loss", "up", "down", correlation=ConstantHeadLoss(zeta=5.0), area=0.25)
    )
    net.add_element(LosslessConnectionElement("exit", "down", "pb"))

    with _no_area_warning():
        net.resolve_all_topology()

    assert net.elements["loss"].area == pytest.approx(0.25)
