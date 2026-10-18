"""A wall-heated compressible channel flags the heating its momentum omits (#551).

The compressible kernel is adiabatic Fanno flow at the inlet Tt; a wall's heat
only raises the downstream node's Tt. To first order Rayleigh flow loses
dPt / Pt = (gamma M^2 / 2) ln(Tt_exit / Tt_in) (cooling gains it), which the
diagnostics now report and the Inspector flags above 0.1% of the inlet Pt.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import pytest

from combaero.network import ChannelElement, FlowNetwork, NetworkSolver, PressureBoundary

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402


def _solve(g: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve(timeout=30)
    assert r["__success__"], r["__message__"]
    return r


@pytest.mark.parametrize("seed", [10, 29])
def test_a_heated_and_a_cooled_duct_report_their_omitted_rayleigh_change(seed: int) -> None:
    """The stress hotcold pairs that failed the second law: the heated side
    omits a loss, the cooled side a gain, both flagged; the estimate is the
    first-order formula on the duct's own Mach numbers and end Tt."""
    g, _ = FAMILIES["hotcold"](seed)
    r = _solve(g)
    d = r["__element_diag__"]
    assert d["cch"]["rayleigh_omitted_frac"] > 1e-2
    assert d["hch"]["rayleigh_omitted_frac"] < -1e-2
    assert d["cch"]["rayleigh_omitted"] == 1.0 and d["hch"]["rayleigh_omitted"] == 1.0
    # The chamber carries the stagnation T it mixes (#576).
    Tt_up, Tt_dn = g.nodes["ci"].Tt, r["cm.T"]
    m2 = 0.5 * (d["cch"]["M_in_duct"] ** 2 + d["cch"]["M_exit_duct"] ** 2)
    expected = 0.5 * d["cch"]["gamma_ref"] * m2 * math.log(Tt_dn / Tt_up)
    assert d["cch"]["rayleigh_omitted_frac"] == pytest.approx(expected, rel=1e-9)


def test_a_duct_without_a_wall_reports_nothing() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2e5, Tt=300.0))
    g.add_node(PressureBoundary("b", Pt=1e5, Tt=300.0))
    g.add_element(ChannelElement("c", "a", "b", length=1.0, diameter=0.01, regime="compressible"))
    d = _solve(g)["__element_diag__"]["c"]
    assert "rayleigh_omitted" not in d
