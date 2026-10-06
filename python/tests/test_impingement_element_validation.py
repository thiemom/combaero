"""The impingement ELEMENT, driven with Florschuetz et al. (1981)'s own flow (#460).

Correlation-level scoring cannot see how an element turns its mass flow into a
jet Reynolds number: the jet-array runner scores Nu against Re_j directly. So
this drives a ChannelElement carrying an ImpingementModel with the flow that
gives the paper's Re_j, and scores the element's Nu = h d / k against Fig. 5's
absolute first-row data (Re_j = 1e4, inline).

Before #460 the element took the jet flow as rho v * ConvectiveSurface.area, so
its Re_j scaled with A_surface / A_cross: on this configuration it missed Fig. 5
by +268% on average (+45% to +762%). The correlation itself, at the paper's
Re_j, misses by 5.6%. The element must now match the correlation.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb
from combaero.network.components import (
    ChannelElement,
    ConvectiveSurface,
    ImpingementModel,
    NetworkMixtureState,
)
from validation.cooling.schema import DATA_ROOT

X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)
T, P = 300.0, 1.2e5
D_JET = 0.00254  # Florschuetz's 0.1 in holes
SPAN = 0.122  # array span [m]; the physics must not depend on it


def _fig5() -> list[tuple[float, float, float, float]]:
    """(xn_d, yn_d, z_d, Nu1) for every digitised Fig. 5 point."""
    out = []
    for path in sorted((DATA_ROOT / "florschuetz1981").glob("fig5_xn_*_zd_*.csv")):
        parts = path.stem.split("_")
        xn, z = float(parts[2]), float(parts[4])
        for line in path.read_text().splitlines()[1:]:
            yn, nu = map(float, line.split(","))
            out.append((xn, yn, z, nu))
    return out


def _element_nu(xn_d: float, yn_d: float, z_d: float, Re_j: float = 1.0e4) -> float:
    """One element = row 1 of a span-wide array, carrying that row's jet flow."""
    mu = cb.complete_state(T, P, X_AIR).transport.mu
    k = cb.complete_state(T, P, X_AIR).transport.k
    xn, yn, z = xn_d * D_JET, yn_d * D_JET, z_d * D_JET
    n_holes = SPAN / yn
    m_dot = n_holes * Re_j * mu * math.pi * D_JET / 4.0  # Re_j = 4 m/(pi d mu)
    # The crossflow channel: z x SPAN. ChannelElement's area comes from its
    # diameter, so choose the diameter that gives that area; Dh is the
    # rectangle's hydraulic diameter.
    area = z * SPAN
    el = ChannelElement(
        id="row1",
        from_node="a",
        to_node="b",
        length=xn,
        diameter=math.sqrt(4.0 * area / math.pi),
        Dh=2.0 * z * SPAN / (z + SPAN),
        surface=ConvectiveSurface(
            area=xn * SPAN,
            model=ImpingementModel(d_jet=D_JET, xn_d=xn_d, yn_d=yn_d, z_d=z_d, row=1),
        ),
        t_hot=400.0,
    )
    state = NetworkMixtureState(T=T, P=P, Pt=P, Tt=T, m_dot=m_dot, Y=Y_AIR)
    return el.htc_and_T(state).h * D_JET / k


def test_fig5_is_all_there() -> None:
    assert len(_fig5()) == 27


def test_the_element_reproduces_fig5_like_the_correlation_does() -> None:
    """Fidelity THROUGH the element: the paper's own absolute data."""
    err = np.array([_element_nu(xn, yn, z) / nu - 1.0 for xn, yn, z, nu in _fig5()])
    assert abs(err.mean()) < 0.03, err.mean()
    assert np.abs(err).mean() < 0.08, np.abs(err).mean()
    assert np.abs(err).max() < 0.20, np.abs(err).max()


def test_the_element_equals_the_correlation_at_the_papers_re_j() -> None:
    """The element adds bookkeeping, not physics: at the paper's Re_j it must
    return the correlation's Nu exactly (row 1, Gc/Gj = 0)."""
    s = cb.florschuetz_1981_inline()
    Pr = cb.complete_state(T, P, X_AIR).transport.Pr
    for xn, yn, z, _ in _fig5():
        nu = cb.jet_array_impingement_nu(s, 1.0e4, 0.0, Pr, xn, yn, z).Nu
        assert _element_nu(xn, yn, z) == pytest.approx(nu, rel=1e-9)


def test_the_answer_does_not_depend_on_the_arbitrary_span() -> None:
    """SPAN only sets how many holes and how much flow; per hole nothing moves.
    Before #460 the element's Nu scaled with it through the area ratio."""
    global SPAN
    base = _element_nu(8.0, 6.0, 2.0)
    saved = SPAN
    try:
        SPAN = 3.0 * saved
        assert _element_nu(8.0, 6.0, 2.0) == pytest.approx(base, rel=1e-9)
    finally:
        SPAN = saved
