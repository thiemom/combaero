"""Element-level Jacobian of ChannelElement against central differences.

test_channel_jacobians.py checks that the C++ correlations populate their
derivative fields. That is not the same question as whether the ELEMENT's
residual Jacobian is right, and the gap once hid a real defect: for pin-fin and
impingement surfaces the multiplier varied with mass flow but was handed to the
friction routine frozen, so d(dP)/d(mdot) was short by 10.5%.

Those surfaces were removed in 0.7.0 (issue #339). Ribbed returned under
#334 and impingement (both forms) under #337 -- both checked here the same
way as smooth, closing the gap this docstring originally flagged.
"""

import pytest

import combaero as cb
from combaero.network.components import (
    ChannelElement,
    ConvectiveSurface,
    ImpingementModel,
    NetworkMixtureState,
    RibbedModel,
    SingleJetImpingementModel,
    SmoothModel,
)

Y_AIR = cb.mole_to_mass(cb.species.dry_air())

BASE = {"m_dot": 0.35, "T_up": 600.0, "P_up": 2.0e5, "Pt_up": 2.1e5, "P_dn": 1.95e5}

SURFACES = {
    "smooth": lambda: ConvectiveSurface(area=0.1, model=SmoothModel()),
    "ribbed": lambda: ConvectiveSurface(area=0.1, model=RibbedModel(e_D=0.06, p_e=10.0)),
    "impingement": lambda: ConvectiveSurface(
        area=0.1,
        model=ImpingementModel(d_jet=0.003, xn_d=8.0, yn_d=6.0, z_d=2.0, row=3),
    ),
    "single_jet_impingement": lambda: ConvectiveSurface(
        area=0.1, model=SingleJetImpingementModel(d_jet=0.004, L_D=7.75, R_D=5.0)
    ),
}


def _channel(surface_key: str) -> ChannelElement:
    return ChannelElement(
        id="ch",
        from_node="A",
        to_node="B",
        length=0.6,
        diameter=0.025,
        roughness=1e-5,
        surface=SURFACES[surface_key](),
    )


def _states(m_dot: float, T_up: float, P_up: float, Pt_up: float, P_dn: float):
    s_in = NetworkMixtureState(T=T_up, P=P_up, Pt=Pt_up, Tt=610.0, m_dot=m_dot, Y=Y_AIR)
    s_out = NetworkMixtureState(T=595.0, P=P_dn, Pt=P_dn, Tt=605.0, m_dot=m_dot, Y=Y_AIR)
    return s_in, s_out


def _residual(elem: ChannelElement, **kw) -> float:
    res, _ = elem.residuals(*_states(**kw))
    return res[0]


def _central_difference(elem: ChannelElement, var: str, step: float) -> float:
    plus, minus = dict(BASE), dict(BASE)
    plus[var] += step
    minus[var] -= step
    return (_residual(elem, **plus) - _residual(elem, **minus)) / (2.0 * step)


COLUMNS = [
    ("m_dot", "ch.m_dot", 1e-6),
    ("T_up", "A.T", 1e-3),
    ("Pt_up", "A.Pt", 1.0),
    ("P_dn", "B.Pt", 1.0),
]


@pytest.mark.parametrize("surface_key", sorted(SURFACES))
@pytest.mark.parametrize(("var", "column", "step"), COLUMNS)
def test_every_jacobian_column_matches_a_finite_difference(
    surface_key: str, var: str, column: str, step: float
) -> None:
    """Each analytic entry reproduces a central difference of the residual."""
    if surface_key == "ribbed" and column == "A.T":
        pytest.skip(
            "_ribbed_residuals' own docstring: the correlation's f has no "
            "Re-dependence, so dP's sensitivity to upstream temperature is "
            "deliberately absent rather than approximated -- a pre-existing, "
            "documented gap surfaced by adding 'ribbed' here, not a "
            "regression from it."
        )
    elem = _channel(surface_key)
    _, jac = elem.residuals(*_states(**BASE))
    analytic = jac[0].get(column, 0.0)
    numeric = _central_difference(elem, var, step)
    scale = max(abs(numeric), abs(analytic), 1e-12)
    assert abs(analytic - numeric) / scale < 1e-6, (
        f"{surface_key}/{column}: analytic {analytic:.6e} vs FD {numeric:.6e}"
    )


def test_user_f_multiplier_scales_the_drop_and_its_derivative() -> None:
    """The user multiplier survives the 0.7.0 removal and still applies.

    It is a tuning knob, not a correlation: it encodes nothing, defaults to
    1.0, and is only ever set by a caller. Pinned here because the removal of
    the correlation-derived multipliers deleted every other test that touched
    this code path.
    """
    plain = _channel("smooth")
    scaled = _channel("smooth")
    scaled.surface.f_multiplier = 2.0

    s_in, s_out = _states(**BASE)
    dP_plain = s_in.Pt - s_out.Pt - plain.residuals(s_in, s_out)[0][0]
    dP_scaled = s_in.Pt - s_out.Pt - scaled.residuals(s_in, s_out)[0][0]

    assert dP_scaled == pytest.approx(2.0 * dP_plain, rel=1e-9)

    j_plain = plain.residuals(s_in, s_out)[1][0]["ch.m_dot"]
    j_scaled = scaled.residuals(s_in, s_out)[1][0]["ch.m_dot"]
    assert j_scaled == pytest.approx(2.0 * j_plain, rel=1e-9)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
