"""Element-level Jacobian of ChannelElement against central differences.

The existing test_channel_jacobians.py checks that the C++ correlations
populate their derivative fields (``assert result.dh_dmdot != 0.0``). That is
not the same question as whether the ELEMENT's residual Jacobian is right, and
the gap hid a real defect: for pin-fin and impingement surfaces the surface
multiplier varies with mass flow, but was handed to the C++ friction routine
as a frozen constant, so d(dP)/d(mdot) was short by 10.5% and d(dP)/dT by
13.1%. Ribbed and dimpled were exact throughout -- their multipliers carry no
Reynolds dependence -- which is what localised the fault.
"""

import pytest

import combaero as cb
from combaero.network.components import (
    ChannelElement,
    ConvectiveSurface,
    DimpledModel,
    ImpingementModel,
    NetworkMixtureState,
    PinFinModel,
    RibbedModel,
    SmoothModel,
)

Y_AIR = cb.mole_to_mass(cb.species.dry_air())

# Kept inside every correlation's validated range so a range warning never
# masks a derivative error.
BASE = {"m_dot": 0.35, "T_up": 600.0, "P_up": 2.0e5, "Pt_up": 2.1e5, "P_dn": 1.95e5}

SURFACES = {
    "smooth": lambda: ConvectiveSurface(area=0.1, model=SmoothModel()),
    "ribbed": lambda: ConvectiveSurface(
        area=0.1, model=RibbedModel(e_D=0.05, pitch_to_height=10.0)
    ),
    "dimpled": lambda: ConvectiveSurface(area=0.1, model=DimpledModel(d_Dh=0.2, h_d=0.2, S_d=1.5)),
    "pin_fin": lambda: ConvectiveSurface(
        area=0.1,
        model=PinFinModel(pin_diameter=0.003, channel_height=0.006, S_D=2.5, X_D=2.5),
    ),
    "impingement": lambda: ConvectiveSurface(
        area=0.1,
        model=ImpingementModel(d_jet=0.002, z_D=3.0, x_D=5.0, y_D=5.0, A_target=0.05, Cd_jet=0.8),
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


# (perturbed variable, jacobian column, FD step)
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
    elem = _channel(surface_key)
    _, jac = elem.residuals(*_states(**BASE))
    analytic = jac[0].get(column, 0.0)
    numeric = _central_difference(elem, var, step)
    scale = max(abs(numeric), abs(analytic), 1e-12)
    assert abs(analytic - numeric) / scale < 1e-6, (
        f"{surface_key}/{column}: analytic {analytic:.6e} vs FD {numeric:.6e}"
    )


@pytest.mark.parametrize("surface_key", ["pin_fin", "impingement"])
def test_array_drop_is_the_correlation_drop(surface_key: str) -> None:
    """A localised array's element drop IS its correlation's dP.

    Pins the formulation: the array correlations return a pressure drop
    outright, so the element takes it rather than converting it into a
    multiplier on pipe friction. Routing it through a locally restated f_base
    cancelled exactly only when C++ picked the same correlation -- with
    friction_model="petukhov" it moved the drop by -24.7%.
    """
    for friction_model in ("haaland", "colebrook", "serghides", "petukhov"):
        elem = _channel(surface_key)
        elem.friction_model = friction_model
        s_in, s_out = _states(**BASE)
        res, _ = elem.residuals(s_in, s_out)
        dP_element = s_in.Pt - s_out.Pt - res[0]
        expected = elem.surface.f_multiplier * elem.htc_and_T(s_in).dP
        assert dP_element == pytest.approx(expected, rel=1e-12), friction_model


def test_ddP_dmdot_is_not_the_element_mass_flow() -> None:
    """ChannelResult.ddP_dmdot uses the ARRAY's flow area, not the element's.

    Chaining it directly is wrong by the area ratio, which is why the element
    chains ddP_dvelocity instead. This test exists so the trap is stated
    somewhere a reader will find it: the two differ by 45x for this geometry.
    """
    elem = _channel("pin_fin")
    s_in, _ = _states(**BASE)
    result = elem.htc_and_T(s_in)

    rho = s_in.density()
    correct = result.ddP_dvelocity / (rho * elem.area)
    naive = result.ddP_dmdot

    # A_cross = channel_height * pin_diameter * (S_D - 1) / S_D, from
    # heat_transfer.cpp -- the pin array's minimum section.
    a_cross = 0.006 * 0.003 * (2.5 - 1.0) / 2.5
    assert naive / correct == pytest.approx(elem.area / a_cross, rel=1e-6)
    assert naive / correct > 40.0, "the trap should be large enough to notice"


@pytest.mark.xfail(
    strict=True,
    reason=(
        "The pin-fin and impingement correlations expose no dP sensitivity to "
        "static pressure. dP scales as 1/rho, so the true column is dP/P -- "
        "about 0.058 here, a first-order term. Left absent rather than "
        "approximated from a pipe-friction stand-in."
    ),
)
def test_array_jacobian_carries_a_static_pressure_column() -> None:
    elem = _channel("pin_fin")
    _, jac = elem.residuals(*_states(**BASE))
    numeric = _central_difference(elem, "P_up", 1.0)
    assert numeric != 0.0, "the residual really does depend on P_up"
    assert jac[0].get("A.P", 0.0) == pytest.approx(numeric, rel=1e-6)


def test_reverse_flow_drop_opposes_the_flow() -> None:
    """Friction acts against the flow, so the drop follows sign(m_dot)."""
    elem = _channel("pin_fin")
    forward = dict(BASE)
    reverse = dict(BASE, m_dot=-BASE["m_dot"])

    s_in_f, s_out_f = _states(**forward)
    s_in_r, s_out_r = _states(**reverse)
    dP_f = s_in_f.Pt - s_out_f.Pt - elem.residuals(s_in_f, s_out_f)[0][0]
    dP_r = s_in_r.Pt - s_out_r.Pt - elem.residuals(s_in_r, s_out_r)[0][0]

    assert dP_f > 0.0
    assert dP_r == pytest.approx(-dP_f, rel=1e-12)
    # Magnitude depends on |m_dot|, so the derivative does not flip sign.
    assert elem.residuals(s_in_f, s_out_f)[1][0]["ch.m_dot"] == pytest.approx(
        elem.residuals(s_in_r, s_out_r)[1][0]["ch.m_dot"], rel=1e-12
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
