"""Pin-fin arrays as a ConvectiveSurface and in ChannelElement (#335).

The C++ tests pin the correlations; these pin the composition the element
adds: Vmax from the array geometry, the fin model, the base-area reference,
the per-row drop and its Jacobian, and the inline-friction refusal.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network.components import (
    ChannelElement,
    ConvectiveSurface,
    NetworkMixtureState,
    PinFinModel,
    SmoothModel,
)

X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)
STATE = {"T": 700.0, "P": 1.5e6, "X": X_AIR}
DH = 0.009  # channel hydraulic diameter [m]
FLOW = {"velocity": 30.0, "diameter": DH, "length": 0.05, "T_hot": 1100.0}


def _surface(**kw) -> ConvectiveSurface:
    return ConvectiveSurface(area=0.01, model=PinFinModel(**kw))


def _run(surface: ConvectiveSurface, **flow):
    return surface.htc_and_T(**STATE, **{**FLOW, **flow})


def test_defaults_are_a_working_example_inside_the_metzger_box() -> None:
    r = _run(_surface())
    assert r.h > 0.0
    assert not r.extrapolated, "the default geometry must sit inside its own sets' box"
    assert 0.0 < r.eta_fin <= 1.0


def test_the_element_composes_the_correlation_exactly() -> None:
    """Vmax, Re_D, Nu and the base-area fin model, recomputed by hand."""
    model = PinFinModel()
    r = _run(_surface())
    geom = model.geometry()
    rho = cb.density(STATE["T"], STATE["P"], STATE["X"])
    tr = cb.complete_state(STATE["T"], STATE["P"], STATE["X"]).transport
    v_max = FLOW["velocity"] / cb.pin_fin_amin_over_afrontal(geom)
    Re = rho * v_max * model.pin_diameter / tr.mu
    Nu = cb.evaluate_pin_fin_nu(cb.metzger_1986_staggered_nu(), geom, Re, tr.Pr).Nu
    h = Nu * tr.k / model.pin_diameter
    frac = cb.pin_fin_area_fractions(geom)
    eta = cb.pin_fin_array_efficiency(
        h, model.k_pin, model.pin_diameter, model.H_D * model.pin_diameter, frac.pin_over_total
    ).eta_fin
    assert r.Re == pytest.approx(Re, rel=1e-9)
    assert r.h_array == pytest.approx(h, rel=1e-9)
    assert r.h == pytest.approx(h * (frac.endwall_exposed + eta * frac.pin), rel=1e-9)
    f = cb.evaluate_pin_fin_friction(cb.metzger_1982_staggered_friction(), geom, Re).f
    assert r.dP == pytest.approx(2.0 * rho * v_max**2 * model.N_rows * f, rel=1e-9)


def test_conducting_pins_lose_effectiveness_and_isothermal_pins_do_not() -> None:
    iso = _run(_surface(k_pin=math.inf))
    alloy = _run(_surface(k_pin=20.0))
    ceramic = _run(_surface(k_pin=2.0))
    assert iso.eta_fin == pytest.approx(1.0, abs=1e-12)
    assert iso.h > alloy.h > ceramic.h


def test_inline_defaults_to_chyu_friction_never_to_staggered() -> None:
    """Chyu (1990) Fig. 6 is the only inline friction source in hand."""
    model = PinFinModel(arrangement=cb.PinArrangement.Inline, N_rows=7)
    assert model.resolved_f_set().name == "chyu_1990_inline_straight_friction"
    r = _run(ConvectiveSurface(area=0.01, model=model))
    assert r.f == pytest.approx(0.1693 / 4.0, rel=1e-12)
    # A staggered friction set on an inline array is refused, not used.
    with pytest.raises(ValueError, match="not substituted"):
        PinFinModel(
            arrangement=cb.PinArrangement.Inline,
            f_set=cb.metzger_1982_staggered_friction(),
        )


def test_inline_with_a_user_friction_set_and_the_chyu_transfer() -> None:
    mine = cb.metzger_1982_staggered_friction()
    mine.name, mine.source = "rig_inline_f", "measured"
    mine.provenance = cb.RibProvenance.User
    mine.arrangement = cb.PinArrangement.Inline
    model = PinFinModel(
        arrangement=cb.PinArrangement.Inline,
        f_set=mine,
        modifier=cb.chyu_1998_inline_over_staggered(),
    )
    assert model.resolved_nu_set().arrangement == cb.PinArrangement.Staggered
    stag = _run(_surface())
    inl = _run(ConvectiveSurface(area=0.01, model=model))
    ratio = cb.evaluate_pin_fin_modifier(
        cb.chyu_1998_inline_over_staggered(), model.geometry(), inl.Re
    ).ratio_Nu
    assert inl.Nu == pytest.approx(stag.Nu * ratio, rel=1e-9)


def test_a_modifier_from_the_wrong_arrangement_is_rejected() -> None:
    mine = cb.metzger_1982_staggered_friction()
    mine.name = "x"
    with pytest.raises(ValueError, match="transfers to"):
        PinFinModel(modifier=cb.chyu_1998_inline_over_staggered(), f_set=mine)


def test_wall_coupling_derivative_matches_finite_difference() -> None:
    """dh/dmdot on the channel mass flow rho v pi Dh^2/4 (channel_smooth's)."""
    s = _surface()
    v0, dv = 30.0, 1e-4
    r0 = _run(s, velocity=v0)
    dh_dv = (_run(s, velocity=v0 + dv).h - _run(s, velocity=v0 - dv).h) / (2 * dv)
    rho = cb.density(STATE["T"], STATE["P"], STATE["X"])
    dh_dmdot_fd = dh_dv / (rho * math.pi / 4.0 * DH * DH)
    assert r0.dh_dmdot == pytest.approx(dh_dmdot_fd, rel=1e-4)
    for field in ("dh_dmdot", "dh_dT", "dT_aw_dmdot", "dT_aw_dT"):
        assert isinstance(getattr(r0, field), float)


# ---------------------------------------------------------------------------
# ChannelElement residual path
# ---------------------------------------------------------------------------


def _element(model) -> ChannelElement:
    return ChannelElement(
        id="ch",
        from_node="A",
        to_node="B",
        length=0.05,
        diameter=DH,
        surface=ConvectiveSurface(area=0.01, model=model),
    )


def _states(m_dot: float, T: float = 700.0, P: float = 1.5e6):
    a = NetworkMixtureState(T=T, P=P, Pt=P + 2e4, Tt=T + 2.0, m_dot=m_dot, Y=Y_AIR)
    b = NetworkMixtureState(T=T - 3.0, P=P - 5e4, Pt=P - 5e4, Tt=T, m_dot=m_dot, Y=Y_AIR)
    return a, b


@pytest.mark.parametrize("m_dot", [0.002, 0.02, 0.08])
def test_drop_jacobian_matches_central_differences(m_dot: float) -> None:
    el = _element(PinFinModel())
    _, jac = el.residuals(*_states(m_dot))

    def r(**kw) -> float:
        return el.residuals(*_states(**{"m_dot": m_dot, **kw}))[0][0]

    h = 1e-7
    assert jac[0]["ch.m_dot"] == pytest.approx(
        (r(m_dot=m_dot + h) - r(m_dot=m_dot - h)) / (2 * h), rel=1e-5
    )
    dT = 1e-2
    assert jac[0]["A.T"] == pytest.approx((r(T=700.0 + dT) - r(T=700.0 - dT)) / (2 * dT), rel=2e-3)
    dPs = 10.0
    fd_P = (r(P=1.5e6 + dPs) - r(P=1.5e6 - dPs)) / (2 * dPs)
    # The residual's Pt terms move with P in _states; strip them (d/dP = 0).
    assert jac[0]["A.P"] == pytest.approx(fd_P, rel=2e-3)


def test_drop_is_odd_and_its_derivative_even_in_the_flow() -> None:
    el = _element(PinFinModel())
    out = {}
    for m in (0.02, -0.02):
        a, b = _states(m)
        res, jac = el.residuals(a, b)
        out[m] = (a.Pt - b.Pt - res[0], jac[0]["ch.m_dot"])
    assert out[0.02][0] == pytest.approx(-out[-0.02][0], rel=1e-12)
    assert out[0.02][1] == pytest.approx(out[-0.02][1], rel=1e-12)


def test_drop_is_the_per_row_definition_not_a_quarter_of_it() -> None:
    """The code #332 removed used rho V^2/2: a factor 4 low against
    Metzger's f = dP / (2 rho Vmax^2 N). Pin the definition."""
    model = PinFinModel()
    el = _element(model)
    m = 0.02
    a, b = _states(m)
    dP = a.Pt - b.Pt - el.residuals(a, b)[0][0]
    rho = cb.density(a.T, a.P, cb.mass_to_mole(Y_AIR))
    a_min = el.area * cb.pin_fin_amin_over_afrontal(model.geometry())
    v_max = m / (rho * a_min)
    mu = cb.complete_state(a.T, a.P, cb.mass_to_mole(Y_AIR)).transport.mu
    Re = m * model.pin_diameter / (a_min * mu)
    f = cb.evaluate_pin_fin_friction(cb.metzger_1982_staggered_friction(), model.geometry(), Re).f
    assert dP == pytest.approx(2.0 * rho * v_max**2 * model.N_rows * f, rel=1e-6)


def test_friction_multiplier_scales_the_drop() -> None:
    a, b = _states(0.02)
    plain = _element(PinFinModel())
    scaled = _element(PinFinModel())
    scaled.surface.f_multiplier = 1.5
    dp0 = a.Pt - b.Pt - plain.residuals(a, b)[0][0]
    dp1 = a.Pt - b.Pt - scaled.residuals(a, b)[0][0]
    assert dp1 == pytest.approx(1.5 * dp0, rel=1e-12)


def test_a_pin_fin_channel_solves_in_a_network() -> None:
    from combaero.network import FlowNetwork, NetworkSolver
    from combaero.network.components import PressureBoundary

    def solve(model) -> dict:
        g = FlowNetwork()
        a = PressureBoundary("A")
        a.Pt, a.Tt, a.Y = 1.52e6, 700.0, Y_AIR
        b = PressureBoundary("B")
        b.Pt, b.Tt, b.Y = 1.50e6, 700.0, Y_AIR
        g.add_node(a)
        g.add_node(b)
        g.add_element(_element(model))
        return NetworkSolver(g).solve()

    smooth = solve(SmoothModel())
    pins = solve(PinFinModel())
    assert smooth["__success__"], smooth.get("__message__")
    assert pins["__success__"], pins.get("__message__")
    assert 0.0 < pins["ch.m_dot"] < smooth["ch.m_dot"]
