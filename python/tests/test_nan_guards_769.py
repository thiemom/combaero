"""Non-finite input is refused, not answered (#769).

A check written ``if (x <= 0) throw`` lets NaN through, because every
comparison with NaN is false: the function answered with NaN or, after a
clamp, a plausible number. Each public function below is called with NaN (and
+inf) in each numeric argument and must raise. The exception type is the one
the function already used for that argument's domain check -- ValueError
(std::invalid_argument), RuntimeError (std::runtime_error) or IndexError
(std::out_of_range, the saturation-pressure range check) -- so callers that
catch it keep working.
"""

import math
from collections.abc import Callable
from typing import Any

import pytest

import combaero as cb
from combaero import _core

Annulus = cb.geometry.Annulus
CanAnnularFlowGeometry = cb.geometry.CanAnnularFlowGeometry
Tube = cb.geometry.Tube

NAN = math.nan
INF = math.inf

# (name, function, valid arguments, exception per argument index or a single
# exception for all, indices skipped with the reason in the comment)
Case = tuple[str, Callable[..., Any], tuple[Any, ...], Any, tuple[int, ...]]

_HUMID_TPR = {0: IndexError, 1: RuntimeError, 2: RuntimeError}


def _method(cls: type, name: str) -> Callable[..., Any]:
    """Build cls(*args) and call its method name."""

    def call(*args: float) -> Any:
        return getattr(cls(*args), name)()

    return call


CASES: list[Case] = [
    # -- humid air ---------------------------------------------------------
    ("saturation_vapor_pressure", cb.saturation_vapor_pressure, (300.0,), IndexError, ()),
    ("vapor_pressure", cb.vapor_pressure, (300.0, 0.5), {0: IndexError, 1: RuntimeError}, ()),
    ("humidity_ratio", cb.humidity_ratio, (300.0, 101325.0, 0.5), _HUMID_TPR, ()),
    (
        "water_vapor_mole_fraction",
        cb.water_vapor_mole_fraction,
        (300.0, 101325.0, 0.5),
        _HUMID_TPR,
        (),
    ),
    ("humid_air_composition", cb.humid_air_composition, (300.0, 101325.0, 0.5), _HUMID_TPR, ()),
    # P is unused by dewpoint (kept for API symmetry), so it is not checked.
    ("dewpoint", cb.dewpoint, (300.0, 101325.0, 0.5), {0: IndexError, 2: ValueError}, (1,)),
    ("wet_bulb_temperature", cb.wet_bulb_temperature, (300.0, 101325.0, 0.5), _HUMID_TPR, ()),
    ("humid_air_enthalpy", cb.humid_air_enthalpy, (300.0, 101325.0, 0.5), _HUMID_TPR, ()),
    (
        "humid_air_enthalpy_nasa9",
        cb.humid_air_enthalpy_nasa9,
        (300.0, 101325.0, 0.5),
        _HUMID_TPR,
        (),
    ),
    ("humid_air_density", cb.humid_air_density, (300.0, 101325.0, 0.5), _HUMID_TPR, ()),
    (
        "HumidAir.set_TP_RH",
        lambda T, P, RH: cb.HumidAir().set_TP_RH(T, P, RH),
        (300.0, 101325.0, 0.5),
        {0: IndexError, 1: RuntimeError, 2: ValueError},
        (),
    ),
    (
        "HumidAir.set_TP_omega",
        lambda T, P, w: cb.HumidAir().set_TP_omega(T, P, w),
        (300.0, 101325.0, 0.01),
        {0: IndexError, 1: ValueError, 2: ValueError},
        (),
    ),
    # -- materials ---------------------------------------------------------
    ("k_inconel718", cb.k_inconel718, (800.0,), RuntimeError, ()),
    ("k_haynes230", cb.k_haynes230, (800.0,), RuntimeError, ()),
    ("k_stainless_steel_310", cb.k_stainless_steel_310, (800.0,), RuntimeError, ()),
    ("k_stainless_steel_316", cb.k_stainless_steel_316, (800.0,), RuntimeError, ()),
    ("k_aluminum_6061", cb.k_aluminum_6061, (400.0,), RuntimeError, ()),
    # hours = +inf is the fully sintered limit, so only NaN is refused there.
    ("k_tbc_ysz", cb.k_tbc_ysz, (1200.0, 100.0), RuntimeError, ()),
    (
        "get_material_conductivity",
        lambda T: cb.get_material_conductivity("inconel718", T),
        (800.0,),
        RuntimeError,
        (),
    ),
    # -- acoustics ---------------------------------------------------------
    (
        "helmholtz_frequency",
        cb.helmholtz_frequency,
        (1e-3, 1e-4, 0.01, 340.0, 0.85),
        ValueError,
        (),
    ),
    ("strouhal", cb.strouhal, (100.0, 0.01, 10.0), ValueError, ()),
    ("frequency_from_strouhal", cb.frequency_from_strouhal, (0.3, 0.01, 10.0), ValueError, ()),
    ("quarter_wave_frequency", cb.quarter_wave_frequency, (0.5, 340.0), ValueError, ()),
    ("half_wave_frequency", cb.half_wave_frequency, (0.5, 340.0), ValueError, ()),
    ("stokes_layer", cb.stokes_layer, (1.5e-5, 100.0), ValueError, ()),
    ("thermal_layer", cb.thermal_layer, (2.2e-5, 100.0), ValueError, ()),
    (
        "effective_viscothermal_layer",
        cb.effective_viscothermal_layer,
        (2e-4, 2.5e-4, 1.4),
        ValueError,
        (),
    ),
    (
        "helmholtz_Q",
        cb.helmholtz_Q,
        (1e-3, 1e-4, 0.01, 1.5e-5, 2.2e-5, 1.4, 100.0),
        ValueError,
        (),
    ),
    ("tube_Q", cb.tube_Q, (0.5, 0.05, 1.5e-5, 2.2e-5, 1.4, 100.0), ValueError, ()),
    ("damping_ratio", cb.damping_ratio, (10.0,), ValueError, ()),
    ("bandwidth", cb.bandwidth, (100.0, 10.0), ValueError, ()),
    ("axial_mode_upstream", cb.axial_mode_upstream, (100.0, 0.1), ValueError, ()),
    ("axial_mode_downstream", cb.axial_mode_downstream, (100.0, 0.1), ValueError, ()),
    ("axial_mode_split", cb.axial_mode_split, (100.0, 0.1), ValueError, ()),
    ("acoustic_properties", cb.acoustic_properties, (100.0, 1.2, 340.0, 1.0, 2e-5), ValueError, ()),
    ("wavelength", cb.wavelength, (100.0, 340.0), ValueError, ()),
    ("frequency_from_wavelength", cb.frequency_from_wavelength, (3.4, 340.0), ValueError, ()),
    ("acoustic_impedance", cb.acoustic_impedance, (1.2, 340.0), ValueError, ()),
    ("sound_pressure_level", cb.sound_pressure_level, (1.0, 2e-5), ValueError, ()),
    ("particle_velocity", cb.particle_velocity, (1.0, 1.2, 340.0), ValueError, ()),
    (
        "tube_axial_modes",
        lambda L, D, c: cb.tube_axial_modes(
            Tube(L, D), c, cb.BoundaryCondition.Closed, cb.BoundaryCondition.Open
        ),
        (0.5, 0.05, 340.0),
        ValueError,
        (),
    ),
    (
        "annulus_axial_modes",
        lambda L, Di, Do, c: cb.annulus_axial_modes(
            Annulus(L, Di, Do), c, cb.BoundaryCondition.Closed, cb.BoundaryCondition.Open
        ),
        (0.5, 0.4, 0.6, 340.0),
        ValueError,
        (),
    ),
    (
        "annulus_azimuthal_modes",
        lambda Di, Do, c: cb.annulus_azimuthal_modes(Annulus(0.5, Di, Do), c),
        (0.4, 0.6, 340.0),
        ValueError,
        (),
    ),
    (
        "orifice_impedance_with_flow",
        cb.orifice_impedance_with_flow,
        (500.0, 1.0, 1.0, 1e-3, 1e-3, 0.05, 0.7, 1.2, 340.0),
        RuntimeError,
        (),
    ),
    ("is_whistling_risk", cb.is_whistling_risk, (500.0, 5.0, 1e-3), RuntimeError, ()),
    # -- geometry ----------------------------------------------------------
    ("hydraulic_diameter", cb.hydraulic_diameter, (1e-3, 0.1), ValueError, ()),
    ("hydraulic_diameter_rect", cb.hydraulic_diameter_rect, (0.02, 0.05), ValueError, ()),
    ("hydraulic_diameter_annulus", cb.hydraulic_diameter_annulus, (0.6, 0.4), ValueError, ()),
    ("circular_area", cb.circular_area, (0.05,), ValueError, ()),
    ("annular_area", cb.annular_area, (0.6, 0.4), ValueError, ()),
    ("cylinder_volume", cb.cylinder_volume, (0.05, 0.5), ValueError, ()),
    ("residence_time", cb.residence_time, (1e-3, 1e-4), ValueError, ()),
    ("residence_time_mdot", cb.residence_time_mdot, (1e-3, 0.1, 1.2), ValueError, ()),
    ("space_velocity", cb.space_velocity, (1e-4, 1e-3), ValueError, ()),
    ("Tube.area", lambda L, D: Tube(L, D).area(), (0.5, 0.05), ValueError, (0,)),
    ("Tube.perimeter", lambda L, D: Tube(L, D).perimeter(), (0.5, 0.05), ValueError, (0,)),
    ("Tube.volume", lambda L, D: Tube(L, D).volume(), (0.5, 0.05), ValueError, ()),
    *[
        (
            f"Annulus.{m}",
            _method(Annulus, m),
            (0.5, 0.4, 0.6),
            ValueError,
            (0,),
        )
        for m in ("D_mean", "gap", "area", "circumference")
    ],
    (
        "Annulus.volume",
        lambda L, Di, Do: Annulus(L, Di, Do).volume(),
        (0.5, 0.4, 0.6),
        ValueError,
        (),
    ),
    *[
        (
            f"CanAnnularFlowGeometry.{m}",
            _method(CanAnnularFlowGeometry, m),
            (0.3, 0.4, 0.6, 0.2, 0.1, 0.05),
            ValueError,
            skip,
        )
        for m, skip in (
            ("D_mean", (0, 3, 4, 5)),
            ("gap", (0, 3, 4, 5)),
            ("area", (0, 3, 4, 5)),
            ("volume", (3, 4, 5)),
            ("circumference", (0, 3, 4, 5)),
            ("area_primary", (0, 1, 2, 3, 5)),
            ("volume_primary", (0, 1, 2, 5)),
            ("volume_transition", (0, 3)),
            ("volume_total", ()),
            ("length_total", (1, 2, 4)),
        )
    ],
    # -- elsewhere: the same 'x <= 0' pattern, made NaN-safe ----------------
    ("orifice_mdot", cb.orifice_mdot, (2e5, 1e5, 1e-4, 0.6, 1.2), ValueError, ()),
    ("orifice_velocity", cb.orifice_velocity, (2e5, 1e5, 1.2), ValueError, ()),
    ("orifice_area", cb.orifice_area, (0.1, 2e5, 1e5, 0.6, 1.2), ValueError, ()),
    # Nu is not range-checked (a correlation's output passes through).
    ("htc_from_nusselt", cb.htc_from_nusselt, (100.0, 0.03, 0.05), ValueError, (0,)),
    ("beta_from_diameters", cb.beta_from_diameters, (0.05, 0.1), ValueError, ()),
    ("friction_haaland e_D", cb.friction_haaland, (1e5, 1e-4), ValueError, (0,)),
]


def _params() -> list[Any]:
    out = []
    for name, fn, args, exc, skip in CASES:
        for i in range(len(args)):
            if i in skip:
                continue
            want = exc[i] if isinstance(exc, dict) else exc
            if want is None:
                continue
            for bad in (NAN, INF):
                if name == "k_tbc_ysz" and i == 1 and bad == INF:
                    continue  # +inf hours: the fully sintered limit
                out.append(pytest.param(fn, args, i, bad, want, id=f"{name}[{i}]={bad}"))
    return out


@pytest.mark.parametrize(("fn", "args", "index", "bad", "exc"), _params())
def test_non_finite_argument_is_refused(
    fn: Callable[..., Any], args: tuple[Any, ...], index: int, bad: float, exc: type
) -> None:
    fn(*args)  # the valid point itself answers
    bad_args = list(args)
    bad_args[index] = bad
    with pytest.raises(exc):
        fn(*bad_args)


# -- follow-up: NaN that a clamp, a root finder or a "no correction" branch
# turned into a plausible finite number. Found by replaying every numeric
# _core call the suite makes with NaN in each argument. Only the listed
# argument indices are refused; finite values keep their existing meaning
# (beta outside (0, 1) is still "no correction", D_h <= 0 still selects the
# equivalent circle, m_ref <= 0 still switches the ramp off).
AIR = cb.dry_air()

# (name, function, valid arguments, refused argument indices)
FollowUp = tuple[str, Callable[..., Any], tuple[Any, ...], tuple[int, ...]]

FOLLOW_UP: list[FollowUp] = [
    (
        "fanno_duct",
        lambda Pt, Tt, G, L, D, e, fm: cb.fanno_duct(Pt, Tt, G, AIR, L, D, e, "haaland", fm),
        (200000.0, 300.0, 278.0, 1.0, 0.02, 1e-5, 1.0),
        (3, 4, 5, 6),
    ),
    (
        "fanno_choked_mass_flux",
        lambda L, D, e, fm: cb.fanno_choked_mass_flux(200000.0, 300.0, AIR, L, D, e, "haaland", fm),
        (1.0, 0.02, 1e-5, 1.0),
        (0, 1, 2, 3),
    ),
    (
        "fanno_rayleigh_choked_mass_flux",
        lambda Pt0, Tt0, Tte, L, D, e, fm: cb.fanno_rayleigh_choked_mass_flux(
            Pt0, Tt0, Tte, AIR, L, D, e, "fixed", fm
        ),
        (200000.0, 600.0, 900.0, 0.5, 0.01, 0.0, 0.02),
        (0, 1, 2, 3, 4, 5, 6),
    ),
    (
        "fanno_rayleigh_march",
        lambda G, Pt0, Tt0, Tte, e, fm: cb.fanno_rayleigh_march(
            G, Pt0, Tt0, Tte, AIR, 0.5, 0.01, e, "fixed", fm
        ),
        (150.0, 200000.0, 600.0, 900.0, 0.0, 0.02),
        (0, 1, 2, 3, 4, 5),
    ),
    (
        "fanno_length_between",
        lambda D, e, fm: cb.fanno_length_between(
            100.0, 300.0, 0.125, 1.0, AIR, D, e, "haaland", fm
        ),
        (0.02, 1e-5, 1.0),
        (0, 1, 2),
    ),
    (
        "fanno_dx_dmach",
        lambda D, e, fm: cb.fanno_dx_dmach(300.0, 300.0, 0.2, AIR, D, e, "haaland", fm),
        (0.02, 1e-5, 1.0),
        (0, 1, 2),
    ),
    (
        "fanno_channel_flow",
        lambda L, D, e, fm: cb.fanno_channel_flow(
            104000.0, 600.0, AIR, 103980.0, True, L, D, e, "haaland", fm
        ),
        (0.019, 0.076, 0.0, 1.0),
        (0, 1, 2, 3),
    ),
    (
        "fanno_rayleigh_channel_flow",
        lambda L, D, e, fm: cb.fanno_rayleigh_channel_flow(
            102000.0, 1200.0, 1184.0, AIR, 100000.0, False, L, D, e, "haaland", fm
        ),
        (0.127, 0.05, 1e-5, 1.0),
        (0, 1, 2, 3),
    ),
    ("velocity_of_approach_factor", cb.velocity_of_approach_factor, (0.5,), (0,)),
    ("expansibility_factor kappa", cb.expansibility_factor, (0.5, 1000.0, 101325.0, 1.4), (3,)),
    ("jobson_force_defect", cb.jobson_force_defect, (0.6, 0.5, 1.0), (0, 1, 2)),
    ("jobson_ci_from_cd", cb.jobson_ci_from_cd, (0.6, 0.2, 1.0), (0, 1, 2)),
    (
        "orifice_jobson_residuals_and_jacobian",
        lambda Cd, A, beta, Cv: _core.orifice_jobson_residuals_and_jacobian(
            0.0, 300000.0, 300000.0, 400.0, AIR, 152730.0, Cd, A, beta, Cv
        ),
        (0.61, 1e-4, 0.0, 0.98),
        (0, 1, 2, 3),
    ),
    (
        "orifice_iso5167_residuals_and_jacobian",
        lambda Cd, A, beta: _core.orifice_iso5167_residuals_and_jacobian(
            0.0, 200000.0, 200000.0, 600.0, AIR, 180000.0, Cd, A, beta
        ),
        (0.61, 1e-4, 0.0),
        (0, 1, 2),
    ),
    (
        "film_effectiveness_baldauf_2002 Tu",
        cb.film_effectiveness_baldauf_2002,
        (20.0, 2.0, 1.2, 30.0, 3.0, 0.015),
        (5,),
    ),
    (
        "film_effectiveness_baldauf_2002_and_derivatives Tu",
        cb.film_effectiveness_baldauf_2002_and_derivatives,
        (20.0, 2.0, 1.2, 30.0, 3.0, 0.015),
        (5,),
    ),
    (
        "film_span_average_baldauf_2002 Tu",
        cb.film_span_average_baldauf_2002,
        (10, 7.5, 0.0, 62.5, 2.1, 1.9, 30.0, 7.5, 0.05),
        (8,),
    ),
    # v = 0 is the branch where the 10 K floor answered for a NaN T_static.
    (
        "T_adiabatic_wall v=0",
        lambda Ts, v: cb.T_adiabatic_wall(Ts, v, 300.0, 101325.0, AIR),
        (300.0, 0.0),
        (0, 1),
    ),
    (
        "T_adiabatic_wall",
        lambda Ts, v: cb.T_adiabatic_wall(Ts, v, 300.0, 101325.0, AIR),
        (300.0, 50.0),
        (0, 1),
    ),
    (
        "area_change_loss conical",
        lambda F_in, F_out, L, D_h: cb.area_change_loss(
            0.264, 103643.0, 300.0, AIR, F_in, F_out, "conical", L, D_h, False
        ),
        (0.00196, 0.00503, 0.05, 0.0),
        (0, 1, 2, 3),
    ),
    (
        "area_change_flow",
        lambda F_in, F_out, D_h: cb.area_change_flow(
            200000.0, 300.0, AIR, 199000.0, F_in=F_in, F_out=F_out, model="sharp", D_h=D_h
        ),
        (0.001, 0.004, 0.0),
        (0, 1, 2),
    ),
    ("impingement_eta_r", cb.impingement_eta_r, (2.71, 5.0, 4.0, 1.0), (1, 2)),
    ("omega22_and_derivative delta_star", cb.omega22_and_derivative, (0.3, 0.0), (1,)),
    (
        "combustor_fraction_loss m_ref",
        cb.combustor_fraction_loss,
        (0.0, 0.0005, 0.001, 0.02, 0.001),
        (4,),
    ),
]


def _follow_up_params() -> list[Any]:
    out = []
    for name, fn, args, refused in FOLLOW_UP:
        for i in refused:
            for bad in (NAN, INF, -INF):
                out.append(pytest.param(fn, args, i, bad, id=f"{name}[{i}]={bad}"))
    return out


@pytest.mark.parametrize(("fn", "args", "index", "bad"), _follow_up_params())
def test_non_finite_argument_is_refused_follow_up(
    fn: Callable[..., Any], args: tuple[Any, ...], index: int, bad: float
) -> None:
    fn(*args)  # the valid point itself answers
    bad_args = list(args)
    bad_args[index] = bad
    with pytest.raises(ValueError):
        fn(*bad_args)


def test_finite_conventions_are_kept() -> None:
    # Finite values outside the formal range keep their documented meaning;
    # only non-finite values are new refusals.
    assert cb.velocity_of_approach_factor(0.0) == 1.0
    assert cb.velocity_of_approach_factor(1.5) == 1.0
    assert cb.expansibility_factor(0.5, 1000.0, 101325.0, 1.0) == 1.0
    a = cb.area_change_flow(200000.0, 300.0, AIR, 199000.0, F_in=0.001, F_out=0.004, D_h=0.0)
    b = cb.area_change_flow(
        200000.0, 300.0, AIR, 199000.0, F_in=0.001, F_out=0.004, D_h=(4 * 0.001 / math.pi) ** 0.5
    )
    assert a.m_dot == pytest.approx(b.m_dot, rel=1e-12)
    off = cb.combustor_fraction_loss(0.0, 0.0005, 0.001, 0.02, 0.0)
    assert off.xi == pytest.approx(0.02)
    assert cb.T_adiabatic_wall(5.0, 0.0, 300.0, 101325.0, AIR) == 10.0  # the 10 K floor


def test_swapped_annulus_is_refused() -> None:
    # A swapped annulus used to return a negative area (-0.306 m^2 here).
    swapped = Annulus(1.0, 0.8, 0.5)
    for method in ("area", "gap", "D_mean", "volume", "circumference"):
        with pytest.raises(ValueError, match="D_outer > D_inner"):
            getattr(swapped, method)()
    with pytest.raises(ValueError, match="D_outer > D_inner"):
        CanAnnularFlowGeometry(0.3, 0.6, 0.4, 0.2, 0.1, 0.05).area()
    with pytest.raises(ValueError, match="D_outer > D_inner"):
        cb.annular_area(0.5, 0.8)


def test_get_material_conductivity_still_answers_finite_T() -> None:
    # Finite T outside the range follows the envelope (#882), flat below
    # 200 K; only NaN, which used to come out as the 1e-3 floor, is refused.
    k_lo = cb.get_material_conductivity("inconel718", 1.0)
    assert k_lo == pytest.approx(cb.k_inconel718(300.0) - 50.0 * 0.0146, rel=1e-6)
    with pytest.raises(RuntimeError):
        cb.get_material_conductivity("inconel718", NAN)


def test_ntu_and_film_coefficient_keep_their_infinite_limits() -> None:
    # +inf NTU and h are physical limits, not bad input; only NaN is refused.
    assert cb.effectiveness_counterflow(INF, 0.5) == pytest.approx(1.0)
    assert cb.thermal_resistance(INF, 1.0) == 0.0
    with pytest.raises(ValueError):
        cb.effectiveness_counterflow(NAN, 0.5)
    with pytest.raises(ValueError):
        cb.thermal_resistance(NAN, 1.0)
