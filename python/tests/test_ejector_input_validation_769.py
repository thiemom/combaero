"""Ejector kernels refuse non-finite and out-of-domain input (#769).

src/ejector.cpp validated nothing: a NaN area ratio gave a finite residual and
ejector_huang_mixing_coefficient(nan) returned 0.84. Every public entry point
now raises ValueError (std::invalid_argument) for:

- non-finite input in any argument;
- parameters outside their physical domain (gamma > 1, r_gas > 0, nozzle
  efficiencies in (0, 1], recovery > 0, 1 < A_p1/A_t <= A_3/A_t, positive
  areas, loss coefficients and eps_frac >= 0, smootherstep windows lo < hi);
- in the standalone closures, a non-positive supply stagnation state.

The two whole-element residual functions take Newton iterates: those are
checked for finiteness only, because a trial point legitimately carries a
negative pressure or a reversed flow (measured on the 'ejector' stress
family, where the element saw p_g down to -3e-3 Pa and P_py to -5e7 Pa).
"""

import math
from collections.abc import Callable
from typing import Any

import pytest

from combaero import _core
from combaero.network.ejector_element import EjectorElement

NAN = math.nan
INF = math.inf

# Valid points (recorded from the suite's own ejector calls).
SUPPLY = (501324.3, 293.15, 101323.2, 293.15)  # p_g, t_g, p_e, t_e
GEOM = (1.8595, 5.3039)  # A_p1/A_t, A_3/A_t
GAS = (1.4009, 287.04)  # gamma, r_gas

CLOSURES_9 = [
    "ejector_critical_back_pressure_and_jacobian",
    "ejector_dead_head_back_pressure_and_jacobian",
    "ejector_huang_critical_back_pressure_and_jacobian",
    "ejector_sudden_expansion_dead_head_and_jacobian",
]

ELEMENT_ARGS = (
    0.00987,  # mp
    0.00692,  # ms
    0.01678,  # mdot_out
    501325.0,  # p_g
    293.15,  # t_g
    101325.0,  # p_e
    293.15,  # t_e
    103132.6,  # p_out
    132900.9,  # p_py
    *GEOM,
    8.553e-06,  # area_throat
    1.5904e-05,  # area_nozzle_exit
    2.9460e-05,  # area_secondary
    *GAS,
    0.95,  # eta_primary
    0.85,  # eta_secondary
    1.0,  # recovery_efficiency
    0.001,  # eps_frac
    0.9,  # s_choke_lo
    0.999,  # s_choke_hi
    0.98,  # s_sub_lo
    1.02,  # s_sub_hi
)

INCOMPRESSIBLE_ARGS = (
    0.3469,  # mp
    0.3469,  # ms
    0.6938,  # mo
    401325.0,  # pt_p
    0.2060,  # v_p
    1325.0,  # pt_s
    1.5432,  # v_s
    403975.0,  # pt_o
    101325.0,  # p3
    1.5904e-05,  # area_nozzle_exit
    2.9460e-05,  # area_secondary
    4.5365e-05,  # area_mixing
    0.0526,  # k_nozzle
    0.1765,  # k_suction
    0.0,  # k_throat
    0.0,  # k_diffuser
    1.0,  # recovery
)

# (name, function, valid arguments)
CASES: list[tuple[str, Callable[..., Any], tuple[Any, ...]]] = [
    *[(n, getattr(_core, n), (*SUPPLY, *GEOM, *GAS, 1.0)) for n in CLOSURES_9],
    (
        "ejector_entrainment_ratio_and_jacobian",
        _core.ejector_entrainment_ratio_and_jacobian,
        (*SUPPLY, *GEOM, GAS[0]),
    ),
    (
        "ejector_choked_mass_flow_and_jacobian",
        _core.ejector_choked_mass_flow_and_jacobian,
        (5e5, 293.15, 3.14e-5, *GAS, 0.95),
    ),
    (
        "ejector_cd_nozzle_mass_flow_and_jacobian",
        _core.ejector_cd_nozzle_mass_flow_and_jacobian,
        (101325.0, 300.0, 40000.0, 3.14e-5, 1e-4, 1.4, 287.0, 0.95, 1e-3),
    ),
    (
        "ejector_cd_nozzle_exit_static_and_jacobian",
        _core.ejector_cd_nozzle_exit_static_and_jacobian,
        (0.005, 101325.0, 300.0, 3.14e-5, 1e-4, 1.4, 287.0, 0.95, 1e-3),
    ),
    (
        "ejector_jetpump_discharge_and_jacobian",
        _core.ejector_jetpump_discharge_and_jacobian,
        (700000.0, 440.0, 22800.0, 280.0, 6840.0, 0.79, 1.401, 1.0, 0.95, 0.85),
    ),
    ("ejector_huang_mixing_coefficient", _core.ejector_huang_mixing_coefficient, (6.0,)),
    (
        "ejector_element_residuals_and_jacobian",
        _core.ejector_element_residuals_and_jacobian,
        ELEMENT_ARGS,
    ),
    (
        "ejector_jetpump_incompressible_residuals_and_jacobian",
        _core.ejector_jetpump_incompressible_residuals_and_jacobian,
        INCOMPRESSIBLE_ARGS,
    ),
]


def _non_finite_params() -> list[Any]:
    out = []
    for name, fn, args in CASES:
        for i in range(len(args)):
            for bad in (NAN, INF, -INF):
                out.append(pytest.param(fn, args, i, bad, id=f"{name}[{i}]={bad}"))
    return out


@pytest.mark.parametrize(("fn", "args", "index", "bad"), _non_finite_params())
def test_non_finite_argument_is_refused(
    fn: Callable[..., Any], args: tuple[Any, ...], index: int, bad: float
) -> None:
    fn(*args)  # the valid point itself answers
    bad_args = list(args)
    bad_args[index] = bad
    with pytest.raises(ValueError):
        fn(*bad_args)


def _with(args: tuple[Any, ...], index: int, value: float) -> list[Any]:
    out = list(args)
    out[index] = value
    return out


CLOSURE = (*SUPPLY, *GEOM, *GAS, 1.0)

# (id, function, valid arguments, index, out-of-domain value)
DOMAIN: list[tuple[str, Callable[..., Any], tuple[Any, ...], int, float]] = [
    # Standalone closures: supply stagnation state must be positive.
    *[
        (f"{n} {lbl}", getattr(_core, n), CLOSURE, i, v)
        for n in CLOSURES_9
        for i, lbl, v in (
            (0, "p_g=0", 0.0),
            (1, "t_g<0", -1.0),
            (2, "p_e<0", -1.0),
            (3, "t_e=0", 0.0),
            (4, "A_p1/A_t=1", 1.0),
            (5, "A_3<A_p1", 1.5),
            (6, "gamma=1", 1.0),
            (7, "r_gas=0", 0.0),
            (8, "recovery=0", 0.0),
        )
    ],
    (
        "choked eta>1",
        _core.ejector_choked_mass_flow_and_jacobian,
        (5e5, 293.15, 3.14e-5, *GAS, 0.95),
        5,
        1.2,
    ),
    (
        "choked area<0",
        _core.ejector_choked_mass_flow_and_jacobian,
        (5e5, 293.15, 3.14e-5, *GAS, 0.95),
        2,
        -1e-5,
    ),
    (
        "cd nozzle A_e<A_t",
        _core.ejector_cd_nozzle_mass_flow_and_jacobian,
        (101325.0, 300.0, 40000.0, 3.14e-5, 1e-4, 1.4, 287.0, 0.95, 1e-3),
        4,
        1e-5,
    ),
    (
        "cd nozzle eps<0",
        _core.ejector_cd_nozzle_mass_flow_and_jacobian,
        (101325.0, 300.0, 40000.0, 3.14e-5, 1e-4, 1.4, 287.0, 0.95, 1e-3),
        8,
        -1e-3,
    ),
    (
        "exit static m_dot<0",
        _core.ejector_cd_nozzle_exit_static_and_jacobian,
        (0.005, 101325.0, 300.0, 3.14e-5, 1e-4, 1.4, 287.0, 0.95, 1e-3),
        0,
        -0.001,
    ),
    (
        "jet pump omega<0",
        _core.ejector_jetpump_discharge_and_jacobian,
        (700000.0, 440.0, 22800.0, 280.0, 6840.0, 0.79, 1.401, 1.0, 0.95, 0.85),
        5,
        -0.1,
    ),
    (
        "jet pump eta_secondary=0",
        _core.ejector_jetpump_discharge_and_jacobian,
        (700000.0, 440.0, 22800.0, 280.0, 6840.0, 0.79, 1.401, 1.0, 0.95, 0.85),
        9,
        0.0,
    ),
    ("phi_m A_3/A_t=1", _core.ejector_huang_mixing_coefficient, (6.0,), 0, 1.0),
    # Element kernel: parameters only.
    ("element A_p1/A_t<1", _core.ejector_element_residuals_and_jacobian, ELEMENT_ARGS, 9, 0.9),
    ("element A_3<A_p1", _core.ejector_element_residuals_and_jacobian, ELEMENT_ARGS, 10, 1.5),
    ("element area_throat=0", _core.ejector_element_residuals_and_jacobian, ELEMENT_ARGS, 11, 0.0),
    ("element gamma<1", _core.ejector_element_residuals_and_jacobian, ELEMENT_ARGS, 14, 0.9),
    ("element eta_p>1", _core.ejector_element_residuals_and_jacobian, ELEMENT_ARGS, 16, 1.05),
    ("element recovery<0", _core.ejector_element_residuals_and_jacobian, ELEMENT_ARGS, 18, -1.0),
    ("element s_choke lo>=hi", _core.ejector_element_residuals_and_jacobian, ELEMENT_ARGS, 20, 1.0),
    ("element s_sub lo>=hi", _core.ejector_element_residuals_and_jacobian, ELEMENT_ARGS, 23, 0.98),
    (
        "incompressible area_mixing=0",
        _core.ejector_jetpump_incompressible_residuals_and_jacobian,
        INCOMPRESSIBLE_ARGS,
        11,
        0.0,
    ),
    (
        "incompressible k_throat<0",
        _core.ejector_jetpump_incompressible_residuals_and_jacobian,
        INCOMPRESSIBLE_ARGS,
        14,
        -0.1,
    ),
    (
        "incompressible recovery=0",
        _core.ejector_jetpump_incompressible_residuals_and_jacobian,
        INCOMPRESSIBLE_ARGS,
        16,
        0.0,
    ),
]


@pytest.mark.parametrize(
    ("fn", "args", "index", "value"),
    [pytest.param(fn, a, i, v, id=name) for name, fn, a, i, v in DOMAIN],
)
def test_out_of_domain_parameter_is_refused(
    fn: Callable[..., Any], args: tuple[Any, ...], index: int, value: float
) -> None:
    fn(*args)
    with pytest.raises(ValueError):
        fn(*_with(args, index, value))


def test_element_kernels_accept_out_of_domain_iterates() -> None:
    # The values the 'ejector' stress family's solves pass at trial points:
    # negative supply and mixing-plane pressures, reversed flows. The element
    # answers (through its C1 floors); refusing would turn a rejected probe
    # into a failed evaluation.
    probes = {0: -0.0023, 1: -1.27, 3: -0.0028, 5: -7.9e4, 7: -5.2e5, 8: -5.2e7}
    for i, v in probes.items():
        rj = _core.ejector_element_residuals_and_jacobian(*_with(ELEMENT_ARGS, i, v))
        assert len(rj.residuals) == 4
    for i, v in {0: -0.29, 1: -7.7, 5: -2.2e6, 7: -5.7e5, 8: -5.3e6}.items():
        rj = _core.ejector_jetpump_incompressible_residuals_and_jacobian(
            *_with(INCOMPRESSIBLE_ARGS, i, v)
        )
        assert all(math.isfinite(r) for r in rj.residuals)


def test_normal_shock_limit_is_kept() -> None:
    # A_3 = A_p1 (no suction annulus) is the closures' normal-shock limit.
    a = GEOM[0]
    res = _core.ejector_sudden_expansion_dead_head_and_jacobian(*SUPPLY, a, a, *GAS, 1.0)
    assert math.isfinite(res.value.p_c_pa)


def _element(**kw: Any) -> EjectorElement:
    args: dict[str, Any] = {
        "id": "ej",
        "primary_node": "p",
        "secondary_node": "s",
        "outlet_node": "o",
        "throat_area": 3.14e-5,
        "nozzle_exit_area": 1e-4,
        "mixing_area": 8e-4,
    }
    args.update(kw)
    return EjectorElement(**args)


@pytest.mark.parametrize(
    "kw",
    [
        {"throat_area": NAN},
        {"throat_area": INF},
        {"nozzle_exit_area": NAN},
        {"nozzle_exit_area": INF},
        {"mixing_area": NAN},
        {"mixing_area": INF},
        {"recovery_efficiency": NAN},
        {"recovery_efficiency": INF},
        {"k_throat": NAN},
        {"k_diffuser": INF},
    ],
    ids=lambda kw: ",".join(f"{k}={v}" for k, v in kw.items()),
)
def test_element_constructor_refuses_non_finite(kw: dict[str, float]) -> None:
    _element()  # the valid element builds
    with pytest.raises(ValueError):
        _element(**kw)


ABOVE_ONE = math.nextafter(1.0, 2.0)

# (id, function, valid arguments, index of the recovery argument)
RECOVERY: list[tuple[str, Callable[..., Any], tuple[Any, ...], int]] = [
    *[(n, getattr(_core, n), CLOSURE, 8) for n in CLOSURES_9],
    (
        "ejector_jetpump_discharge_and_jacobian",
        _core.ejector_jetpump_discharge_and_jacobian,
        (700000.0, 440.0, 22800.0, 280.0, 6840.0, 0.79, 1.401, 1.0, 0.95, 0.85),
        7,
    ),
    (
        "ejector_element_residuals_and_jacobian",
        _core.ejector_element_residuals_and_jacobian,
        ELEMENT_ARGS,
        18,
    ),
    (
        "ejector_jetpump_incompressible_residuals_and_jacobian",
        _core.ejector_jetpump_incompressible_residuals_and_jacobian,
        INCOMPRESSIBLE_ARGS,
        16,
    ),
]


@pytest.mark.parametrize(
    ("fn", "args", "index"), [pytest.param(fn, a, i, id=n) for n, fn, a, i in RECOVERY]
)
def test_recovery_above_one_is_refused(
    fn: Callable[..., Any], args: tuple[Any, ...], index: int
) -> None:
    # Recovery multiplies the lossless mixed stagnation pressure: above 1 the
    # ejector would gain pressure beyond lossless mixing. 1 itself answers.
    fn(*_with(args, index, 1.0))
    with pytest.raises(ValueError, match="recovery"):
        fn(*_with(args, index, ABOVE_ONE))


def test_element_constructor_caps_recovery_at_one() -> None:
    assert _element(recovery_efficiency=1.0).recovery_efficiency == 1.0
    with pytest.raises(ValueError, match="recovery_efficiency"):
        _element(recovery_efficiency=ABOVE_ONE)


def test_gui_schema_caps_recovery_at_one() -> None:
    from pydantic import ValidationError

    from gui.backend.schemas import EjectorData

    assert EjectorData(recovery_efficiency=1.0).recovery_efficiency == 1.0
    for bad in (ABOVE_ONE, 0.0, NAN):
        with pytest.raises(ValidationError, match="recovery_efficiency"):
            EjectorData(recovery_efficiency=bad)
