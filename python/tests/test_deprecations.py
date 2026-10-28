"""Every name deprecated in 0.7.0 is gone in 0.8.0; every 0.8.0 deprecation warns.

0.7.0 made each deprecated name warn and name its removal release; 0.8.0
removes them. This is the inventory: each removed function, alias, element,
keyword and strategy must fail the way a missing name fails (AttributeError,
ImportError, TypeError, ValueError), not keep working silently, and the
replacements must stay.

The second part is the inventory for the names deprecated in 0.8.0 (the
stagnation spellings, #782, removed in 0.9.0): each warns, names 0.9.0 and
its replacement, and agrees with the replacement bit for bit.
"""

from __future__ import annotations

import warnings

import pytest

import combaero as cb
import combaero.network as net
from combaero import _core, _renamed, _solver_tools
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
    VortexElement,
    pressure_loss,
)


@pytest.mark.parametrize(
    "name",
    ["standard_dry_air_composition", "smr_wgs_equilibrium", "smr_wgs_equilibrium_adiabatic"],
)
def test_the_removed_functions_are_gone(name: str) -> None:
    assert name not in cb.__all__
    with pytest.raises(AttributeError):
        getattr(cb, name)
    with pytest.raises(AttributeError):
        getattr(_core, name)


def test_the_replacement_functions_remain() -> None:
    assert len(cb.species.dry_air()) == cb.num_species()
    assert callable(cb.reforming_equilibrium)
    assert callable(cb.reforming_equilibrium_adiabatic)


@pytest.mark.parametrize("alias", ["ConstantLoss", "LinearThetaLoss"])
@pytest.mark.parametrize("module", [net, pressure_loss], ids=["network", "pressure_loss"])
def test_the_pressure_loss_aliases_are_gone(module: object, alias: str) -> None:
    with pytest.raises(AttributeError):
        getattr(module, alias)
    assert net.ConstantFractionLoss is pressure_loss.ConstantFractionLoss
    assert net.LinearThetaFractionLoss is pressure_loss.LinearThetaFractionLoss


def test_the_border_carnot_loss_element_is_gone() -> None:
    assert "BorderCarnotLossElement" not in net.__all__
    with pytest.raises(AttributeError):
        net.BorderCarnotLossElement  # noqa: B018
    with pytest.raises(ImportError):
        from combaero.network.components import BorderCarnotLossElement  # noqa: F401
    for module in (_core, _solver_tools):
        assert not hasattr(module, "border_carnot_loss_residual_and_jacobian")
    assert not hasattr(_core, "BorderCarnotLossResult")


@pytest.mark.parametrize(
    "kwargs",
    [{"omega_rpm": 3000.0}, {"r_c": 0.02}, {"n": 2.0}],
    ids=["omega_rpm", "r_c", "n"],
)
def test_the_removed_vortex_keywords_raise(kwargs: dict[str, float]) -> None:
    with pytest.raises(TypeError):
        VortexElement("v", "a", "b", **kwargs)


def test_the_vortex_omega_rpm_attribute_is_gone() -> None:
    el = VortexElement("v", "a", "b", shaft_speed_rpm=3000.0)
    assert el.shaft_speed_rpm == 3000.0
    with pytest.raises(AttributeError):
        el.omega_rpm  # noqa: B018


def test_the_vatistas_functions_remain() -> None:
    """Only the element's r_c/n went; the free functions are unchanged."""
    assert callable(cb.vatistas_v0_bar)


def test_incompressible_warmstart_raises_and_names_the_valid_strategies() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=200000.0, Tt=300.0))
    g.add_node(PressureBoundary("out", Pt=100000.0, Tt=300.0))
    g.add_element(OrificeElement("orf", "in", "out", 0.6, 0.05, correlation="fixed"))
    with pytest.raises(ValueError, match="outletref_warmstart") as info:
        NetworkSolver(g).solve(init_strategy="incompressible_warmstart", timeout=10.0)
    assert "'incompressible_warmstart'" in str(info.value)  # the rejected value
    for valid in ("default", "analytical_pt_prop", "homotopy", "continuation"):
        assert valid in str(info.value)


def test_the_warmstart_maxfev_keyword_went_with_it() -> None:
    with pytest.raises(TypeError):
        NetworkSolver(FlowNetwork()).solve(warmstart_maxfev=150)  # type: ignore[call-arg]


def test_an_unknown_network_attribute_still_raises() -> None:
    with pytest.raises(AttributeError):
        net.NoSuchElement  # noqa: B018


def test_a_gui_file_saved_with_incompressible_warmstart_still_loads() -> None:
    """The library refuses the strategy; a saved GUI file (or a headless
    payload) gets its outlet-referenced successor, as the frontend maps it."""
    from gui.backend.schemas import SolverSettings

    s = SolverSettings(init_strategy="incompressible_warmstart")  # type: ignore[arg-type]
    assert s.init_strategy == "outletref_warmstart"


# ---------------------------------------------------------------------------
# Deprecated in 0.8.0, removed in 0.9.0: the stagnation spellings (#782).
#
# P/T are static and Pt/Tt total. Each old spelling (a keyword, a function
# name, a struct field) must warn, name 0.9.0 and its replacement, and give
# the replacement's answer bit for bit. The inventory is combaero._renamed;
# every entry in it has a case here.
# ---------------------------------------------------------------------------


REMOVAL_782 = "removed in 0.9.0"

_X = list(cb.species.dry_air())
_Y = list(cb.mole_to_mass(_X))
_TT, _PT = 640.0, 3.1e5
_NOZZLE_GEOM = {
    "A_inlet": 3e-3,
    "A_throat": 1e-3,
    "A_exit": 2.5e-3,
    "x_throat": 0.05,
    "x_exit": 0.15,
}
_UP = {"m_dot": 0.4, "P_total_up": _PT, "P_static_up": 3.0e5, "Tt_up": _TT, "Y_up": _Y}


def _calls() -> dict[str, dict[str, object]]:
    """A valid all-keyword call, in the NEW names, for every renamed-kwarg function."""
    tp = {"Tt": _TT, "Pt": _PT}
    return {
        "nozzle_flow": {**tp, "P_back": 2.2e5, "A_eff": 1.3e-3, "X": _X},
        "solve_A_eff_from_mdot": {**tp, "P_back": 2.2e5, "mdot_target": 0.7, "X": _X},
        "solve_P_back_from_mdot": {**tp, "A_eff": 1.3e-3, "mdot_target": 0.5, "X": _X},
        "solve_Pt_from_mdot": {
            "Tt": _TT,
            "P_back": 2.2e5,
            "A_eff": 1.3e-3,
            "mdot_target": 0.7,
            "X": _X,
        },
        "critical_pressure_ratio": {**tp, "X": _X},
        "mach_from_pressure_ratio": {**tp, "P": 2.4e5, "X": _X},
        "mass_flux_isentropic": {**tp, "P": 2.4e5, "X": _X},
        "nozzle_cd": {**tp, "P_exit": 0.4e5, **_NOZZLE_GEOM, "X": _X},
        "nozzle_thrust_cd": {**tp, "P_design": 0.4e5, "P_amb": 1e5, **_NOZZLE_GEOM, "X": _X},
        "duct_static_state": {**tp, "G": 300.0, "X": _X},
        "duct_htc_static_state": {**tp, "G": 300.0, "X": _X},
        "duct_htc_static_state_dY": {
            **tp,
            "G": 300.0,
            "X": _X,
            "state": cb.duct_htc_static_state(_TT, _PT, 300.0, _X),
        },
        "sonic_static_pressure": {"Tt": _TT, "G": 700.0, "Y": _Y},
        "orifice_compressible_mdot_and_jacobian": {
            **tp,
            "P_back": 2.0e5,
            "X": _X,
            "Cd": 0.62,
            "area": 1e-3,
        },
        "T_from_stagnation": {"Tt": _TT, "M": 0.6, "X": _X},
        "P_from_stagnation": {**tp, "M": 0.6, "X": _X},
        "approach_static_pressure": {**tp, "G": 250.0, "Y": _Y, "incompressible": False},
        "orifice_residuals_and_jacobian": {**_UP, "P_static_down": 2.5e5, "Cd": 0.62, "area": 1e-3},
        "orifice_compressible_residuals_and_jacobian": {
            **_UP,
            "P_static_down": 2.5e5,
            "Cd": 0.62,
            "area": 1e-3,
        },
        "channel_residuals_and_jacobian": {
            **_UP,
            "P_static_down": 2.5e5,
            "L": 1.2,
            "D": 0.03,
            "roughness": 1e-5,
            "friction_model": "haaland",
        },
        "heat_entropy_flux": {
            "m_dot": 0.8,
            "Tt_in": 600.0,
            "Tt_out": 900.0,
            "Y": _Y,
            "M_in": 0.1,
            "M_out": 0.3,
        },
    }


_OLD_FIELDS = {old for mapping, _ in _renamed.FIELDS.values() for old in mapping}


def _snapshot(obj: object, depth: int = 0) -> object:
    """Every number a result carries, for an exact comparison."""
    if obj is None or isinstance(obj, bool | int | float | str):
        return obj
    if isinstance(obj, list | tuple):
        return [_snapshot(v, depth + 1) for v in obj]
    if depth > 3:
        return repr(type(obj))
    out = {}
    for name in dir(obj):
        if name.startswith("_") or name in _OLD_FIELDS:
            continue
        value = getattr(obj, name)
        if callable(value):
            continue
        out[name] = _snapshot(value, depth + 1)
    return out


def _renamed_call(fn_name: str) -> dict[str, object]:
    new_kwargs = _calls()[fn_name]
    back = {new: old for old, new in _renamed.KWARGS[fn_name].items()}
    return {back.get(k, k): v for k, v in new_kwargs.items()}


def _function(name: str) -> object:
    return getattr(cb, name, None) or getattr(_solver_tools, name)


def test_the_782_inventory_is_covered() -> None:
    assert set(_calls()) == set(_renamed.KWARGS)


@pytest.mark.parametrize("fn_name", sorted(_renamed.KWARGS))
def test_an_old_keyword_warns_and_agrees_bit_for_bit(fn_name: str) -> None:
    fn = _function(fn_name)
    new = fn(**_calls()[fn_name])
    old_kwargs = _renamed_call(fn_name)
    for old, new_name in _renamed.KWARGS[fn_name].items():
        assert old in old_kwargs
        with pytest.warns(DeprecationWarning, match=REMOVAL_782) as rec:
            got = fn(**old_kwargs)
        messages = " ".join(str(w.message) for w in rec)
        assert f"{old}=" in messages and f"{new_name}=" in messages
    assert _snapshot(got) == _snapshot(new)


@pytest.mark.parametrize("fn_name", sorted(_renamed.KWARGS))
def test_old_and_new_keyword_together_raise(fn_name: str) -> None:
    kwargs = _calls()[fn_name]
    old, new = next(iter(_renamed.KWARGS[fn_name].items()))
    with pytest.raises(TypeError, match=old):
        _function(fn_name)(**kwargs, **{old: kwargs[new]})


@pytest.mark.parametrize("fn_name", sorted(_renamed.KWARGS))
def test_the_new_keywords_do_not_warn(fn_name: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        _function(fn_name)(**_calls()[fn_name])


def test_every_access_path_reaches_the_same_wrapped_function() -> None:
    for name in _renamed.KWARGS:
        core_fn = getattr(_core, name)
        assert getattr(core_fn, "__renamed_kwargs__", None) == _renamed.KWARGS[name]
        for module in (cb, _solver_tools):
            if hasattr(module, name):
                assert getattr(module, name) is core_fn


_FUNCTION_ARGS: dict[str, tuple[object, ...]] = {
    "T0_from_static": (560.0, 0.6, _X),
    "T0_from_static_v": (560.0, 280.0, _X),
    "P0_from_static": (2.5e5, 560.0, 0.6, _X),
    "solve_P0_from_mdot": (_TT, 2.2e5, 1.3e-3, 0.7, _X),
    "T0_from_static_and_jacobian_M": (560.0, 0.6, _X),
    "P0_from_static_and_jacobian_M": (2.5e5, 560.0, 0.6, _X),
}


def test_the_782_function_inventory_is_covered() -> None:
    assert set(_FUNCTION_ARGS) == set(_renamed.FUNCTIONS)


@pytest.mark.parametrize(("old", "new"), sorted(_renamed.FUNCTIONS.items()))
def test_an_old_function_name_warns_and_agrees_bit_for_bit(old: str, new: str) -> None:
    args = _FUNCTION_ARGS[old]
    with pytest.warns(DeprecationWarning, match=REMOVAL_782) as rec:
        got = _function(old)(*args)
    assert new in str(rec[0].message)
    assert _snapshot(got) == _snapshot(_function(new)(*args))


def test_the_old_function_names_stay_exported_until_0_9() -> None:
    for old in ("T0_from_static", "T0_from_static_v", "P0_from_static", "solve_P0_from_mdot"):
        assert old in cb.__all__
    for new in ("Tt_from_static", "Tt_from_static_v", "Pt_from_static", "solve_Pt_from_mdot"):
        assert new in cb.__all__


def _instances() -> dict[str, object]:
    s1 = cb.MassStream(0.8, 600.0, 3.0e5, _Y)
    s2 = cb.MassStream(0.5, 400.0, 3.1e5, _Y)
    mix = _solver_tools.mixer_from_streams_and_jacobians([s1, s2], 1500.0, 0.0)
    calls = _calls()
    return {
        "MixtureState": cb.MixtureState(2.5e5, 3.0e5, 560.0, 640.0, 0.7, _Y),
        "MassStream": s1,
        "StreamJacobian": mix.dT_mix_d_stream[0],
        "NozzleSolution": cb.nozzle_cd(**calls["nozzle_cd"]),
        "DuctStaticState": cb.duct_static_state(_TT, _PT, 300.0, _X),
        "SonicStatic": _solver_tools.sonic_static_pressure(_TT, 700.0, _Y),
        "OrificeResult": _solver_tools.orifice_residuals_and_jacobian(
            **calls["orifice_residuals_and_jacobian"]
        ),
        "ChannelSolverResult": _solver_tools.channel_residuals_and_jacobian(
            **calls["channel_residuals_and_jacobian"]
        ),
    }


def test_the_782_field_inventory_is_covered() -> None:
    assert set(_instances()) == set(_renamed.FIELDS)


_FIELD_CASES = [
    (cls, old, new, writable)
    for cls, (mapping, writable) in sorted(_renamed.FIELDS.items())
    for old, new in sorted(mapping.items())
]


@pytest.mark.parametrize(
    ("cls", "old", "new", "writable"), _FIELD_CASES, ids=[f"{c}.{o}" for c, o, _, _ in _FIELD_CASES]
)
def test_an_old_field_warns_and_reads_the_new_one(
    cls: str, old: str, new: str, writable: bool
) -> None:
    obj = _instances()[cls]
    assert type(obj).__name__ == cls
    with pytest.warns(DeprecationWarning, match=REMOVAL_782) as rec:
        got = getattr(obj, old)
    assert f"{cls}.{new}" in str(rec[0].message)
    assert _snapshot(got) == _snapshot(getattr(obj, new))
    if writable:
        with pytest.warns(DeprecationWarning, match=REMOVAL_782):
            setattr(obj, old, 123.25)
        assert getattr(obj, new) == 123.25
    else:
        with pytest.raises(AttributeError):
            setattr(obj, old, 1.0)


def test_mixture_state_old_constructor_keywords_warn_and_agree() -> None:
    new = cb.MixtureState(P=2.5e5, Pt=3.0e5, T=560.0, Tt=640.0, m_dot=0.7, Y=_Y)
    with pytest.warns(DeprecationWarning, match=REMOVAL_782):
        old = cb.MixtureState(P=2.5e5, P_total=3.0e5, T=560.0, T_total=640.0, m_dot=0.7, Y=_Y)
    assert (old.P, old.Pt, old.T, old.Tt, old.m_dot) == (new.P, new.Pt, new.T, new.Tt, new.m_dot)
    assert old.density() == new.density()
    with pytest.raises(TypeError, match="P_total"):
        cb.MixtureState(Pt=3.0e5, P_total=3.0e5)


def test_compressible_nozzle_flow_old_keywords_warn_and_agree() -> None:
    from combaero import compressible

    new = compressible.nozzle_flow(Tt=_TT, Pt=_PT, X=_X, P_back=2.2e5, A_eff=1.3e-3)
    with pytest.warns(DeprecationWarning, match=REMOVAL_782):
        old = compressible.nozzle_flow(T=_TT, P=_PT, X=_X, P_back=2.2e5, A_eff=1.3e-3)
    assert _snapshot(old) == _snapshot(new)
    with pytest.raises(TypeError):
        compressible.nozzle_flow(T=_TT, Tt=_TT, Pt=_PT, X=_X, P_back=2.2e5, A_eff=1.3e-3)


def test_the_new_782_names_do_not_warn() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        for cls, obj in _instances().items():
            for new in _renamed.FIELDS[cls][0].values():
                getattr(obj, new)
        for new, args in zip(_renamed.FUNCTIONS.values(), _FUNCTION_ARGS.values(), strict=True):
            _function(new)(*args)
        cb.MixtureState(Pt=3.0e5, Tt=640.0)
