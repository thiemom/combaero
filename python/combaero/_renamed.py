"""Deprecated spellings of the stagnation names (#782); removed in 0.9.0.

The rule settled in #773/#777 for pressure now holds for temperature too:
``P``/``T`` are static, ``Pt``/``Tt`` are total (stagnation). The public API
used ``T0``/``P0``, ``T_total``/``P_total`` and a bare ``T`` or ``T_up`` for
the same total quantities. Those spellings still work until 0.9.0, each with
a ``DeprecationWarning`` that names its replacement:

- a keyword argument (``duct_static_state(T0=...)``): the call is forwarded
  with the new keyword. Positional calls are untouched, so the new and old
  spellings reach the same compiled function with the same values;
- a function name (``T0_from_static``): a forwarding wrapper;
- a struct field (``MixtureState.T_total``): a property on the pybind11 class
  that reads (and, where the field is writable, writes) the new field.

Passing an old keyword together with its replacement is a ``TypeError``.

pybind11 has no native keyword aliases. A second ``(*args, **kwargs)``
overload in C++ would catch the old keywords without a Python frame on the
new path, but it also catches every malformed call and replaces pybind11's
"incompatible function arguments" message, which lists the signatures, with
its own. A thin Python wrapper keeps that message and costs one frame per
call. ``install`` patches the ``_core`` module in place, before
``combaero/__init__.py`` re-exports from it, so every path that reaches these
functions (``combaero.X``, ``combaero._core.X``, ``combaero._solver_tools.X``)
sees the same object.
"""

from __future__ import annotations

import functools
import warnings
from collections.abc import Callable, Mapping
from types import ModuleType
from typing import Any

#: The release that removes every name deprecated here.
REMOVAL = "0.9.0"


def _warn(old: str, new: str, stacklevel: int = 3) -> None:
    warnings.warn(
        f"{old} is deprecated and will be removed in {REMOVAL}: use {new}.",
        DeprecationWarning,
        stacklevel=stacklevel,
    )


def _remap(qualname: str, mapping: Mapping[str, str], kwargs: dict[str, Any]) -> None:
    for old, new in mapping.items():
        if old in kwargs:
            if new in kwargs:
                raise TypeError(
                    f"{qualname}() got both '{old}' and its replacement '{new}'; pass '{new}' only"
                )
            _warn(f"{qualname}({old}=...)", f"{qualname}({new}=...)", stacklevel=4)
            kwargs[new] = kwargs.pop(old)


def renamed_kwargs(
    fn: Callable[..., Any], qualname: str, mapping: Mapping[str, str]
) -> Callable[..., Any]:
    """Wrap ``fn`` so the old keyword names in ``mapping`` warn and forward."""
    olds = frozenset(mapping)

    @functools.wraps(fn)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        if kwargs and not olds.isdisjoint(kwargs):
            _remap(qualname, mapping, kwargs)
        return fn(*args, **kwargs)

    wrapper.__renamed_kwargs__ = dict(mapping)  # type: ignore[attr-defined]
    return wrapper


def renamed_function(fn: Callable[..., Any], old: str, new: str) -> Callable[..., Any]:
    """A forwarding wrapper for a function's old name."""

    @functools.wraps(fn)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        _warn(old, new)
        return fn(*args, **kwargs)

    wrapper.__name__ = old
    wrapper.__qualname__ = old
    wrapper.__doc__ = f"Deprecated alias of {new}; removed in {REMOVAL}.\n\n{fn.__doc__ or ''}"
    return wrapper


def renamed_field(cls: type, old: str, new: str, writable: bool) -> None:
    """Give ``cls`` a property ``old`` that warns and reads/writes ``new``."""
    name = cls.__name__

    def fget(self: Any) -> Any:
        _warn(f"{name}.{old}", f"{name}.{new}")
        return getattr(self, new)

    def fset(self: Any, value: Any) -> None:
        _warn(f"{name}.{old}", f"{name}.{new}")
        setattr(self, new, value)

    doc = f"Deprecated alias of {new}; removed in {REMOVAL}."
    setattr(cls, old, property(fget, fset if writable else None, doc=doc))


def renamed_init_kwargs(cls: type, mapping: Mapping[str, str]) -> None:
    """Accept the old constructor keywords of ``cls`` with a warning."""
    init = cls.__init__
    olds = frozenset(mapping)
    qualname = cls.__name__

    @functools.wraps(init)
    def __init__(self: Any, *args: Any, **kwargs: Any) -> None:  # noqa: N807
        if kwargs and not olds.isdisjoint(kwargs):
            _remap(qualname, mapping, kwargs)
        init(self, *args, **kwargs)

    cls.__init__ = __init__  # type: ignore[misc]


_TP = {"T0": "Tt", "P0": "Pt"}
_EJ = {"p0": "Pt", "t0": "Tt"}

#: Keyword renames: function -> {old keyword: new keyword}.
KWARGS: dict[str, dict[str, str]] = {
    "nozzle_flow": _TP,
    "solve_A_eff_from_mdot": _TP,
    "solve_P_back_from_mdot": _TP,
    "solve_Pt_from_mdot": {"T0": "Tt"},
    "critical_pressure_ratio": _TP,
    "mach_from_pressure_ratio": _TP,
    "mass_flux_isentropic": _TP,
    "nozzle_cd": _TP,
    "nozzle_thrust_cd": _TP,
    "duct_static_state": _TP,
    "duct_htc_static_state": _TP,
    "duct_htc_static_state_dY": _TP,
    "sonic_static_pressure": {"T0": "Tt"},
    "orifice_compressible_mdot_and_jacobian": _TP,
    "T_from_stagnation": {"T0": "Tt"},
    "P_from_stagnation": _TP,
    "approach_static_pressure": {"P_total": "Pt", "T_total": "Tt"},
    "orifice_residuals_and_jacobian": {"T_up": "Tt_up"},
    "channel_residuals_and_jacobian": {"T_up": "Tt_up"},
    "orifice_compressible_residuals_and_jacobian": {"T_up": "Tt_up"},
    "heat_entropy_flux": {"T0_in": "Tt_in", "T0_out": "Tt_out"},
    # The head loss is total-total (Lefebvre; #782 audit).
    "combustor_head_loss": {"T_in": "Tt_in", "P_in": "Pt_in"},
    # The #782 audit's stragglers: totals still spelled the old way.
    "chamber_merge_face_state": {"T_main": "Tt_main", "T": "Tt"},
    "momentum_chamber_residual_and_jacobian": {"P_total": "Pt", "T": "Tt"},
    "orifice_tap_pressure_and_jacobian": {"P_total_up": "Pt_up", "P_total_down": "Pt_down"},
    "ejector_choked_mass_flow_and_jacobian": _EJ,
    "ejector_cd_nozzle_mass_flow_and_jacobian": _EJ,
    "ejector_cd_nozzle_exit_static_and_jacobian": _EJ,
}

#: Function renames: old name -> new name.
FUNCTIONS: dict[str, str] = {
    "T0_from_static": "Tt_from_static",
    "T0_from_static_v": "Tt_from_static_v",
    "P0_from_static": "Pt_from_static",
    "solve_P0_from_mdot": "solve_Pt_from_mdot",
    "T0_from_static_and_jacobian_M": "Tt_from_static_and_jacobian_M",
    "P0_from_static_and_jacobian_M": "Pt_from_static_and_jacobian_M",
}

#: Field renames: class -> ({old field: new field}, writable).
FIELDS: dict[str, tuple[dict[str, str], bool]] = {
    "MixtureState": ({"P_total": "Pt", "T_total": "Tt"}, True),
    "MassStream": ({"T": "Tt", "P_total": "Pt"}, True),
    "StreamJacobian": ({"d_T": "d_Tt", "d_P_total": "d_Pt"}, False),
    "NozzleSolution": ({"T0": "Tt", "P0": "Pt"}, False),
    "DuctStaticState": (
        {
            "dT_dT0": "dT_dTt",
            "dT_dP0": "dT_dPt",
            "dP_dT0": "dP_dTt",
            "dP_dP0": "dP_dPt",
            "drho_dT0": "drho_dTt",
            "drho_dP0": "drho_dPt",
            "dGeff_dT0": "dGeff_dTt",
            "dGeff_dP0": "dGeff_dPt",
        },
        False,
    ),
    "SonicStatic": ({"dP_dT0": "dP_dTt"}, False),
    "OrificeResult": ({"d_mdot_dT_up": "d_mdot_dTt_up"}, False),
    "ChannelSolverResult": ({"d_dP_dT_up": "d_dP_dTt_up"}, False),
    "CombustorLossResult": ({"dxi_dT_in": "dxi_dTt_in", "dxi_dP_in": "dxi_dPt_in"}, False),
    "MixerResult": (
        {
            "T_mix": "Tt_mix",
            "P_total_mix": "Pt_mix",
            "dT_mix_d_stream": "dTt_mix_d_stream",
            "dP_total_mix_d_stream": "dPt_mix_d_stream",
            "dT_mix_d_delta_h": "dTt_mix_d_delta_h",
            "dT_mix_dQ": "dTt_mix_dQ",
        },
        False,
    ),
}

#: Constructor keyword renames: class -> {old keyword: new keyword}.
INIT_KWARGS: dict[str, dict[str, str]] = {
    "MixtureState": {"P_total": "Pt", "T_total": "Tt"},
}


def install(core: ModuleType) -> None:
    """Patch the deprecated spellings into the compiled module, once."""
    if getattr(core, "__renamed_installed__", False):
        return
    for name, mapping in KWARGS.items():
        setattr(core, name, renamed_kwargs(getattr(core, name), name, mapping))
    for old, new in FUNCTIONS.items():
        setattr(core, old, renamed_function(getattr(core, new), old, new))
    for cls_name, (mapping, writable) in FIELDS.items():
        cls = getattr(core, cls_name)
        for old, new in mapping.items():
            renamed_field(cls, old, new, writable)
    for cls_name, mapping in INIT_KWARGS.items():
        renamed_init_kwargs(getattr(core, cls_name), mapping)
    core.__renamed_installed__ = True  # type: ignore[attr-defined]
