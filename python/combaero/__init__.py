"""Public Python API for the combaero package.

This module re-exports selected functions from the compiled _core extension.
It is written to work both when running from the source tree and when
combaero is installed as a wheel.
"""

try:
    # Preferred: local extension in the same package (installed wheel or
    # in-tree build where _core was successfully built next to this file).
    from ._core import mixture_h, adiabatic_T_wgs  # type: ignore[attr-defined]
except ModuleNotFoundError:
    # Fallback: attempt to import from an installed combaero package that
    # already has _core available, then re-export the symbols.
    import importlib

    _core = importlib.import_module("combaero._core")
    mixture_h = _core.mixture_h
    adiabatic_T_wgs = _core.adiabatic_T_wgs

__all__ = ["mixture_h", "adiabatic_T_wgs"]
