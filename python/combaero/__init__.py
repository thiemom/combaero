"""Public Python API for the combaero package.

This module re-exports selected functions from the compiled _core extension.
It is written to work both when running from the source tree and when
combaero is installed as a wheel.
"""

from __future__ import annotations

from importlib import import_module

try:  # Python 3.8+
    from importlib.metadata import PackageNotFoundError, version as _pkg_version
except ModuleNotFoundError:  # pragma: no cover - very old Python
    from importlib_metadata import (  # type: ignore[import-not-found]
        PackageNotFoundError,
        version as _pkg_version,
    )


def _load_version() -> str:
    """Return the installed distribution version, or a safe fallback.

    When running from a source tree without an installed wheel, the
    distribution metadata may not be available; in that case we expose a
    placeholder string instead of raising.
    """

    try:
        return _pkg_version("combaero")
    except PackageNotFoundError:
        return "0.0.0+local"


__version__: str = _load_version()


try:
    # Preferred: local extension in the same package (installed wheel or
    # in-tree build where _core was successfully built next to this file).
    from ._core import mixture_h, adiabatic_T_wgs  # type: ignore[attr-defined]
except ModuleNotFoundError:
    # Fallback: attempt to import from an installed combaero package that
    # already has _core available, then re-export the symbols.
    _core = import_module("combaero._core")
    mixture_h = _core.mixture_h
    adiabatic_T_wgs = _core.adiabatic_T_wgs


__all__ = ["mixture_h", "adiabatic_T_wgs", "__version__"]
