import importlib


def test_import_combaero_package() -> None:
    """Basic smoke test that the top-level combaero package is importable."""
    mod = importlib.import_module("combaero")
    assert mod is not None


def test_import_core_module() -> None:
    """Basic smoke test that the internal _core extension can be imported."""
    core = importlib.import_module("combaero._core")
    assert core is not None


def test_no_module_imports_a_dev_only_dependency() -> None:
    """Every module of the package imports without sympy (a dev-only
    dependency; the MSVC smoke job installs numpy and scipy only). Run in a fresh interpreter so another test's import of sympy
    cannot hide a leak."""
    import subprocess
    import sys

    probe = (
        "import importlib, pkgutil, sys\n"
        "import combaero\n"
        "for m in pkgutil.walk_packages(combaero.__path__, 'combaero.'):\n"
        "    importlib.import_module(m.name)\n"
        "print(sorted(n for n in ('sympy', 'mpmath') if n in sys.modules))\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, check=True
    ).stdout.strip()
    assert out == "[]", f"dev-only modules imported by the package: {out}"
