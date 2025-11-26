import importlib


def test_import_combaero_package() -> None:
    """Basic smoke test that the top-level combaero package is importable."""
    mod = importlib.import_module("combaero")
    assert mod is not None


def test_import_core_module() -> None:
    """Basic smoke test that the internal _core extension can be imported."""
    core = importlib.import_module("combaero._core")
    assert core is not None
