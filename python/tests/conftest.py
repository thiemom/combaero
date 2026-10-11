"""Every converged network solve in the suite is also a conservation test (#481).

``NetworkSolver.solve`` is wrapped for the session: a solve that reports
success must close the global energy balance (see _closure_check.py). A test
that legitimately builds a non-closing network marks itself
``@pytest.mark.no_closure_check``.
"""

from __future__ import annotations

import os

import pytest

_ACTIVE = {"on": True}


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "no_closure_check: skip the global energy/mass closure assertion"
    )
    from _closure_check import closure

    from combaero.network import solver as solver_mod

    original = solver_mod.NetworkSolver.solve
    if getattr(original, "_closure_checked", False):
        return

    def solve(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        if (
            _ACTIVE["on"]
            and result.get("__success__", False)
            and "__x_solution__" in result
            and not getattr(self, "_in_closure_check", False)
        ):
            self._in_closure_check = True
            try:
                c = closure(self, result)
            finally:
                self._in_closure_check = False
            if c is not None and abs(c["energy"]) > c["energy_tol"]:
                raise AssertionError(
                    f"converged solve does not conserve energy: imbalance "
                    f"{c['energy']:.6g} W (tolerance {c['energy_tol']:.3g} W, scale "
                    f"{c['scale']:.6g} W, node mass imbalance {c['mass']:.3g} kg/s) "
                    f"in {os.environ.get('PYTEST_CURRENT_TEST', '?')}"
                )
        return result

    solve._closure_checked = True
    solver_mod.NetworkSolver.solve = solve


@pytest.fixture(autouse=True)
def _closure_marker(request: pytest.FixtureRequest):
    _ACTIVE["on"] = request.node.get_closest_marker("no_closure_check") is None
    yield
    _ACTIVE["on"] = True
