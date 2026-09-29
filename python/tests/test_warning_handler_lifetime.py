"""A Python warning handler must not outlive the interpreter.

`set_warning_handler` stores its argument in a function-local static on the
C++ side. When that argument is a Python callable, the static owns a
`py::object`, and a function-local static is destroyed by `__cxa_atexit` --
which runs AFTER `Py_Finalize()`. Releasing a Python reference on a finalised
interpreter segfaults.

The symptom is invisible to an ordinary test: the process does all its work,
passes every assertion, prints its output, and only then dies during shutdown.
`pytest` itself exits 0 throughout, because it never installs a handler and
leaves it there. So these tests run real subprocesses and assert on the EXIT
CODE, which is the only place the defect shows.

Measured before the fix (exit 139, SIGSEGV):

    set_warning_handler(lambda m: None)          <- the docstring's own example
    get_warning_handler() then set it back
    with suppress_warnings(): pass

`set_warning_handler(None)` was always clean, which is what localised it to a
Python callable being owned by the static rather than to the handler
machinery in general.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap

import pytest

# Each case installs a Python callable and leaves it in place at exit, which
# is the shape that crashed. The names are what pytest prints on failure.
LEAVES_A_PYTHON_HANDLER_INSTALLED = {
    "set_lambda": """
        import combaero as cb
        cb.set_warning_handler(lambda msg: None)
        """,
    "set_named_function": """
        import combaero as cb

        def handler(msg):
            pass

        cb.set_warning_handler(handler)
        """,
    "get_then_restore": """
        import combaero as cb
        previous = cb.get_warning_handler()
        cb.set_warning_handler(lambda msg: None)
        cb.set_warning_handler(previous)
        """,
    "suppress_warnings": """
        import combaero as cb
        with cb.suppress_warnings():
            pass
        """,
    "suppress_warnings_nested": """
        import combaero as cb
        with cb.suppress_warnings():
            with cb.suppress_warnings():
                pass
        """,
    "handler_that_captures_state": """
        import combaero as cb
        seen = []
        cb.set_warning_handler(seen.append)
        """,
}


def run_script(body: str) -> subprocess.CompletedProcess[str]:
    """Run `body` in a fresh interpreter and return the completed process."""
    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(body)],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )


@pytest.mark.parametrize("case", sorted(LEAVES_A_PYTHON_HANDLER_INSTALLED))
def test_installed_handler_does_not_crash_at_shutdown(case: str) -> None:
    """A handler left installed must not take the process down at exit."""
    result = run_script(LEAVES_A_PYTHON_HANDLER_INSTALLED[case])
    assert result.returncode == 0, (
        f"{case!r} exited {result.returncode} "
        f"({'signal ' + str(result.returncode - 128) if result.returncode > 128 else 'error'}). "
        "A Python warning handler is being released after Py_Finalize(); the "
        "module teardown capsule in _core.cpp should have cleared it.\n"
        f"stderr: {result.stderr}"
    )


def test_the_documented_reset_path_still_works() -> None:
    """`set_warning_handler(None)` was the one clean path before the fix.

    Kept so a future change cannot fix the cases above by breaking this one.
    """
    result = run_script("""
        import combaero as cb
        cb.set_warning_handler(lambda msg: None)
        cb.set_warning_handler(None)
        """)
    assert result.returncode == 0, result.stderr


def test_the_handler_still_receives_warnings() -> None:
    """The teardown capsule must not disarm the handler early.

    Clearing the handler at module teardown would be a silent regression if it
    happened at any other moment -- warnings would simply stop arriving, and
    no crash would announce it. This runs a real warning end to end.
    """
    result = run_script("""
        import combaero as cb

        seen = []
        cb.set_warning_handler(seen.append)
        # Re = 100 is laminar, well below Dittus-Boelter's validated range,
        # so the correlation warns. Chosen because it has routed through
        # warn() for far longer than this fix, keeping the test independent
        # of any one caller.
        cb.nusselt_dittus_boelter(100.0, 0.7, True)

        assert len(seen) == 1, seen
        assert "dittus_boelter" in seen[0], seen
        print("OK")
        """)
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_suppress_warnings_actually_suppresses_and_then_restores() -> None:
    """The context manager's contract, across a process boundary.

    In-process coverage of this lives in test_heat_transfer.py; this variant
    exists because the restore path is what installs a Python callable into
    the static, so it belongs to the lifetime story too.
    """
    result = run_script("""
        import combaero as cb

        seen = []
        cb.set_warning_handler(seen.append)

        with cb.suppress_warnings():
            cb.nusselt_dittus_boelter(100.0, 0.7, True)
        assert seen == [], f"warning leaked out of suppress_warnings: {seen}"

        cb.nusselt_dittus_boelter(100.0, 0.7, True)
        assert len(seen) == 1, f"handler was not restored: {seen}"
        print("OK")
        """)
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout
