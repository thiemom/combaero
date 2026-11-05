"""A result's keys follow the network, never the operating point.

Every palette network of ``scripts/gui_result_keys.py`` is solved through the
GUI's ``/solve`` at its own point, at zero flow, reversed and choked (the
driving pressure difference, or the fed mass flow, scaled by 0.3, 0, -0.5 and 30).
The dropped node's numeric keys -- a value, or null for an undefined one (the
backend writes NaN as null) -- must be the same set at every point that
converges. A key that exists only when its value is nonzero, or only in one
regime, makes the Probe, the CSV columns and the key list depend on the last
bit of a solved flow: ``Q_withheld`` came and went with ``m * (1/m)``
rounding, and differed between macOS and Linux.

The test asserts keys, not convergence. Whether a hard point converges
depends on wall clock (the stall detector's window, the timeout): CI's
runners stalled two impingement arrays at 30x that converge on a desktop,
and a loaded VM converged the listed zero-flow orifices in a pipe that fail
serially. So every point that converges is compared, listed or not, and only
the network's own point plus one other must converge for the check to mean
anything.

An undefined value is NaN in the Python result, null in the GUI's JSON.
"""

from __future__ import annotations

import copy
import sys
import warnings
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "scripts"))

import gui_result_keys as keys  # noqa: E402
import test_gui_http_routes as http  # noqa: E402

#: own point, part load (so a network with a declared direction or no
#: zero-flow root still has a second point), zero flow, reversed, choked
FACTORS = (1.0, 0.3, 0.0, -0.5, 30.0)

#: (kind, variant, factor) solves that do not converge, and why. Known
#: failures, documented, and given a short budget; one that converges anyway
#: (a faster or slower runner) is compared like any other point.
EXPECTED_FAILURES: dict[tuple[str, str, float], str] = {
    # A combustor with no oxidizer has nothing to burn with (#892): refused.
    ("fuel_boundary", "palette", 0.0): "no oxidizer",
    ("fuel_boundary", "palette", -0.5): "no oxidizer",
    # Issue to file: chamber_static_state returns NaN at the start point of an
    # orifice in a pipe at zero flow (the bare palette orifice solves).
    ("orifice", "in a pipe", 0.0): "issue to file: chamber_static_state NaN",
    ("orifice", "in a pipe, jobson1955", 0.0): "issue to file: chamber_static_state NaN",
    ("orifice", "in a pipe, Lichtarowicz", 0.0): "issue to file: chamber_static_state NaN",
    ("orifice", "in a pipe, IdelchikThick", 0.0): "issue to file: chamber_static_state NaN",
    # Converges outside pytest only through the outlet-referenced retry, which
    # COMBAERO_SOLVER_STRICT (conftest) turns off: the same NaN.
    ("orifice", "in a pipe, McGreehanSchotsch", 0.0): "issue to file: chamber_static_state NaN",
    # A turbomachine outside its operating range (#893: power > 0, a pressure
    # ratio, no reverse flow).
    ("compressor", "palette", -0.5): "reverse flow through a compressor",
    ("compressor", "palette", 30.0): "pressure ratio beyond the compressor's",
    ("turbine", "palette", 0.0): "a turbine with no flow",
    ("turbine", "palette", -0.5): "reverse flow through a turbine",
    # Impingement arrays do not run backwards; 30x the drive does converge.
    ("impingement_array", "palette", -0.5): "reverse flow through a jet plate",
    ("impingement_array", "bypass from a boundary", -0.5): "reverse flow through a jet plate",
    ("impingement_array", "thermal wall", -0.5): "reverse flow through a jet plate",
    # A Mynard / constant-K tee has a declared flow direction (#629): no
    # reversal, and no determined split at zero flow.
    ("mpce_tee", "palette", 0.0): "tee at zero flow",
    ("mpce_tee", "palette", -0.5): "tee against its declared direction",
    ("mpce_tee", "palette", 30.0): "tee: choked junction",
    ("mpce_tee", "branch, mynard", 0.0): "tee at zero flow",
    ("mpce_tee", "branch, mynard", -0.5): "tee against its declared direction",
    ("mpce_tee", "branch, mynard", 30.0): "tee: choked junction",
    ("mpce_tee", "merge, mynard", 0.0): "tee at zero flow",
    ("mpce_tee", "branch, constant_k", 0.0): "tee at zero flow",
    ("mpce_tee", "branch, constant_k", -0.5): "tee against its declared direction",
    ("mpce_tee", "branch, constant_k", 30.0): "tee: choked junction",
    ("mpce_tee", "merge, constant_k", 0.0): "tee at zero flow",
    ("mpce_tee", "branch, ideal_momentum", 30.0): "tee: choked junction",
    # Reverse flow through an ejector is out of scope (Huang's backflow).
    ("ejector", "palette", -0.5): "ejector backflow",
    ("ejector", "kracik_dvorak", -0.5): "ejector backflow",
    ("ejector", "incompressible (jet pump)", -0.5): "ejector backflow",
    # A Tt_target heat on a node with no through-flow is undetermined: refused.
    ("heat_source", "outlet Tt (total)", 0.0): "Tt_target with no through-flow",
    # An effusion plate with no coolant and no gas flow has no cooling path:
    # its wall T is undetermined, refused (#998).
    ("effusion_plate", "palette", 0.0): "no cooling path",
    ("effusion_plate", "on a gas-side chamber", 0.0): "no cooling path",
    ("effusion_plate", "on a gas-side chamber, baldauf_sellers film", 0.0): "no cooling path",
    ("effusion_liner", "palette", 0.0): "no cooling path",
    ("effusion_liner", "baldauf_sellers film", 0.0): "no cooling path",
}

#: [s] the solver budget at a listed failing point: it only has to show
#: that the point still fails, and its retries ran to the networks' 20 s.
#: The one place a time limit enters here, and only for the listed points.
EXPECTED_FAILURE_TIMEOUT = 3.0

_NETWORKS = list(keys.networks())


def _scaled(graph: dict, f: float) -> dict:
    """The driving pressure difference above the lowest boundary, or the fed
    mass flow, times f."""
    g = copy.deepcopy(graph)
    pts = [
        n["data"]["Pt"]
        for n in g["nodes"]
        if n["type"] == "pressure_boundary" and "Pt" in n["data"]
    ]
    p_lo = min(pts) if pts else None
    for n in g["nodes"]:
        d = n["data"]
        if n["type"] == "pressure_boundary" and p_lo is not None and "Pt" in d:
            d["Pt"] = p_lo + (d["Pt"] - p_lo) * f
        if n["type"] == "mass_boundary" and "m_dot" in d:
            d["m_dot"] = d["m_dot"] * f
    return g


@pytest.mark.slow
@pytest.mark.parametrize(
    ("kind", "variant", "graph"),
    _NETWORKS,
    ids=[f"{k}[{v}]" for k, v, _ in _NETWORKS],
)
def test_the_keys_are_the_same_at_every_operating_point(
    kind: str, variant: str, graph: dict
) -> None:
    sets: dict[float, set[str]] = {}
    for f in FACTORS:
        point = _scaled(graph, f)
        if (kind, variant, f) in EXPECTED_FAILURES:
            point.setdefault("solver_settings", {})["timeout"] = EXPECTED_FAILURE_TIMEOUT
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = http._post("/solve", point)
        if not (r.status_code == 200 and r.json()["success"]):
            assert f != 1.0, f"its own point did not converge: {r.text[:300]}"
            continue
        body = r.json()
        res = body["element_results"].get("x") or body["node_results"].get("x")
        assert res is not None
        sets[f] = keys._numeric_keys(res)
    assert len(sets) >= 2, f"only x{sorted(sets)} converged: nothing to compare"
    ref = sets[1.0]
    for f, got in sets.items():
        assert got == ref, {
            "at": f"x{f}",
            "absent here": sorted(ref - got),
            "only here": sorted(got - ref),
        }
