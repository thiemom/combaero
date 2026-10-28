"""A solve does not depend on PYTHONHASHSEED.

Python randomises str hashing per process, so iterating a set of node or
unknown names follows a different order in every process. Three near-zero-flow
blends (orifice reverse-flow, port discharge, incompressible channel reversal)
built their Jacobian dicts from a set union; the key order then set how
assemble() summed coinciding terms, a 1-ulp change that moved the Newton
path. The stress networks random:0 and random:3 each draw such an element:
the same network solved bit-differently from one process to the next (same
root, within tolerance). Same class as the #289 BFS seed (solver.py).
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

# Four short subprocesses (import + setup + one Jacobian each).
pytestmark = pytest.mark.slow

STRESS = Path(__file__).resolve().parents[2] / "scripts" / "stress"

_PROBE = """
import hashlib, os, sys, warnings
import numpy as np
sys.path.insert(0, {stress!r})
warnings.simplefilter("ignore")
from networks import FAMILIES
from combaero.network import NetworkSolver

orig = NetworkSolver._residuals_and_jacobian
def first_jacobian(self, x, *a, **k):
    r = orig(self, x, *a, **k)
    if r[1] is not None:
        J = r[1].toarray() if hasattr(r[1], "toarray") else np.asarray(r[1])
        h = hashlib.sha1(np.ascontiguousarray(x).tobytes() + J.tobytes()).hexdigest()
        print(h, J.shape[0], int(np.count_nonzero(J)), flush=True)
        os._exit(0)  # not an exception: the solver's safety net would catch it
    return r
NetworkSolver._residuals_and_jacobian = first_jacobian
net, _ = FAMILIES[{family!r}]({seed})
NetworkSolver(net).solve(timeout=120.0)
"""


def _first_jacobian_in_process(family: str, seed: int, hash_seed: int) -> str:
    """Hash of (x, J) at the solve's first Jacobian, in a fresh process.

    The first Jacobian, not the whole solve: the solver's stall detector and
    retry budgets are wall-clock, so on a loaded runner two processes can
    end a long solve differently for reasons unrelated to hash order. The
    first Jacobian is where the hash-order dependence entered (main differed
    there by 1 ulp in 4 of 70 nonzeros for random:3), and it is reached
    before any budget can act.
    """
    code = _PROBE.format(stress=str(STRESS), family=family, seed=seed)
    env = {**os.environ, "PYTHONHASHSEED": str(hash_seed)}
    done = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
    )
    return done.stdout.strip().splitlines()[-1]


@pytest.mark.parametrize(("family", "seed"), [("random", 0), ("random", 3)])
def test_the_first_jacobian_is_bit_identical_across_hash_seeds(family: str, seed: int) -> None:
    a = _first_jacobian_in_process(family, seed, 0)
    b = _first_jacobian_in_process(family, seed, 1)
    assert a == b, f"(x, J) hash differs: {a} vs {b}"
