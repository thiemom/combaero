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

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

# Two networks x two processes, about 30 s: CI's parallel slow job.
pytestmark = pytest.mark.slow

STRESS = Path(__file__).resolve().parents[2] / "scripts" / "stress"

_PROBE = """
import json, sys, warnings
sys.path.insert(0, {stress!r})
warnings.simplefilter("ignore")
from networks import FAMILIES
from combaero.network import NetworkSolver
net, _ = FAMILIES[{family!r}]({seed})
r = NetworkSolver(net).solve(timeout=120.0)
out = {{k: float(v).hex() for k, v in r.items()
        if isinstance(v, float | int) and not isinstance(v, bool)}}
print(json.dumps(out, sort_keys=True))
"""


def _solve_in_process(family: str, seed: int, hash_seed: int) -> dict[str, str]:
    code = _PROBE.format(stress=str(STRESS), family=family, seed=seed)
    env = {**os.environ, "PYTHONHASHSEED": str(hash_seed)}
    done = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
    )
    return json.loads(done.stdout.strip().splitlines()[-1])


@pytest.mark.parametrize(("family", "seed"), [("random", 0), ("random", 3)])
def test_the_solve_is_bit_identical_across_hash_seeds(family: str, seed: int) -> None:
    a = _solve_in_process(family, seed, 0)
    b = _solve_in_process(family, seed, 1)
    assert a.keys() == b.keys()
    differ = sorted(k for k in a if a[k] != b[k])
    assert not differ, f"{len(differ)} outputs differ, e.g. {differ[:5]}"
