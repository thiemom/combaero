"""The stress harness (scripts/stress) keeps building valid networks.

A dev tool, not shipped: it builds every element family cold and over its
range, which is how #614 and #619-#624 were found. This only guards that its
generators still construct and validate when an element's API changes; the
sweeps themselves are run by hand (scripts/stress/run.py).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))

from networks import FAMILIES  # noqa: E402


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_every_family_builds_valid_networks(family: str) -> None:
    for seed in range(5):
        net, tags = FAMILIES[family](seed)
        net.validate()
        assert net.elements and tags


def test_a_seed_is_reproducible() -> None:
    """A recorded finding (family, seed) must rebuild the same network."""
    for family, build in FAMILIES.items():
        a, b = build(3)[0], build(3)[0]
        assert sorted(a.elements) == sorted(b.elements), family
        assert sorted(a.nodes) == sorted(b.nodes), family


def test_ejector_draws_stay_inside_their_operating_range() -> None:
    """Incompressible ejector draws put the outlet below the jet pump's
    shut-off head (#619), dPt_o = 2 (Pt_p - Pt_s) b (1 - b/2), b = A_e/A_3,
    and every port chamber is in the draw's regime (it used to default to
    compressible)."""
    seen = 0
    for seed in range(40):
        net, tags = FAMILIES["ejector"](seed)
        regime = tags[1]
        assert {net.nodes[n].regime for n in ("mp", "ms", "mo")} == {regime}
        if regime != "incompressible":
            continue
        seen += 1
        ej = net.elements["ej"]
        b = ej.nozzle_exit_area / ej.mixing_area
        Pp, Ps, Po = (net.nodes[n].Pt for n in ("bp", "bs", "bo"))
        assert Ps < Po < Ps + 2.0 * (Pp - Ps) * b * (1.0 - 0.5 * b)
    assert seen > 5
