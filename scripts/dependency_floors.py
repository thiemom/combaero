"""The runtime dependency floors that combaero and combaero-gui declare.

Dependabot keeps the lock file at the newest versions, so every PR tests the
top of each declared range. This script names the bottom: one requirement
line per runtime dependency of both packages (core dependencies, the
`examples` and `generator` extras, and combaero-gui), which the nightly
`dependency-floors` job installs with `--resolution lowest-direct` and tests.

A requirement without a lower bound is an error: lowest-direct would then
resolve it to the first release ever published (fastapi 0.1.0, uvicorn
0.0.1), which tests nothing and, for the GUI's web stack, ships known
vulnerabilities. A floor is the lowest version that (a) has no known
advisory and (b) passes the suite -- see the floor comments in each
pyproject.toml.

    uv run scripts/dependency_floors.py            # print the requirements
    uv run scripts/dependency_floors.py --check    # fail on an unbounded one
    uv run scripts/dependency_floors.py --verify   # installed == floor, each
    uv run scripts/dependency_floors.py --osv      # no floor has an advisory

`--verify` exists because an install that leaves an already-satisfied
newer version in place passes every test while testing nothing (measured:
`uv pip install --resolution lowest-direct` over a synced env kept numpy
2.5.4). The job therefore compiles exact pins first and verifies after.

`--osv` asks OSV (api.osv.dev) for advisories that affect each floor
version. Dependabot moves only the top of each range, so an advisory
published against an old release reaches the floor and nothing else; run
it with every Dependabot round (and nightly, in the floors job), and raise
a flagged floor to the first fixed release that passes the suite.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
import tomllib
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from packaging.requirements import Requirement
from packaging.version import Version

ROOT = Path(__file__).resolve().parent.parent
CORE = ROOT / "pyproject.toml"
GUI = ROOT / "gui" / "pyproject.toml"
CORE_EXTRAS = ("examples", "generator")
OWN = {"combaero", "combaero-gui"}
OSV_BATCH = "https://api.osv.dev/v1/querybatch"


def requirements() -> list[tuple[str, Requirement]]:
    """(source, requirement) for every runtime dependency of both packages."""
    out: list[tuple[str, Requirement]] = []
    core = tomllib.loads(CORE.read_text())["project"]
    out += [("pyproject.toml", Requirement(r)) for r in core["dependencies"]]
    for extra in CORE_EXTRAS:
        out += [(f"pyproject.toml [{extra}]", Requirement(r)) for r in core["optional-dependencies"][extra]]
    gui = tomllib.loads(GUI.read_text())["project"]
    out += [("gui/pyproject.toml", Requirement(r)) for r in gui["dependencies"]]
    return [(src, r) for src, r in out if r.name not in OWN]


def unbounded(reqs: list[tuple[str, Requirement]]) -> list[str]:
    """Requirements that declare no lower bound (>=, ==, ~= or >)."""
    return [
        f"{src}: {r}"
        for src, r in reqs
        if not any(s.operator in (">=", "==", "~=", ">") for s in r.specifier)
    ]


def floor(r: Requirement) -> str | None:
    """The version a lowest-direct resolve installs for `r`."""
    for s in r.specifier:
        if s.operator in (">=", "==", "~="):
            return s.version
    return None


def floors(reqs: list[tuple[str, Requirement]]) -> dict[str, str | None]:
    """The effective floor per package: the highest one any requirement names."""
    return {
        name: max((floor(r) for _, r in reqs if r.name == name and floor(r)), key=Version, default=None)
        for name in dict.fromkeys(r.name for _, r in reqs)
    }


def not_at_floor(reqs: list[tuple[str, Requirement]]) -> list[str]:
    """Installed versions that differ from the declared floor."""
    out = []
    for name, want in floors(reqs).items():
        try:
            have = version(name)
        except PackageNotFoundError:
            have = "not installed"
        if want is not None and have != "not installed" and Version(have) == Version(want):
            continue
        out.append(f"{name}: installed {have}, floor {want}")
    return out


def advisories(reqs: list[tuple[str, Requirement]]) -> dict[str, list[str]]:
    """OSV advisory ids per floor, as "name==version"; [] when clean."""
    pinned = {name: v for name, v in floors(reqs).items() if v is not None}
    queries = [{"package": {"ecosystem": "PyPI", "name": n}, "version": v} for n, v in pinned.items()]
    body = json.dumps({"queries": queries}).encode()
    req = urllib.request.Request(OSV_BATCH, data=body, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        results = json.load(resp)["results"]
    return {
        f"{n}=={v}": sorted(x["id"] for x in res.get("vulns", []))
        for (n, v), res in zip(pinned.items(), results, strict=True)
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="fail if a runtime dependency has no floor")
    ap.add_argument("--verify", action="store_true", help="fail unless each one is installed at its floor")
    ap.add_argument("--osv", action="store_true", help="fail if OSV lists an advisory for any floor version")
    args = ap.parse_args()
    reqs = requirements()
    if args.osv:
        found = advisories(reqs)
        for pin, ids in found.items():
            print(f"{pin:30s} {'ADVISORY ' + ', '.join(ids) if ids else 'ok'}")
        flagged = sum(bool(ids) for ids in found.values())
        print(f"{len(found)} floors, {flagged} with an advisory")
        return 1 if flagged else 0
    if args.verify:
        off = not_at_floor(reqs)
        for o in off:
            print(f"not at floor: {o}")
        print(f"{len(dict.fromkeys(r.name for _, r in reqs))} runtime dependencies, {len(off)} not at their floor")
        return 1 if off else 0
    bad = unbounded(reqs)
    if args.check:
        for b in bad:
            print(f"no floor: {b}")
        print(f"{len(reqs)} runtime requirements, {len(bad)} without a floor")
        return 1 if bad else 0
    if bad:
        print("\n".join(f"no floor: {b}" for b in bad), file=sys.stderr)
        return 1
    # A name listed twice (pandas: core [generator] and the GUI) prints twice;
    # the resolver intersects the two specifiers.
    print("\n".join(str(r) for _, r in reqs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
