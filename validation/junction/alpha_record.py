"""Record: what the joining alpha stands in for (reported, nothing tuned).

Scores the frozen closure set of ``joining_variants.MODELS`` at kernel level
on every joining source the harness has, at the q each network solve reached:

- Bassett 2001 K7/K8/K11/K12 measured (role: selected),
- Idelchik 1966 joining tables (role: tuned -- alpha 0.2 is its optimum),
- Wang 2014 (cross-source; unsealed once for this question, so it is spent
  for any closure-form choice made from these numbers),
- Mynard 2015 Ref3D joining samples (laminar CFD, reported only).

Blaisdell & Manson 1963 has its own pre-registered scorer, bm1963_score.py.

    uv run python -m validation.junction.alpha_record
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

import numpy as np

from validation.junction import mynard_fidelity as mf
from validation.junction.joining_variants import MODELS, tee
from validation.junction.models.mpce_network import MPCENetwork
from validation.junction.network_runner import run_network

NAMES = list(MODELS)
JOIN = {"K7", "K8", "K11", "K12"}

# A row: (label fields..., K_reference, {model: K})
Row = tuple


def _bassett_idelchik_geometry(r) -> tuple[float, float, int, float]:
    """(q_lat, lateral deflection in deg, kernel slot, psi) as mpce_network builds it.

    Bassett's type 4 (K7/K8) has the lateral opposing the outlet, so its axis
    is 180 - theta; K8/K11 are straight-branch coefficients (slot 1), K7/K12
    lateral (slot 2), and the straight-branch K is quoted against 1 - q.
    """
    qf = r.q_converged if r.q_converged is not None else r.q
    th = r.theta_deg if r.theta_deg is not None else 90.0
    psi = r.psi or 1.0
    if r.paper == "bassett2001":
        if r.K_id in ("K7", "K8"):
            return (1 - qf if r.K_id == "K8" else qf), 180.0 - th, (1 if r.K_id == "K8" else 2), psi
        return (1 - qf if r.K_id == "K11" else qf), th, (1 if r.K_id == "K11" else 2), psi
    return qf, th, (1 if r.K_id == "K11" else 2), psi


def _models_at(q: float, psi: float, th_deg: float, slot: int) -> dict[str, float]:
    U, A, t = tee(q, psi, th_deg)
    return {m: f(U, A, t)[slot] for m, f in MODELS.items()}


def network_rows() -> tuple[list[Row], list[Row]]:
    recs = [
        r
        for r in run_network(MPCENetwork(), topologies=("imposed_q",))
        if r.converged and r.K_extracted is not None
    ]
    bi, wang = [], []
    for r in recs:
        if (r.paper == "bassett2001" and r.K_id in JOIN) or r.paper == "idelchik1966":
            if r.K_measured_at_q_converged is None:
                continue
            q, th, slot, psi = _bassett_idelchik_geometry(r)
            bi.append(
                (
                    r.paper,
                    r.K_id,
                    psi,
                    r.theta_deg,
                    r.K_measured_at_q_converged,
                    _models_at(q, psi, th, slot),
                )
            )
        elif r.paper == "wang2014":
            q = r.q_converged if r.q_converged is not None else r.q
            slot = 2 if r.K_id == "K_13" else 1
            interior = 0 < r.q < 1
            wang.append((r.psi, r.mach, interior, r.K_measured, _models_at(q, r.psi, 45.0, slot)))
    return bi, wang


def mynard_rows() -> list[Row]:
    spec = {s.key: s for s in mf.PANELS}
    out = []
    for p in mf.samples():
        s = spec[p.panel]
        if s.flow_type not in (4, 6):
            continue
        kw = {"lam": 0.5, "phi_deg": None, "area": 1.0}
        if s.sweep == "lam":
            kw["lam"] = p.x
        elif s.sweep == "phi":
            kw["phi_deg"] = p.x
        elif s.sweep == "area":
            kw["area"] = p.x
        Q, A, th, slots = mf._junction(s.flow_type, **kw)
        U = Q / A
        ks = {m: f(U, A, th)[slots[s.branch] + 1] for m, f in MODELS.items()}
        out.append((s.flow_type, s.sweep, p.x, p.K_ref3d, ks))
    return out


def line(label: str, sel: Sequence[Row], ref: Callable[[Row], float]) -> None:
    if not sel:
        return
    out = f"  {label:30s} {len(sel):4d}"
    for m in NAMES:
        e = np.array([r[-1][m] - ref(r) for r in sel])
        out += f"  {np.abs(e).mean():6.3f} ({e.mean():+.2f})"
    print(out)


def main() -> None:
    bi, wang = network_rows()
    my = mynard_rows()
    print(
        "MAE (bias), kernel level at the converged q" + "".join(f"  {m[:15]:>15s}" for m in NAMES)
    )
    for paper, role in (
        ("bassett2001", "selected"),
        ("idelchik1966", "tuned: alpha 0.2 is its optimum"),
    ):
        P = [r for r in bi if r[0] == paper]
        print(f"{paper} ({role})")
        ref = lambda r: r[4]  # noqa: E731
        line("psi = 1", [r for r in P if r[2] == 1], ref)
        line("psi > 1", [r for r in P if r[2] > 1], ref)
        for K in sorted({r[1] for r in P}):
            line(f"{K} psi > 1", [r for r in P if r[1] == K and r[2] > 1], ref)
        if paper == "idelchik1966":
            for th in (30, 45, 90):
                line(f"theta {th} psi > 1", [r for r in P if r[3] == th and r[2] > 1], ref)
    print("wang2014 (CROSS-SOURCE; spent for closure-form choices)")
    ref = lambda r: r[3]  # noqa: E731
    line("interior", [r for r in wang if r[2]], ref)
    line("interior psi = 1", [r for r in wang if r[2] and r[0] == 1], ref)
    line("interior psi > 1", [r for r in wang if r[2] and r[0] > 1], ref)
    line("interior psi > 1, M < 0.3", [r for r in wang if r[2] and r[0] > 1 and r[1] < 0.3], ref)
    print("mynard2015 Ref3D joining (laminar CFD, reported only)")
    ref = lambda r: r[3]  # noqa: E731
    line("Fig 9 T", [r for r in my if r[0] == 4], ref)
    line(
        "Fig 9 area sweep, psi != 1",
        [r for r in my if r[0] == 4 and r[1] == "area" and not math.isclose(r[2], 1.0)],
        ref,
    )
    line("Fig 11 Y", [r for r in my if r[0] == 6], ref)


if __name__ == "__main__":
    main()
