"""Score the ideal momentum tee (#517) against its pre-registration.

The predictions were written into issue #517 before any of this ran:

- S1, dividing straight run: against Bassett's measured K5 (Fig 7c) the model
  OVER-predicts, the bias growing as the straight fraction q falls: its
  kappa = 1 station gives K5 = (1 - q)^2, Bassett's own fit
  K5 = q^2 - 1.5 q + 0.5 is the kappa = 0.75 one, so the ideal model sits up
  to ~+0.4 above it near q = 0.2.
- S2, dividing branch: against Bassett's measured K6 (Fig 7a/7b), MAE and
  bias reported per cell, size not predicted.
- M1, merge: the network reproduces the scored ``BM theory`` variant at the
  operating point it reached. An implementation check, not new evidence.

Everything is scored at NETWORK level, through ``IdealMomentumTeeElement``
on the harness's imposed-flow networks, i.e. the element production would
run; the shipped Mynard tee is scored on the same networks for context.
Nothing is tuned.

Roles, printed with every block (docs/VALIDATION_POLICY.md): Bassett 2001 is
SELECTED -- it chose the shipped closure's eta and its dividing-streamline
recovery -- so it is evidence about the ideal model but not an independent
test of the shipped one. Idelchik 1966 joining is TUNED (alpha 0.2 is its
optimum). Wang 2014 is cross-source but SPENT for closure-form choices.

    uv run python -m validation.junction.ideal_tee_score
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence

import numpy as np

from validation.junction.joining_variants import MODELS, tee
from validation.junction.models import bassett2001
from validation.junction.models.ideal_network import IdealMomentumNetwork
from validation.junction.models.mpce_network import MPCENetwork
from validation.junction.models._network_builder import _M_DOT_REF
from validation.junction.network_runner import NetworkRecord, run_network

#: The joining coefficients M1 checks, by source.
JOIN = {"K7", "K8", "K11", "K12"}

#: Network K against kernel K, M1, relative to max(|K|, 1). The network's
#: port heads close on the isentropic stagnation state, the kernel's on
#: rho u^2 / 2; they part by ~M^2/4 of the head, so the tolerance holds only
#: where every port is slow. The harness's reference flow puts the common leg
#: at M ~ 0.025 and the branch at q psi times that.
M1_TOLERANCE = 0.014
M1_MACH_LIMIT = 0.15
_M_COMMON_REF = 0.025

def _key(r: NetworkRecord) -> tuple[str, float, float]:
    return (r.K_id, float(r.psi or 1.0), float(r.theta_deg or 90.0))


def _scored(recs: Sequence[NetworkRecord]) -> dict[tuple[str, float, float, float], NetworkRecord]:
    """Records that converged and landed on the digitised curve, by (K, psi, theta, q)."""
    out = {}
    for r in recs:
        if r.converged and r.error is not None:
            out[(*_key(r), round(r.q, 6))] = r
    return out


def _bassett_fit(K_id: str, q: float, psi: float, theta_deg: float) -> float:
    """Bassett's own calculated coefficient, on the file's q axis."""
    if K_id == "K5":
        return bassett2001.K5(q)
    return bassett2001.K6(q, psi, math.radians(theta_deg))


def _stats(e: Sequence[float]) -> str:
    a = np.asarray(e, dtype=float)
    return f"{np.abs(a).mean():6.3f} ({a.mean():+.3f})" if a.size else "     -        "


def dividing_report(
    ideal: Sequence[NetworkRecord], shipped: Sequence[NetworkRecord]
) -> dict[str, dict[str, float]]:
    """S1 and S2: per measured slice, MAE (bias) of ideal, shipped, Bassett's fit.

    Returns, per coefficient, the all-slice ``n`` and the bias of each model.
    """
    I, S = _scored(ideal), _scored(shipped)
    cells: dict[tuple[str, float, float], list[tuple]] = defaultdict(list)
    for k, r in I.items():
        if r.paper != "bassett2001" or r.K_id not in ("K5", "K6"):
            continue
        s = S.get(k)
        fit = _bassett_fit(r.K_id, r.q, *_key(r)[1:])
        cells[_key(r)].append(
            (r.q, r.K_measured_at_q_converged, r.K_extracted, s.K_extracted if s else None, fit)
        )
    head = f"  {'cell':24s} {'n':>3s}  {'ideal':>15s}  {'shipped':>15s}  {'Bassett fit':>15s}"
    summary: dict[str, dict[str, float]] = {}
    for K_id, title in (
        ("K5", "S1  K5, dividing straight run (Fig 7c), q = straight fraction"),
        ("K6", "S2  K6, dividing branch (Fig 7a/7b), q = branch fraction"),
    ):
        print(f"{title}\n    Bassett 2001 measured -- role: SELECTED")
        print(head)
        all_rows = []
        for key in sorted(k for k in cells if k[0] == K_id):
            rows = cells[key]
            all_rows += rows
            print(
                f"  psi {key[1]:g}, theta {key[2]:>3g}{'':10s} {len(rows):3d}  "
                f"{_stats([i - m for _, m, i, _, _ in rows])}  "
                f"{_stats([s - m for _, m, _, s, _ in rows if s is not None])}  "
                f"{_stats([f - m for _, m, _, _, f in rows])}"
            )
        print(
            f"  {'all':24s} {len(all_rows):3d}  "
            f"{_stats([i - m for _, m, i, _, _ in all_rows])}  "
            f"{_stats([s - m for _, m, _, s, _ in all_rows if s is not None])}  "
            f"{_stats([f - m for _, m, _, _, f in all_rows])}"
        )
        if K_id == "K5":
            # The pre-registered sign and location: ideal minus Bassett's fit.
            gap = [(q, i - f) for q, _, i, _, f in all_rows]
            q_max, g_max = max(gap, key=lambda t: t[1])
            low = [m_err for q, m_err in ((q, i - m) for q, m, i, _, _ in all_rows) if q < 0.5]
            high = [m_err for q, m_err in ((q, i - m) for q, m, i, _, _ in all_rows) if q >= 0.5]
            print(
                f"  pre-registered: ideal - Bassett fit peaks at {g_max:+.3f} (q = {q_max:.2f}); "
                f"bias vs measured {np.mean(low):+.3f} for q < 0.5, {np.mean(high):+.3f} for q >= 0.5"
            )
        summary[K_id] = {
            "n": float(len(all_rows)),
            "bias_ideal": float(np.mean([i - m for _, m, i, _, _ in all_rows])),
            "bias_shipped": float(np.mean([s - m for _, m, _, s, _ in all_rows if s is not None])),
        }
        print()
    return summary


def merge_check(ideal: Sequence[NetworkRecord]) -> float:
    """M1: network K against the kernel ``BM theory`` K at the reached q.

    Returns the worst relative difference among records whose branch stays
    below ``M1_MACH_LIMIT``. Above it the difference is compressibility, and
    the report demonstrates that rather than asserting it: the worst
    fast-branch record is re-solved at a tenth of the flow, where an M^2
    effect must fall ~100x.
    """
    bm = MODELS["BM theory"]
    rows = []
    for r in ideal:
        if not (r.converged and r.K_extracted is not None):
            continue
        if r.paper == "bassett2001" and r.K_id in JOIN:
            q = r.q_converged if r.q_converged is not None else r.q
            th = r.theta_deg if r.theta_deg is not None else 90.0
            if r.K_id in ("K7", "K8"):
                th = 180.0 - th
            straight = r.K_id in ("K8", "K11")
            q_lat = 1.0 - q if straight else q
        elif r.paper == "idelchik1966" and r.K_id in ("K11", "K12"):
            q_lat = r.q_converged if r.q_converged is not None else r.q
            th = r.theta_deg if r.theta_deg is not None else 90.0
            straight = r.K_id == "K11"
        else:
            continue
        if not 0.0 < q_lat < 1.0:
            continue
        psi = float(r.psi or 1.0)
        U, A, t = tee(q_lat, psi, th)
        k_ref = bm(U, A, t)[1 if straight else 2]
        rel = (r.K_extracted - k_ref) / max(abs(k_ref), 1.0)
        mach_branch = q_lat * psi * _M_COMMON_REF
        rows.append((r.paper, rel, mach_branch, q_lat, psi, th, straight))
    print("M1  merge: network K vs kernel 'BM theory' K at the reached q (implementation check)")
    print(f"  {'source':14s} {'branch Mach':>12s} {'n':>5s}  {'max |rel diff|':>15s}")
    for paper in sorted({r[0] for r in rows}):
        for label, sel in (
            (f"< {M1_MACH_LIMIT:g}", lambda m: m < M1_MACH_LIMIT),
            (f">= {M1_MACH_LIMIT:g}", lambda m: m >= M1_MACH_LIMIT),
        ):
            d = [abs(r[1]) for r in rows if r[0] == paper and sel(r[2])]
            if d:
                print(f"  {paper:14s} {label:>12s} {len(d):5d}  {max(d):15.2e}")
    slow = [abs(r[1]) for r in rows if r[2] < M1_MACH_LIMIT]
    worst = max(slow) if slow else 0.0
    verdict = "PASS" if worst <= M1_TOLERANCE else "FAIL"
    print(f"  tolerance {M1_TOLERANCE:g} below branch Mach {M1_MACH_LIMIT:g}: {verdict} (worst {worst:.2e})")
    fast = [r for r in rows if r[2] >= M1_MACH_LIMIT]
    if fast:
        _, rel, mach, q_lat, psi, th, straight = max(fast, key=lambda r: abs(r[1]))
        net = IdealMomentumNetwork()
        res = net._joining(q_lat, psi, math.radians(th), "imposed_q", m_in=_M_DOT_REF / 10.0)
        k_slow = res.K_straight if straight else res.K_lateral
        U, A, t = tee(q_lat, psi, th)
        k_ref = bm(U, A, t)[1 if straight else 2]
        rel_slow = (k_slow - k_ref) / max(abs(k_ref), 1.0)
        print(
            f"  worst fast-branch record (psi {psi:g}, theta {th:g}, q {q_lat:.3f}, branch M "
            f"{mach:.2f}): {rel:+.2e}; at a tenth of the flow {rel_slow:+.2e} "
            f"(ratio {rel / rel_slow:.0f}, ~100 for an M^2 effect)"
        )
    return worst


def main() -> None:
    ideal = run_network(IdealMomentumNetwork(), topologies=("imposed_q",))
    shipped = run_network(MPCENetwork(), topologies=("imposed_q",))
    n_i = sum(r.converged for r in ideal)
    n_s = sum(r.converged for r in shipped)
    print(
        f"imposed-flow networks: ideal converged {n_i}/{len(ideal)}, "
        f"shipped {n_s}/{len(shipped)}\n"
    )
    dividing_report(ideal, shipped)
    merge_check(ideal)


if __name__ == "__main__":
    main()
