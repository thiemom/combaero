"""Scan a solver-facing function for the hazards that stall Newton.

Blends, clamps, floors and validity limits are how physical correlations get
made usable, and each one is a chance to put a flat spot or a corner into a
residual. Reading the code does not reliably find them: the expansion-factor
work (#382) had already smoothed two bounds deliberately and still shipped a
hard ``min(S, 1)`` that put a kink at the no-flow boundary. A scan found it.

Three hazards, distinguished by how the derivative behaves as the sampling
interval shrinks -- which is what separates them from each other and from an
honest steep gradient:

  FLOOR       the derivative is exactly zero over a run of points. A Newton
              step there sees no sensitivity at all.
  KINK        the derivative is discontinuous. The change across the point
              does not shrink when the interval does.
  DIVERGENCE  the derivative is unbounded. The change across the point GROWS
              as the interval shrinks.

A smooth point's derivative change falls off with the interval; that is the
null hypothesis and everything else is reported.

Deliberate saturation is not a defect and is not reported as one: a derivative
that decays smoothly towards zero without reaching it is what a physically
saturating term should do. Only an exactly-zero run counts as a FLOOR by
default. Callers who need "not merely non-zero but with room in it" pass
``floor_atol``.

Usage in a test::

    from validation.solver_smoothness import assert_smooth
    assert_smooth(lambda x: cb.some_correlation(x), 0.1, 10.0, log=True)

Usage as a report::

    python -m validation.solver_smoothness
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

__all__ = [
    "Hazard",
    "Finding",
    "scan",
    "assert_smooth",
    "render",
]


class Hazard(Enum):
    FLOOR = "floor"
    KINK = "kink"
    DIVERGENCE = "divergence"


@dataclass(frozen=True)
class Finding:
    hazard: Hazard
    x: float
    detail: str
    #: For a FLOOR, the end of the flat run; otherwise the same as ``x``.
    x_end: float | None = None

    def __str__(self) -> str:
        where = (
            f"[{self.x:.6g}, {self.x_end:.6g}]"
            if self.x_end is not None and self.x_end != self.x
            else f"{self.x:.6g}"
        )
        return f"{self.hazard.value.upper():<10} at {where}: {self.detail}"


def _emit_floor(
    findings: "list[Finding]",
    start: float,
    end: float,
    span_t: float,
    floor_atol: float,
    min_floor_frac: float,
    pos: "Callable[[float], float]",
) -> None:
    """Record a flat run, unless it is too narrow to be anything but an
    isolated stationary point.

    ``pos`` maps to the coordinate the scan walks, so widths on a log sweep
    are judged in log space.
    """
    width = (pos(end) - pos(start)) / span_t
    if width < min_floor_frac:
        return
    findings.append(
        Finding(
            Hazard.FLOOR,
            start,
            f"|f'| <= {floor_atol:g} over {width * 100:.1f}% of the range",
            end,
        )
    )


def _sample_points(lo: float, hi: float, n: int, log: bool) -> list[float]:
    if log:
        if lo <= 0.0:
            raise ValueError("log scan needs lo > 0")
        return [lo * (hi / lo) ** (i / n) for i in range(n + 1)]
    return [lo + (hi - lo) * i / n for i in range(n + 1)]


def scan(
    f: Callable[[float], float],
    lo: float,
    hi: float,
    *,
    n: int = 2000,
    log: bool = False,
    floor_atol: float = 0.0,
    min_floor_frac: float = 0.01,
    probe_rel: float = 2.0e-3,
    deriv_rel: float = 1.0e-5,
    magnitude_rel: float = 1.0e-3,
    merge_rel: float = 5.0e-3,
) -> list[Finding]:
    """Scan ``f`` over ``[lo, hi]`` and report smoothness hazards.

    ``floor_atol`` -- a derivative at or below this counts as floored.
    Defaults to 0.0, i.e. only an exactly-zero derivative is a FLOOR, so that
    a term saturating smoothly towards zero is not reported. Pass a small
    positive value to require headroom instead.

    ``min_floor_frac`` -- a flat run must span at least this fraction of the
    range to be reported. On a ``log`` scan the fraction is measured in log
    space, so a flat decade of a three-decade sweep counts as a third of the
    range rather than the ~0.1% it occupies linearly. This excludes isolated stationary points, which are
    not hazards: ``x**2`` has a zero derivative at exactly one point and is
    perfectly well behaved. Note that a region where the function saturates
    into floating point -- ``1 - exp(-8x)`` at large ``x``, say -- IS reported,
    and correctly so: ``f`` is exactly constant there and a solver really does
    see no sensitivity.

    ``probe_rel`` -- half-width, relative to the scan span, used to look
    across a candidate point. ``deriv_rel`` sets the central-difference step.
    ``magnitude_rel`` gates reporting against the largest derivative seen, so
    numerical noise in a flat region is not reported as a corner.
    """
    if hi <= lo:
        raise ValueError("hi must exceed lo")
    xs = _sample_points(lo, hi, n, log)
    span = hi - lo

    # Widths are judged in the coordinate the scan actually walks. On a log
    # sweep a flat decade is a large part of the range even though it is a
    # sliver of the linear span -- measuring linearly hid exactly that, the
    # Re floor below 1e4 across a 1e3..1e6 sweep.
    def pos(x: float) -> float:
        return math.log(x) if log else x

    span_t = pos(hi) - pos(lo)

    def deriv(x: float, h: float | None = None) -> float:
        if h is None:
            h = max(abs(x) * deriv_rel, span * deriv_rel, 1e-12)
        return (f(x + h) - f(x - h)) / (2.0 * h)

    derivs = [deriv(x) for x in xs]
    peak = max((abs(d) for d in derivs), default=0.0)
    if peak == 0.0:
        return [
            Finding(Hazard.FLOOR, lo, "derivative is zero across the whole range", hi)
        ]

    findings: list[Finding] = []

    # FLOOR: contiguous runs at or below the tolerance.
    run_start: float | None = None
    prev = lo
    for x, d in zip(xs, derivs, strict=True):
        if abs(d) <= floor_atol:
            if run_start is None:
                run_start = x
        else:
            if run_start is not None:
                _emit_floor(findings, run_start, prev, span_t, floor_atol,
                            min_floor_frac, pos)
                run_start = None
        prev = x
    if run_start is not None:
        _emit_floor(findings, run_start, hi, span_t, floor_atol,
                    min_floor_frac, pos)

    # KINK / DIVERGENCE: how the derivative change scales with the interval.
    raw: list[Finding] = []
    for x in xs:
        # The probe straddles x, and on a log sweep it has to be
        # MULTIPLICATIVE. An additive width with a span-based floor makes the
        # probe enormous relative to x at the small end, which silently
        # defeated the check there -- caught because sqrt(x), a reference
        # divergence, came back clean on a log scan.
        if log:
            k = 1.0 + probe_rel * 20.0
            lo_w, hi_w = x / k, x * k
            lo_n, hi_n = x / k ** 0.25, x * k ** 0.25
        else:
            wide = max(abs(x) * probe_rel, span * probe_rel)
            # Shrink near a boundary rather than skipping. A hazard sitting ON
            # a domain edge is the kind that matters most.
            wide = min(wide, x - lo, hi - x)
            lo_w, hi_w = x - wide, x + wide
            lo_n, hi_n = x - wide / 4.0, x + wide / 4.0
        step = max(abs(x) * deriv_rel, span * deriv_rel, 1e-12)
        if (hi_w - lo_w) <= 4.0 * step or lo_w < lo or hi_w > hi:
            continue
        j_wide = abs(deriv(hi_w) - deriv(lo_w))
        j_narrow = abs(deriv(hi_n) - deriv(lo_n))
        if j_wide <= magnitude_rel * peak:
            continue
        ratio = j_narrow / j_wide
        # Smooth: ~0.25 (the interval was quartered). Kink: ~1. Divergent: >1.
        if ratio > 1.5:
            raw.append(
                Finding(
                    Hazard.DIVERGENCE,
                    x,
                    f"derivative change grows {j_wide:.4g} -> {j_narrow:.4g} "
                    f"when the interval is quartered (smooth would fall ~4x)",
                )
            )
        elif ratio > 0.5:
            raw.append(
                Finding(
                    Hazard.KINK,
                    x,
                    f"derivative change holds {j_wide:.4g} -> {j_narrow:.4g} "
                    f"when the interval is quartered (smooth would fall ~4x)",
                )
            )

    # Endpoints get their own check. An interior probe cannot straddle them,
    # and a hazard sitting exactly ON a domain edge is the kind that matters
    # most -- the crossflow term's derivative is unbounded at U1/Vi = 0, which
    # is both the low end of any sweep and the default operating point. Walk
    # inwards geometrically: if |f'| keeps growing as the distance shrinks,
    # the derivative is unbounded at the edge.
    #
    # "Still changing" is not "unbounded", and the difference is the SHAPE of
    # the growth. A power-law blow-up multiplies |f'| by a constant factor
    # each time the distance is quartered; a function whose derivative merely
    # converges to a finite limit shows ratios that decay towards 1. Testing
    # growth alone flagged log(x) at x = 1, where f' is simply settling to 1.
    for edge, inward in ((lo, +1.0), (hi, -1.0)):
        d = span * probe_rel
        mags = []
        for _ in range(7):
            x = edge + inward * d
            # The difference step has to shrink with the offset, or once the
            # offset drops below a fixed step the stencil straddles the edge
            # and the estimate is meaningless.
            if lo <= x <= hi:
                mags.append(abs(deriv(x, h=d * 1.0e-2)))
            d /= 4.0
        if len(mags) < 5:
            continue
        ratios = [b / a for a, b in zip(mags, mags[1:], strict=False) if a > 0.0]
        if len(ratios) < 4:
            continue
        sustained = ratios[-1] >= 1.3 and ratios[-1] / max(ratios[0], 1e-30) >= 0.6
        if sustained and mags[-1] > 4.0 * peak * magnitude_rel:
            raw.append(
                Finding(
                    Hazard.DIVERGENCE,
                    edge,
                    f"|f'| grows without bound approaching the edge: "
                    f"{mags[0]:.4g} -> {mags[-1]:.4g} as the distance falls "
                    f"{4 ** (len(mags) - 1)}x",
                )
            )

    # One report per location, keeping the most severe.
    for cand in raw:
        near = [
            g
            for g in findings
            if g.hazard in (Hazard.KINK, Hazard.DIVERGENCE)
            and abs(g.x - cand.x) < span * merge_rel
        ]
        if near:
            if cand.hazard is Hazard.DIVERGENCE and all(
                g.hazard is Hazard.KINK for g in near
            ):
                for g in near:
                    findings.remove(g)
                findings.append(cand)
            continue
        findings.append(cand)

    return sorted(findings, key=lambda g: (g.x, g.hazard.value))


def render(findings: list[Finding], label: str = "") -> str:
    head = f"{label}: " if label else ""
    if not findings:
        return f"{head}no floors, kinks or divergences"
    return head + "\n".join("  " + str(g) for g in findings)


def assert_smooth(
    f: Callable[[float], float],
    lo: float,
    hi: float,
    *,
    label: str = "function",
    allow: "set[Hazard] | None" = None,
    **kwargs,
) -> None:
    """Fail unless ``f`` is free of smoothness hazards over ``[lo, hi]``.

    ``allow`` lists hazards that are known and accepted, so a deliberate,
    documented limitation does not have to be silenced by dropping the check
    altogether -- the rest of the range stays guarded.
    """
    findings = [g for g in scan(f, lo, hi, **kwargs) if not (allow and g.hazard in allow)]
    if findings:
        raise AssertionError(
            f"{label} is not solver-smooth over [{lo:g}, {hi:g}]:\n"
            + "\n".join("  " + str(g) for g in findings)
        )


def main() -> None:  # pragma: no cover - reporting entry point
    import combaero as cb

    print("Smoothness scan of the McGreehan-Schotsch orifice surface\n")

    checks = [
        (
            "Y vs Cd (eps = 0, the paper's hard clamp)",
            lambda c: cb.mcgreehan_schotsch_1988_expansion_factor(c, 0.7, 1.4, 0.0),
            0.60,
            1.05,
            {},
        ),
        (
            "Y vs Cd (default smoothing)",
            lambda c: cb.mcgreehan_schotsch_1988_expansion_factor(c, 0.7, 1.4),
            0.60,
            1.05,
            {},
        ),
        (
            "Y vs pressure ratio",
            lambda s: cb.mcgreehan_schotsch_1988_expansion_factor(0.90, s, 1.4),
            0.05,
            1.20,
            {},
        ),
        (
            "Cd vs Re",
            lambda r: cb.mcgreehan_schotsch_1988_cd(r, 0.0, 1.0, 0.0),
            1.0e3,
            1.0e6,
            {"log": True},
        ),
        (
            "Cd vs U1/Vi",
            lambda u: cb.mcgreehan_schotsch_1988_cd(3.2e4, 0.0, 1.0, u),
            0.0,
            3.0,
            {},
        ),
    ]
    for label, fn, lo, hi, kw in checks:
        print(render(scan(fn, lo, hi, **kw), label))
        print()


if __name__ == "__main__":  # pragma: no cover
    main()
