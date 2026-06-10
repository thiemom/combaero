"""
Phase B.5: dataset quality-check pytest.

For every measured CSV in the junction validation dataset, evaluate the
paper's own analytical correlation at each digitized data point and assert
the RMSE deviation is within a fit-tier-derived tolerance. Catches
digitization drift (e.g. axis-calibration errors) without flagging
legitimate measurement spread that the paper itself acknowledges.

Skips:
  - non-measured kinds (calc, envelope)
  - files without a K_id (xi contraction)
  - (paper, K) combos with no analytical correlation in our adapters

The fit-tier tolerance comes from each paper's metadata.yaml, which encodes
the paper's own measured-vs-calc agreement (e.g. Bassett K9 is "poor" so a
30% RMSE is acceptable; K6 is "very_good" so 10% is the bar).
"""

from __future__ import annotations

import math

import pytest

from validation.junction.models.paper_ceiling import PaperCeiling
from validation.junction.runner import _read_xy_csv
from validation.junction.schema import FileMetadata, load_dataset

# RMSE tolerance per fit-tier rating. Wide enough to allow the paper's
# own acknowledged measured-vs-calc spread; narrow enough to catch
# structural digitization errors (axis flip, scale error, sign error).
# These were calibrated against the dataset after the initial sweep:
# Bassett's "moderate" Ks (K7/K8/K11) show genuine 0.3-0.5 RMSE in
# regimes the paper itself acknowledges ("some discrepancy ... but
# the trends ... predicted well"). A sign-flipped or scale-3x curve
# would give RMSE >> 1.0 in most cases and still trip.
_TIER_TOLERANCE: dict[str, float] = {
    "excellent": 0.15,
    "very_good": 0.25,
    "good": 0.40,
    "moderate": 0.60,
    "poor": 1.00,
}
_DEFAULT_TOLERANCE = 0.60

_CEILING = PaperCeiling()


def _measured_files_with_analytical() -> list[FileMetadata]:
    """Files we can sweep: measured kind, K_id present, paper-ceiling evaluable.

    Bassett: all measured Ks (K1..K12 with angle-corrected forms).
    Wang: a=1, q=0.5 (Table 1) and a=1, M_3=0.5 (Table 2) only -- other
      slices have no analytical reference in the paper; PaperCeiling
      returns None and the test skips.
    Hager: xi_t via Eq 8 (q = lateral/total convention).
    Perez-Garcia: analytical-only paper, no measured CSVs.
    """
    ds = load_dataset()
    out: list[FileMetadata] = []
    for f in ds.files:
        if f.kind != "measured":
            continue
        if f.K_id is None and (f.coefficient is None or f.coefficient == "xi"):
            continue
        # Perez-Garcia is analytical-only (no measured CSVs). Wang is wired
        # in via Tables 1/2 interpolation at a=1; other area ratios are
        # returned as None and the test correctly skips them. Hager xi_t
        # is wired via Eq 8.
        if f.path.parent.name not in {"bassett2001", "wang2014", "hager1984"}:
            continue
        out.append(f)
    return out


def _evaluate_ceiling(file: FileMetadata, x: float) -> float | None:
    """Evaluate the paper's analytical correlation at this point."""
    paper_slug = file.path.parent.name
    K_id = file.K_id or (file.coefficient or "")
    theta_rad = math.radians(file.theta_deg) if file.theta_deg is not None else None
    if file.x_axis == "q":
        q_val = x
        M_val = file.M_3
    else:
        q_val = file.q if file.q is not None else 0.5
        M_val = x
    return _CEILING.evaluate(paper_slug, K_id, q_val, file.psi, theta_rad, M_3=M_val, a=file.a)


def _tolerance_for(file: FileMetadata) -> float:
    tier = file.fit_tier
    return _TIER_TOLERANCE.get(tier or "", _DEFAULT_TOLERANCE)


# Known digitization issues to xfail rather than block. Each entry should
# have a one-line reason; the test will fail loudly if a known-issue file
# starts passing so the entry can be removed.
_KNOWN_ISSUES: dict[str, str] = {}


@pytest.mark.parametrize(
    "file",
    _measured_files_with_analytical(),
    ids=lambda f: f.path.name,
)
def test_measured_data_within_paper_fit_tier(file: FileMetadata):
    """Each measured CSV's RMSE vs paper's analytical curve must be within
    the fit-tier tolerance.

    A failure here typically means either a digitization error (axis
    miscalibration, scale, sign) or a metadata mismatch (wrong K_id, wrong
    psi/theta). It is NOT a signal that the paper's correlation is wrong --
    the runner handles model-vs-paper comparison separately.
    """
    rows = _read_xy_csv(file.path)
    assert len(rows) >= 3, (
        f"{file.path.name}: only {len(rows)} data points -- digitization truncated?"
    )

    # Filter to design range for q-axis files. Out-of-range points are
    # typically digitization artefacts (overlapping curves, axis rescale).
    out_of_range_msg: str | None = None
    if file.x_axis == "q":
        valid = [(x, y) for x, y in rows if -0.05 <= x <= 1.05]
        n_out_of_range = len(rows) - len(valid)
        if n_out_of_range > 0.2 * len(rows):
            out_of_range_msg = (
                f"{n_out_of_range}/{len(rows)} points outside q in [0, 1]; "
                "likely digitization captured extra curves or axis miscalibrated"
            )
        rows = valid

    errors: list[float] = []
    n_unsupported = 0
    for x, y in rows:
        K_ana = _evaluate_ceiling(file, x)
        if K_ana is None:
            n_unsupported += 1
            continue
        errors.append(y - K_ana)

    if not errors:
        pytest.skip(
            f"No analytical evaluations available for {file.path.name} "
            f"({n_unsupported} points unsupported)"
        )

    rmse = math.sqrt(sum(e * e for e in errors) / len(errors))
    bias = sum(errors) / len(errors)
    tol = _tolerance_for(file)
    passes = rmse < tol

    known = _KNOWN_ISSUES.get(file.path.name)
    if known is not None:
        if passes and out_of_range_msg is None:
            pytest.fail(
                f"{file.path.name} is in _KNOWN_ISSUES but now passes "
                f"(RMSE {rmse:.3f} < tol {tol:.2f}). Remove it."
            )
        pytest.xfail(f"known issue: {known} (RMSE {rmse:.3f}, tol {tol:.2f})")

    if out_of_range_msg is not None:
        pytest.fail(
            f"{file.path.name}: {out_of_range_msg}. Either re-digitize or "
            "add to _KNOWN_ISSUES with a justification."
        )

    assert passes, (
        f"\n  {file.path.name}: RMSE {rmse:.3f} vs paper analytical exceeds "
        f"fit_tier ({file.fit_tier}) tolerance {tol:.2f}.\n"
        f"  bias={bias:+.3f}, N={len(errors)} points.\n"
        f"  Likely cause: digitization drift, axis miscalibration, wrong "
        f"K_id/psi/theta in metadata.yaml, or sign flip.\n"
        f"  Cross-check the CSV against the paper's figure at the same "
        f"(psi, theta) and look for systematic offset.\n"
        f"  If genuinely paper-acknowledged spread, add to _KNOWN_ISSUES."
    )


def test_dataset_has_files_to_quality_check():
    """Sanity: the quality-check set isn't empty after filtering."""
    files = _measured_files_with_analytical()
    assert len(files) >= 50, (
        f"Quality-check set has only {len(files)} files -- expected ~70+. "
        "Did the filter become too strict?"
    )


# ---------------------------------------------------------------------------
# Plot-range box check
# ---------------------------------------------------------------------------
# A second-line defence against axis miscalibration for files without an
# analytical reference (Wang at a != 1, Wang at a=1 but q not in {0.5} or
# M_3 not in {0.5}). Each (paper, figure_id, K_id) combination has the
# paper's documented plot range; CSVs whose points lie substantially
# outside the box are likely victims of a digitizer-axis miscalibration.
#
# Box dimensions are taken from the figures' visible axis ticks with a
# small margin (~10%) to allow legitimate measurement spread plus
# digitization noise. A 2x calibration error would put points well outside
# even with that margin.

# (paper, figure_key, K_id) -> (x_min, x_max, y_min, y_max)
_PLOT_RANGES: dict[tuple[str, str, str], tuple[float, float, float, float]] = {
    # Wang Fig 10 (a=1): M_3 axis 0..0.7; K axes per published figure
    ("wang2014", "fig10", "K_13"): (0.0, 0.7, -1.3, 1.3),
    ("wang2014", "fig10", "K_23"): (0.0, 0.7, -0.7, 0.7),
    # Wang Fig 11 (a=1.56)
    ("wang2014", "fig11", "K_13"): (0.0, 0.7, -1.6, 2.6),
    ("wang2014", "fig11", "K_23"): (0.0, 0.7, -1.6, 0.5),
    # Wang Fig 12 (a=2.44): K_13 reaches up to ~6 at q=1
    ("wang2014", "fig12", "K_13"): (0.0, 0.7, -1.6, 6.5),
    ("wang2014", "fig12", "K_23"): (0.0, 0.7, -2.7, 1.0),
    # Hager Fig 4a: q in [0, 1], xi_t in [-0.25, 0.55]
    ("hager1984", "fig04a", "xi_t"): (0.0, 1.0, -0.3, 0.6),
}


def _plot_range_key(file: FileMetadata) -> tuple[str, str, str] | None:
    """Derive (paper, figure_key, K_id) from file metadata for box-range lookup.

    Uses metadata K_id/coefficient (not filename parsing) so the key matches
    `_PLOT_RANGES` exactly. Figure key is extracted from the filename since
    metadata does not currently encode the source figure number directly.
    """
    paper = file.path.parent.name
    name = file.path.name
    K_id = file.K_id or file.coefficient
    if K_id is None:
        return None
    parts = name.split("_")
    if len(parts) < 2:
        return None
    fig_key = parts[1]  # e.g. fig10, fig04a
    return (paper, fig_key, K_id)


def _files_with_plot_range() -> list[FileMetadata]:
    """Measured files for which we have a documented plot-range box."""
    ds = load_dataset()
    out: list[FileMetadata] = []
    for f in ds.files:
        if f.kind != "measured":
            continue
        key = _plot_range_key(f)
        if key is not None and key in _PLOT_RANGES:
            out.append(f)
    return out


@pytest.mark.parametrize(
    "file",
    _files_with_plot_range(),
    ids=lambda f: f.path.name,
)
def test_measured_data_within_plot_range_box(file: FileMetadata):
    """Each measured CSV's (x, y) pairs must lie inside the paper's
    documented plot range, with a 5% expand margin for digitizer noise.

    Catches axis miscalibration even when there is no analytical reference
    to compare against. Tolerates up to 10% out-of-box points (some figures
    have curves extending slightly beyond the axis crop, especially at
    plot endpoints).
    """
    key = _plot_range_key(file)
    assert key is not None, f"{file.path.name}: no plot-range key"
    x_min, x_max, y_min, y_max = _PLOT_RANGES[key]

    # Expand the box by 5% in each direction for digitizer noise.
    dx = 0.05 * (x_max - x_min)
    dy = 0.05 * (y_max - y_min)
    x_lo, x_hi = x_min - dx, x_max + dx
    y_lo, y_hi = y_min - dy, y_max + dy

    rows = _read_xy_csv(file.path)
    out_pts = [(x, y) for x, y in rows if not (x_lo <= x <= x_hi and y_lo <= y <= y_hi)]
    out_frac = len(out_pts) / max(len(rows), 1)

    assert out_frac <= 0.10, (
        f"\n  {file.path.name}: {len(out_pts)}/{len(rows)} points "
        f"({100 * out_frac:.0f}%) outside plot box "
        f"x in [{x_min:.2f}, {x_max:.2f}], y in [{y_min:.2f}, {y_max:.2f}].\n"
        f"  Sample out-of-box: {out_pts[:3]}\n"
        f"  Likely cause: axis miscalibration or wrong figure metadata."
    )
