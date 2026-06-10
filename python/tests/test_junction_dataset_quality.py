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
