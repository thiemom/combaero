"""
One-shot migration: docs/bassett/ + docs/junction/ CSVs ->
validation/junction/data/<paper>/ with uniform naming.

Naming convention:
    paper_figXX[subfig]_Kid[_param1=val1][_param2=val2]_<measured|calc>[_extra].csv

Run from repo root:
    uv run python scripts/migrate_junction_validation_data.py

The script is idempotent: rerunning after a partial run will skip files
that have already been moved.

After verifying the moves, this script can be deleted or kept as a record.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DST_ROOT = REPO_ROOT / "validation" / "junction" / "data"


# ---------------------------------------------------------------------------
# Bassett 2001 (some in docs/bassett, some in docs/junction)
# ---------------------------------------------------------------------------

BASSETT_MAPPINGS: dict[str, str] = {
    # Fig 5: contraction coefficient xi vs q (analytical Eq 22), theta sweep
    "docs/bassett/bassett_fig5_zeta_q_theta=0.csv":   "bassett2001/bassett_fig05_xi_theta=0_calc.csv",
    "docs/bassett/bassett_fig5_zeta_q_theta=30.csv":  "bassett2001/bassett_fig05_xi_theta=30_calc.csv",
    "docs/bassett/bassett_fig5_zeta_q_theta=45.csv":  "bassett2001/bassett_fig05_xi_theta=45_calc.csv",
    "docs/bassett/bassett_fig5_zeta_q_theta=60.csv":  "bassett2001/bassett_fig05_xi_theta=60_calc.csv",
    "docs/bassett/bassett_fig5_zeta_q_theta=90.csv":  "bassett2001/bassett_fig05_xi_theta=90_calc.csv",
    "docs/bassett/bassett_fig5_zeta_q_theta=120.csv": "bassett2001/bassett_fig05_xi_theta=120_calc.csv",
    "docs/bassett/bassett_fig5_zeta_q_theta=180.csv": "bassett2001/bassett_fig05_xi_theta=180_calc.csv",

    # Fig 7a: K6 measured, equal-area (psi=1), theta sweep
    "docs/junction/bassett_fig7a_eq_area_k6_q_theta_45.csv":  "bassett2001/bassett_fig07a_K6_theta=45_psi=1_measured.csv",
    "docs/junction/bassett_fig7a_eq_area_k6_q_theta_60.csv":  "bassett2001/bassett_fig07a_K6_theta=60_psi=1_measured.csv",
    "docs/junction/bassett_fig7a_eq_area_k6_q_theta_90.csv":  "bassett2001/bassett_fig07a_K6_theta=90_psi=1_measured.csv",
    "docs/junction/bassett_fig7a_eq_area_k6_q_theta_120.csv": "bassett2001/bassett_fig07a_K6_theta=120_psi=1_measured.csv",

    # Fig 7b: K6 measured at theta=45, psi sweep
    "docs/junction/bassett_fig7b_theta_45_k6_q_area_ratio_psi_1.csv": "bassett2001/bassett_fig07b_K6_theta=45_psi=1_measured.csv",
    "docs/junction/bassett_fig7b_theta_45_k6_q_area_ratio_psi_3.csv": "bassett2001/bassett_fig07b_K6_theta=45_psi=3_measured.csv",

    # Fig 7c: K5 measured (theta- and psi-independent by Eq 15)
    "docs/junction/bassett_fig7c_k5_q_theta_45_psi_1.csv": "bassett2001/bassett_fig07c_K5_theta=45_psi=1_measured.csv",
    "docs/junction/bassett_fig7c_k5_q_theta_45_psi_3.csv": "bassett2001/bassett_fig07c_K5_theta=45_psi=3_measured.csv",
    "docs/junction/bassett_fig7c_k5_q_theta_60_psi_1.csv": "bassett2001/bassett_fig07c_K5_theta=60_psi=1_measured.csv",
    "docs/junction/bassett_fig7c_k5_q_theta_90_psi_1.csv": "bassett2001/bassett_fig07c_K5_theta=90_psi=1_measured.csv",

    # Fig 8: K4 measured, equal-area, theta sweep
    "docs/junction/bassett_fig8_k4_q_triangles_theta_45.csv": "bassett2001/bassett_fig08_K4_theta=45_psi=1_measured.csv",
    "docs/junction/bassett_fig8_k4_q_squares_theta_90.csv":   "bassett2001/bassett_fig08_K4_theta=90_psi=1_measured.csv",
    "docs/junction/bassett_fig8_k4_q_plus_theta_135.csv":     "bassett2001/bassett_fig08_K4_theta=135_psi=1_measured.csv",

    # Fig 10a: K12 raw (no angle correction), equal-area, theta sweep
    "docs/bassett/bassett_fig10a_k12_q_30deg_x.csv":         "bassett2001/bassett_fig10a_K12_theta=30_psi=1_measured.csv",
    "docs/bassett/bassett_fig10a_k12_q_45deg_triangles.csv": "bassett2001/bassett_fig10a_K12_theta=45_psi=1_measured.csv",
    "docs/bassett/bassett_fig10a_k12_q_60deg_diamonds.csv":  "bassett2001/bassett_fig10a_K12_theta=60_psi=1_measured.csv",
    "docs/bassett/bassett_fig10a_k12_q_90deg_squares.csv":   "bassett2001/bassett_fig10a_K12_theta=90_psi=1_measured.csv",
    "docs/bassett/bassett_fig10a_k12_q_90deg_model.csv":     "bassett2001/bassett_fig10a_K12_theta=90_psi=1_calc.csv",

    # Fig 10b: K12 with angle correction (Eq 33), equal-area, theta sweep
    "docs/bassett/bassett_fig10b_k12_q_30deg_x.csv":         "bassett2001/bassett_fig10b_K12_theta=30_psi=1_measured.csv",
    "docs/bassett/bassett_fig10b_k12_q_45deg_triangles.csv": "bassett2001/bassett_fig10b_K12_theta=45_psi=1_measured.csv",
    "docs/bassett/bassett_fig10b_k12_q_60deg_diamonds.csv":  "bassett2001/bassett_fig10b_K12_theta=60_psi=1_measured.csv",
    "docs/bassett/bassett_fig10b_k12_q_90deg_squares.csv":   "bassett2001/bassett_fig10b_K12_theta=90_psi=1_measured.csv",
    "docs/bassett/bassett_fig10b_k12_q_90deg_model.csv":     "bassett2001/bassett_fig10b_K12_theta=90_psi=1_calc.csv",

    # Fig 10c: K12 at theta=45, psi sweep with corresponding model lines
    "docs/bassett/bassett_fig10c_k12_q_psi=1_triangles.csv":  "bassett2001/bassett_fig10c_K12_theta=45_psi=1_measured.csv",
    "docs/bassett/bassett_fig10c_k12_q_psi=2_squares.csv":    "bassett2001/bassett_fig10c_K12_theta=45_psi=2_measured.csv",
    "docs/bassett/bassett_fig10c_k12_q_psi=3_x.csv":          "bassett2001/bassett_fig10c_K12_theta=45_psi=3_measured.csv",
    "docs/bassett/bassett_fig10c_k12_q_psi=4_diamonds.csv":   "bassett2001/bassett_fig10c_K12_theta=45_psi=4_measured.csv",
    "docs/bassett/bassett_fig10c_k12_q_psi=1_model_line.csv": "bassett2001/bassett_fig10c_K12_theta=45_psi=1_calc.csv",
    "docs/bassett/bassett_fig10c_k12_q_psi=3_model_line.csv": "bassett2001/bassett_fig10c_K12_theta=45_psi=3_calc.csv",
    "docs/bassett/bassett_fig10c_k12_q_psi=4_model_line.csv": "bassett2001/bassett_fig10c_K12_theta=45_psi=4_calc.csv",

    # Fig 11b: K11 at theta=45, psi sweep (measured only)
    "docs/bassett/bassett_fig11b_k11_q_psi=1.csv": "bassett2001/bassett_fig11b_K11_theta=45_psi=1_measured.csv",
    "docs/bassett/bassett_fig11b_k11_q_psi=2.csv": "bassett2001/bassett_fig11b_K11_theta=45_psi=2_measured.csv",
    "docs/bassett/bassett_fig11b_k11_q_psi=3.csv": "bassett2001/bassett_fig11b_K11_theta=45_psi=3_measured.csv",
    "docs/bassett/bassett_fig11b_k11_q_psi=4.csv": "bassett2001/bassett_fig11b_K11_theta=45_psi=4_measured.csv",

    # Fig 12a: K8 measured, equal-area (psi=1), theta sweep
    "docs/bassett/bassett_fig12a_k8_q_theta30.csv": "bassett2001/bassett_fig12a_K8_theta=30_psi=1_measured.csv",
    "docs/bassett/bassett_fig12a_k8_q_theta45.csv": "bassett2001/bassett_fig12a_K8_theta=45_psi=1_measured.csv",
    "docs/bassett/bassett_fig12a_k8_q_theta60.csv": "bassett2001/bassett_fig12a_K8_theta=60_psi=1_measured.csv",
    "docs/bassett/bassett_fig12a_k8_q_theta90.csv": "bassett2001/bassett_fig12a_K8_theta=90_psi=1_measured.csv",

    # Fig 12b: K8 at theta=?, psi sweep (psi=3 absent in paper)
    "docs/bassett/bassett_fig12b_k8_q_psi=1.csv": "bassett2001/bassett_fig12b_K8_psi=1_measured.csv",
    "docs/bassett/bassett_fig12b_k8_q_psi=2.csv": "bassett2001/bassett_fig12b_K8_psi=2_measured.csv",
    "docs/bassett/bassett_fig12b_k8_q_psi=4.csv": "bassett2001/bassett_fig12b_K8_psi=4_measured.csv",

    # Fig 13b: K7 at theta=?, psi sweep
    "docs/bassett/bassett_fig13b_k7_q_psi=1.csv": "bassett2001/bassett_fig13b_K7_psi=1_measured.csv",
    "docs/bassett/bassett_fig13b_k7_q_psi=2.csv": "bassett2001/bassett_fig13b_K7_psi=2_measured.csv",
    "docs/bassett/bassett_fig13b_k7_q_psi=4.csv": "bassett2001/bassett_fig13b_K7_psi=4_measured.csv",

    # Fig 14: K9 at delta=90, psi=1 -- Ito-Imai measured data + 2 paper calc lines
    "docs/junction/bassett_fig14_ito_imai_90_data_k9_q_squares.csv":   "bassett2001/bassett_fig14_K9_theta=90_psi=1_measured_ito_imai.csv",
    "docs/junction/bassett_fig14_90deg_calc_solid.csv":                "bassett2001/bassett_fig14_K9_theta=90_psi=1_calc_raw.csv",
    "docs/junction/bassett_fig14_90deg_calc_angle_corr_dashed.csv":    "bassett2001/bassett_fig14_K9_theta=90_psi=1_calc_corr.csv",
}


# ---------------------------------------------------------------------------
# Hager 1984 (currently in docs/bassett/ since user grouped them there)
# ---------------------------------------------------------------------------

HAGER_MAPPINGS: dict[str, str] = {
    # Fig 4a: xi_t measured at delta=90, two datasets (single, manifold) from Hager [3]
    "docs/bassett/hager_fig4a_zt_q_single.csv":   "hager1984/hager_fig04a_xi_t_delta=90_single_measured.csv",
    "docs/bassett/hager_fig4a_zt_q_manifold.csv": "hager1984/hager_fig04a_xi_t_delta=90_manifold_measured.csv",
    # Fig 4b: xi_t arbitrary branch geometry envelope (mean + upper + lower)
    "docs/bassett/hager_fig4b_zt_q_range_mean.csv":   "hager1984/hager_fig04b_xi_t_envelope=mean.csv",
    "docs/bassett/hager_fig4b_zt_q_range_ubound.csv": "hager1984/hager_fig04b_xi_t_envelope=upper.csv",
    "docs/bassett/hager_fig4b_zt_q_range_lbound.csv": "hager1984/hager_fig04b_xi_t_envelope=lower.csv",
}


# ---------------------------------------------------------------------------
# Wang 2014 -- clean up the bulk-rename k12_q_ artifact and add area-ratio
# ---------------------------------------------------------------------------

def _wang_mappings() -> dict[str, str]:
    """Build Wang mappings by listing the source directory.

    The user's bulk-rename left a `k12_q_` artifact and one stray uppercase
    `M3` variant (`wang_fig11_k12_q_k13_M3_q=0.2.csv`). macOS APFS is
    case-insensitive so the same physical file can be matched twice via
    different string casings -- dedup by canonical lowercase basename.
    """
    import re

    area_for_fig = {10: "1", 11: "1.56", 12: "2.44"}
    pattern = re.compile(
        r"^wang_fig(?P<fig>\d+)_k12_q_(?P<K_id>k1[3]|k2[3])_(?:m3|M3)_q=(?P<q>[\d.]+)\.csv$"
    )
    mappings: dict[str, str] = {}
    seen_lower: set[str] = set()
    src_dir = REPO_ROOT / "docs" / "junction"
    for p in sorted(src_dir.glob("wang_*.csv")):
        if p.name.lower() in seen_lower:
            continue
        m = pattern.match(p.name)
        if not m:
            continue
        fig = int(m.group("fig"))
        K_id = m.group("K_id").upper()  # k13 -> K13
        q = m.group("q")
        a = area_for_fig[fig]
        dst = f"wang2014/wang_fig{fig}_{K_id}_a={a}_q={q}_measured.csv"
        mappings[f"docs/junction/{p.name}"] = dst
        seen_lower.add(p.name.lower())
    return mappings


# ---------------------------------------------------------------------------
# Migration driver
# ---------------------------------------------------------------------------


def _git_tracked(rel: Path) -> bool:
    res = subprocess.run(
        ["git", "ls-files", "--error-unmatch", str(rel)],
        cwd=REPO_ROOT, capture_output=True, text=True
    )
    return res.returncode == 0


def migrate(mappings: dict[str, str], *, dry_run: bool) -> tuple[int, int, list[str]]:
    moved = 0
    skipped = 0
    missing: list[str] = []
    for src_rel, dst_rel in mappings.items():
        src = REPO_ROOT / src_rel
        dst = DST_ROOT / dst_rel
        if dst.exists():
            skipped += 1
            continue
        if not src.exists():
            missing.append(src_rel)
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dry_run:
            print(f"  DRY: {src_rel} -> {dst.relative_to(REPO_ROOT)}")
        else:
            tracked = _git_tracked(src.relative_to(REPO_ROOT))
            if tracked:
                # Use git mv to preserve history
                subprocess.run(
                    ["git", "mv", str(src.relative_to(REPO_ROOT)), str(dst.relative_to(REPO_ROOT))],
                    cwd=REPO_ROOT, check=True
                )
            else:
                shutil.move(str(src), str(dst))
            print(f"  mv: {src_rel} -> {dst.relative_to(REPO_ROOT)}")
        moved += 1
    return moved, skipped, missing


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="Preview without moving files")
    args = parser.parse_args()

    all_mappings: dict[str, str] = {}
    all_mappings.update(BASSETT_MAPPINGS)
    all_mappings.update(HAGER_MAPPINGS)
    all_mappings.update(_wang_mappings())

    print(f"Planning to migrate {len(all_mappings)} files...")
    if args.dry_run:
        print("DRY RUN -- no files will be moved.\n")
    print()
    moved, skipped, missing = migrate(all_mappings, dry_run=args.dry_run)
    print()
    print(f"Result: {moved} migrated, {skipped} skipped (already at dst), {len(missing)} src missing.")
    if missing:
        print("\nMissing sources (typo in mapping?):")
        for m in missing:
            print(f"  - {m}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
