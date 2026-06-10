"""
Junction validation runner -- correlation mode.

Iterates over the dataset's measured CSVs, evaluates a candidate model
plus the paper's published correlation (ceiling) at each data point,
emits per-row records. The scorecard aggregates.
"""

from __future__ import annotations

import csv
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Protocol

from validation.junction.equivalences import canonical_K
from validation.junction.schema import Dataset, FileMetadata, load_dataset


class CorrelationModel(Protocol):
    """Any candidate that can evaluate a published K id at given (q, psi, theta)."""

    name: str

    def evaluate(
        self,
        paper: str,
        K_id: str,
        q: float,
        psi: float | None,
        theta_rad: float | None,
        M_3: float | None = None,
        **kwargs: float,
    ) -> float | None:
        """Return K predicted by the model, or None if unsupported."""
        ...


@dataclass
class Record:
    """One (model evaluation, measured data point) pair."""

    paper: str
    file: str
    K_id: str
    canonical_K: str
    psi: float | None
    theta_deg: float | None
    q: float
    M_3: float | None
    K_measured: float
    K_model: float | None
    K_ceiling: float | None
    wall_time_s: float = 0.0


def _read_xy_csv(path: Path) -> list[tuple[float, float]]:
    """Read a digitized CSV with x, y columns. Skips an optional header row."""
    rows: list[tuple[float, float]] = []
    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 2:
                continue
            try:
                x = float(row[0].strip())
                y = float(row[1].strip())
            except ValueError:
                continue
            rows.append((x, y))
    return rows


def _is_supported(file: FileMetadata) -> bool:
    """Skip measured files we can't yet score (no K_id, e.g. Hager envelope)."""
    return file.kind == "measured" and (file.K_id is not None or file.coefficient is not None)


def iter_records(
    model: CorrelationModel,
    ceiling: CorrelationModel,
    dataset: Dataset,
) -> Iterator[Record]:
    """For each measured CSV in the dataset, yield one Record per data point."""
    for file in dataset.files:
        if not _is_supported(file):
            continue
        paper_slug = file.path.parent.name
        K_id = file.K_id or (file.coefficient or "")
        rows = _read_xy_csv(file.path)
        theta_rad = math.radians(file.theta_deg) if file.theta_deg is not None else None
        for x, y in rows:
            # Determine q and M_3 from x_axis convention.
            if file.x_axis == "q":
                q_val = x
                M_val = file.M_3
            else:  # M_3
                q_val = file.q if file.q is not None else 0.5
                M_val = x
            t0 = time.perf_counter()
            K_model = model.evaluate(
                paper_slug, K_id, q_val, file.psi, theta_rad, M_3=M_val, a=file.a
            )
            wall_time = time.perf_counter() - t0
            K_ceil = ceiling.evaluate(
                paper_slug, K_id, q_val, file.psi, theta_rad, M_3=M_val, a=file.a
            )
            yield Record(
                paper=paper_slug,
                file=file.path.name,
                K_id=K_id,
                canonical_K=canonical_K(paper_slug, K_id),
                psi=file.psi,
                theta_deg=file.theta_deg,
                q=q_val,
                M_3=M_val,
                K_measured=y,
                K_model=K_model,
                K_ceiling=K_ceil,
                wall_time_s=wall_time,
            )


def run(model: CorrelationModel, ceiling: CorrelationModel) -> list[Record]:
    """Convenience: load default dataset and return all records."""
    return list(iter_records(model, ceiling, load_dataset()))
