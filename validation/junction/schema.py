"""
Data-loading schema for the junction validation dataset.

Each paper folder has a metadata.yaml describing every CSV. This module
parses those yaml files into typed dataclasses and exposes a single
load_dataset() entrypoint.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import yaml

DATA_ROOT = Path(__file__).parent / "data"

Kind = Literal["measured", "calc", "envelope", "tabulated"]
FlowRegime = Literal["separating", "joining", "combining"]
FitTier = Literal["excellent", "very_good", "good", "moderate", "poor"]
QConvention = Literal["bassett", "hager", "wang", "perez_garcia", "idelchik"]


@dataclass(frozen=True)
class PaperMetadata:
    name: str
    citation: str
    doi: str | None
    q_convention: QConvention
    default_uncertainty_K: float | None
    fit_tiers: dict[str, FitTier]


@dataclass(frozen=True)
class FileMetadata:
    """Metadata for one CSV. K_id is the paper's own notation
    (K1..K12 for Bassett, xi_t/xi_l for Hager, K_13/K_23 for Wang,
    K_hat_1/K_hat_2 for Perez-Garcia)."""

    path: Path
    paper: PaperMetadata
    K_id: str | None  # None if file is e.g. xi contraction coefficient (no K)
    coefficient: str | None  # e.g. 'xi_t' for Hager; defaults to K_id otherwise
    psi: float | None
    theta_deg: float | None
    a: float | None  # Wang area ratio
    q: float | None  # fixed q when x_axis != q
    M_3: float | None  # fixed M_3 when relevant
    x_axis: Literal["q", "M_3"]
    flow_regime: FlowRegime | None
    kind: Kind
    calc_variant: str | None  # 'raw' / 'angle_corrected' for paper calc lines
    envelope_role: str | None  # 'mean' / 'upper' / 'lower' for Hager Fig 4b
    dataset: str | None  # 'single_laterals' / 'manifold_laterals' for Hager Fig 4a
    source: str | None  # external source within paper (e.g. 'Ito and Imai 1973')

    @property
    def fit_tier(self) -> FitTier | None:
        key = self.K_id or self.coefficient
        if key is None:
            return None
        return self.paper.fit_tiers.get(key)

    @property
    def uncertainty_K(self) -> float | None:
        return self.paper.default_uncertainty_K


def _load_one_paper(paper_dir: Path) -> tuple[PaperMetadata, list[FileMetadata]]:
    yaml_path = paper_dir / "metadata.yaml"
    with open(yaml_path) as f:
        raw = yaml.safe_load(f)
    paper = PaperMetadata(
        name=raw["paper"]["name"],
        citation=raw["paper"]["citation"],
        doi=raw["paper"].get("doi"),
        q_convention=raw["paper"]["q_convention"],
        default_uncertainty_K=raw["paper"].get("default_uncertainty_K"),
        fit_tiers=raw.get("fit_tiers", {}),
    )
    files: list[FileMetadata] = []
    for fname, entry in (raw.get("files") or {}).items():
        e: dict[str, Any] = entry or {}
        files.append(
            FileMetadata(
                path=paper_dir / fname,
                paper=paper,
                K_id=e.get("K_id"),
                coefficient=e.get("coefficient") or e.get("K_id"),
                psi=e.get("psi"),
                theta_deg=e.get("theta_deg"),
                a=e.get("a"),
                q=e.get("q"),
                M_3=e.get("M_3"),
                x_axis=e.get("x_axis", "q"),
                flow_regime=e.get("flow_regime"),
                kind=e.get("kind", "measured"),
                calc_variant=e.get("calc_variant"),
                envelope_role=e.get("envelope_role"),
                dataset=e.get("dataset"),
                source=e.get("source"),
            )
        )
    return paper, files


@dataclass(frozen=True)
class Dataset:
    papers: dict[str, PaperMetadata]
    files: list[FileMetadata]

    def by_paper(self, paper_slug: str) -> list[FileMetadata]:
        return [f for f in self.files if f.path.parent.name == paper_slug]

    def measured(self) -> list[FileMetadata]:
        return [f for f in self.files if f.kind == "measured"]

    def by_K(self, K_id: str) -> list[FileMetadata]:
        return [f for f in self.files if f.K_id == K_id]


def load_dataset(root: Path = DATA_ROOT) -> Dataset:
    """Load every paper's metadata.yaml under data/."""
    papers: dict[str, PaperMetadata] = {}
    files: list[FileMetadata] = []
    for paper_dir in sorted(root.iterdir()):
        if not paper_dir.is_dir() or not (paper_dir / "metadata.yaml").exists():
            continue
        paper, paper_files = _load_one_paper(paper_dir)
        papers[paper_dir.name] = paper
        files.extend(paper_files)
    return Dataset(papers=papers, files=files)
