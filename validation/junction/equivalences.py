"""
Cross-paper K-coefficient equivalence table.

Each paper uses its own notation and q convention. This table maps
(paper_slug, K_id) -> (canonical_K, q_transform) so the scorecard can
roll up the same underlying physics across multiple sources.

q_transform takes the paper's q and returns the "canonical" q used by
the corresponding canonical K (Bassett convention: q = m_other / m_com,
where m_other is the flow leg defining the K).
"""

from __future__ import annotations

from typing import Callable

# Canonical K names match the Bassett K1..K12 set as the reference taxonomy.
EQUIV: dict[tuple[str, str], tuple[str, Callable[[float], float]]] = {
    # ----- Separating straight (Bassett K2 = K5, Hager xi_t) -----
    # Hager q is lateral/total -> Bassett K5 q is straight/total = 1 - q_Hager.
    ("bassett2001", "K2"): ("K_straight_sep", lambda q: q),
    ("bassett2001", "K5"): ("K_straight_sep", lambda q: q),
    ("hager1984", "xi_t"): ("K_straight_sep", lambda q: 1.0 - q),
    # ----- Separating lateral (Bassett K6, Hager xi_l at psi=1) -----
    ("bassett2001", "K6"): ("K_lateral_sep", lambda q: q),
    ("hager1984", "xi_l"): ("K_lateral_sep", lambda q: q),
    # ----- Joining lateral / branch-to-outlet (Bassett K12, Wang K_13) -----
    # Wang q = lateral / common = Bassett q for K12; identity transform.
    ("bassett2001", "K12"): ("K_lateral_join", lambda q: q),
    ("wang2014", "K_13"): ("K_lateral_join", lambda q: q),
    # ----- Joining straight (Bassett K11, Wang K_23) -----
    # Wang K_23 normalises straight-branch loss against common dynamic head,
    # same as Bassett K11. q is straight / common per Bassett but Wang uses
    # q = m_1 / m_3 = lateral / common; straight contribution = 1 - q.
    ("bassett2001", "K11"): ("K_straight_join", lambda q: q),
    ("wang2014", "K_23"): ("K_straight_join", lambda q: 1.0 - q),
    # K1, K3, K4, K7, K8, K9, K10 are separating-variant or type-4/5 joining
    # that have no cross-paper equivalent in the current dataset. They are
    # canonicalized as themselves (no rollup).
    ("bassett2001", "K1"): ("K1", lambda q: q),
    ("bassett2001", "K3"): ("K3", lambda q: q),
    ("bassett2001", "K4"): ("K4", lambda q: q),
    ("bassett2001", "K7"): ("K7", lambda q: q),
    ("bassett2001", "K8"): ("K8", lambda q: q),
    ("bassett2001", "K9"): ("K9", lambda q: q),
    ("bassett2001", "K10"): ("K10", lambda q: q),
}


def canonical_K(paper: str, K_id: str) -> str:
    """Return canonical K name for grouping in the scorecard. Falls back to
    "<paper>:<K_id>" if no entry; runner can warn on missing entries."""
    return EQUIV.get((paper, K_id), (f"{paper}:{K_id}", lambda q: q))[0]


def q_transform(paper: str, K_id: str) -> Callable[[float], float]:
    """Transform paper's q into the canonical-K q convention."""
    return EQUIV.get((paper, K_id), (None, lambda q: q))[1]
