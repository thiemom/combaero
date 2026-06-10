"""
Paper-published analytical-correlation ceiling reference.

Dispatches each (paper, K_id) to the paper's own correlation module so the
runner can compute Delta_vs_ceiling for any measured datum. This is what
defines "best 1D physics can do" per the literature.
"""

from __future__ import annotations

import math

from validation.junction.models import (
    bassett2001,
    hager1984,
    perez_garcia2010,
    wang2014,
)


class PaperCeiling:
    """CorrelationModel that delegates per paper-K_id to its source module.

    For Bassett joining K's, uses the body-text angle-corrected forms
    (Eqs 33/34/35/37) which Bassett himself shows in Figs 10b/11/12b/14
    as the "preferred" curves.
    """

    name = "paper_analytical_ceiling"

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
        if paper == "bassett2001":
            return self._bassett(K_id, q, psi, theta_rad)
        if paper == "hager1984":
            return self._hager(K_id, q, theta_rad)
        if paper == "wang2014":
            return self._wang(K_id, q, M_3)
        if paper == "perez_garcia2010":
            return self._perez_garcia(K_id, q, M_3)
        return None

    def _bassett(
        self, K_id: str, q: float, psi: float | None, theta_rad: float | None
    ) -> float | None:
        if K_id == "xi":
            # Fig 5 contraction coefficient: not a loss K.
            return None
        if psi is None:
            psi = 1.0
        if theta_rad is None:
            theta_rad = math.pi / 2.0
        try:
            return bassett2001.evaluate(K_id, q, psi, theta_rad, angle_corrected=True)
        except ValueError:
            return None

    def _hager(self, K_id: str, q: float, theta_rad: float | None) -> float | None:
        if K_id == "xi_t":
            return hager1984.xi_t(q)
        if K_id == "xi_l":
            if theta_rad is None:
                return None
            return hager1984.xi_l(q, theta_rad)
        return None

    def _wang(self, K_id: str, q: float, M_3: float | None) -> float | None:
        if M_3 is None:
            return None
        try:
            if K_id == "K_13":
                return (
                    wang2014.lookup_table1(round(M_3, 6), "K_13") if abs(M_3 - 0.5) > 0.5 else None
                )
            if K_id == "K_23":
                return None  # Wang's tables are at fixed q=0.5 or fixed M=0.5 only
        except ValueError:
            return None
        return None

    def _perez_garcia(self, K_id: str, q: float, M_3: float | None) -> float | None:
        if M_3 is None:
            return None
        # Map paper K_id to Perez-Garcia's K_hat_1/K_hat_2 family.
        # Conservative: only return for explicit K_hat_1/K_hat_2 ids.
        if K_id not in {"K_hat_1", "K_hat_2"}:
            return None
        try:
            return perez_garcia2010.K_hat("C2", K_id, q, M_3)
        except (ValueError, NotImplementedError):
            return None
