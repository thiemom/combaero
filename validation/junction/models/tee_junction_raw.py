"""
Adapter for the in-tree C++ tee_junction K1..K12 functions.

Exposes the *raw* Bassett K math compiled into the C++ engine -- no
network solver, no MCN propagation. Use to validate whether the code's
K transcriptions match the paper.

This adapter is expected to surface the K8 and K12 transcription bugs
flagged in include/tee_junction.h (missing psi factor in K12 bracket,
(1-q)^2/psi vs (1-q)*psi in K8).
"""

from __future__ import annotations

import math

import combaero as cb


class TeeJunctionRaw:
    """Wraps `cb._core` tee_junction K1..K12 calls."""

    name = "tee_junction_raw_cpp"

    # Only 4 of 12 Bassett K's are currently bound to Python via cb._core:
    #     tee_K5, tee_K6, tee_K11, tee_K12.
    # K1, K3, K4, K7, K8, K9, K10 are inline-only in include/tee_junction.h
    # and not exposed via pybind11. Adding those bindings is a separate
    # follow-up; the bound 4 cover the K12 ψ² bug verdict (the headline
    # case) so we proceed.
    _CORE_K_FNS = {
        "K5": "tee_K5",
        "K6": "tee_K6",
        "K11": "tee_K11",
        "K12": "tee_K12",
    }

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
        if paper != "bassett2001":
            # C++ raw K's are Bassett's. No cross-paper canonicalization here.
            return None
        if K_id not in self._CORE_K_FNS:
            return None
        if psi is None:
            psi = 1.0
        if theta_rad is None:
            theta_rad = math.pi / 2.0
        fn = getattr(cb._core, self._CORE_K_FNS[K_id], None)
        if fn is None:
            return None
        # K2/K5 take only q in the C++ API.
        if K_id in {"K2", "K5"}:
            return float(fn(q))
        return float(fn(q, psi, theta_rad))
