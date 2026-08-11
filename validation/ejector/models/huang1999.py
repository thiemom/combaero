"""1-D critical-mode supersonic ejector model: research/validation narrative.

The canonical implementation lives in `combaero.network._ejector_huang1999`
-- `EjectorElement` depends on it at runtime, so it must ship inside the
`combaero` wheel rather than this dev-only `validation/` tree (same
packaging lesson as `validation/junction/models/mynard2010.py`: PyPI install
of combaero-gui 0.4.0 hit `ModuleNotFoundError: No module named
'validation'` before that move; see CHANGELOG). This module re-exports it so
existing validation scripts/tests keep their import path unchanged, and
carries the full paper-by-paper research narrative that belongs here, not in
production code.

References
----------
[Huang 1999] B.J. Huang, J.M. Chang, C.P. Wang, V.A. Petrenko,
    "A 1-D analysis of ejector performance",
    International Journal of Refrigeration 22 (1999) 354-364.
    Local copy: docs/ejector/Huang_1d_analysis_ejector.pdf
    Source of `entrainment_ratio` (Eqs. 1-8) and of the original P_c* chain
    (Eqs. 9-18: phi_m-weighted momentum + Rankine-Hugoniot shock + isentropic
    diffuser) that `critical_back_pressure` no longer uses -- see below.

[Kracik & Dvorak 2016] J. Kracik, V. Dvorak,
    "Development of an Analytical Method for Predicting Flow in a Supersonic
    Air Ejector", EPJ Web of Conferences 114, 02059 (2016),
    DOI: 10.1051/epjconf/201611402059.
    Local copy: docs/ejector/Development_of_an_Analytical_Method_for.pdf
    Source of the mixing closure `critical_back_pressure` implements (their
    Eqs. 7-13, a lambda/q(lambda)/z(lambda) formalism after Christianovic,
    Kiselev and Abramovich): mass, momentum and energy are solved together
    directly from the y-y state to the fully-mixed, subsonic state, with NO
    separate loss coefficient on momentum and NO explicit shock -- conjugate
    supersonic/subsonic roots of the same conservation laws play the shock's
    role, the same way a Fanno/Rayleigh line has two roots.

[Akbarnejad & Ziabasharhagh 2025] S. Akbarnejad, M. Ziabasharhagh,
    "A novel 1D model for the analysis of double-choked ejectors validated by
    CFD simulations", J. Braz. Soc. Mech. Sci. Eng. 47:253 (2025),
    DOI: 10.1007/s40430-025-05536-7.
    Local copy: docs/ejector/art_10.1007_s40430-025-05536-7.pdf
    Independent CFD confirmation that a naive momentum balance across the
    mixing region (as used by Huang and predecessors) is unreliable: their
    Psi_5 metric (relative error of assuming momentum is simply conserved)
    reaches up to 30% at critical back pressure -- close to the ~25-35% gap
    that motivated moving off Huang's original P_c* chain (below). Their own
    proposed fix (drop momentum entirely; close on an externally-imposed
    exact-sonic condition) does not transplant onto a given-geometry forward
    analysis like this one -- it is a geometry-free DESIGN tool, and the
    paper says so directly -- so it was evaluated and not adopted; see
    validation/ejector/README.md's Accuracy section for the comparison.

[Lienhard/McGovern] R.K. McGovern, K.V. Bulusu, M.A. Antar, J.H. Lienhard V,
    "One-dimensional Model of an Optimal Ejector and Parametric Study of
    Ejector Efficiency" (MIT/GWU/KFUPM).
    Local copy: docs/ejector/Lienhard_One dimensional.pdf
    A third, independent source using the same momentum-retaining,
    shock-free mixing pattern as Kracik & Dvorak (their Eqs. 8-11): further
    corroboration, not the equations implemented here.

Why the P_c* chain changed from Huang's original Eqs. 9-18: tested against
all 31 digitized Table 3/4 rows, Huang's own phi_m + shock chain is
systematically low (mean 25%, max 35%) against each row's T_c*-derived
target -- confirmed NOT a gamma effect (P_c* is monotonic in gamma over its
whole valid range, gamma -> 1 limit still short; re-anchoring gamma from
primary stagnation to the local mixed state moves P_c* under 1%) and not an
equation transcription error (hand-verified, and Eqs. 16-17 reproduce
classical normal-shock-table values exactly). Kracik & Dvorak's closure,
applied downstream of the same y-y state, cuts this to mean 6.2%/max 13.2%
-- see README.md's Accuracy section for the full comparison including the
Lienhard/McGovern and Akbarnejad alternatives that were also tested.

Reproduction of the paper's Tables 3-4 (theory column) is exercised in
python/tests/test_ejector_huang.py. Using a defined ideal-gas gamma for R141b
evaluated at the entrained-flow choking plane from the CoolProp EOS (~1.111; see
huang1999_tables.py) the model matches the published theoretical omega to a mean
of 0.8% (max ~2%) across both primary pressures, both suction conditions, and
both nozzles. P_c* (recovery_efficiency=1.0) matches each row's T_c*-derived
target to a mean of 6.2% (max ~13.2%) -- see README.md's Accuracy section.
"""

from __future__ import annotations

from combaero.network._ejector_huang1999 import (
    ETA_P,
    ETA_S,
    PHI_P,
    CriticalPressureResult,
    EjectorGeometry,
    EjectorResult,
    choked_mass_flow,
    critical_back_pressure,
    entrainment_ratio,
)

__all__ = [
    "ETA_P",
    "ETA_S",
    "PHI_P",
    "CriticalPressureResult",
    "EjectorGeometry",
    "EjectorResult",
    "choked_mass_flow",
    "critical_back_pressure",
    "entrainment_ratio",
]
