"""Digitized Figs. 5a and 7a of Chen et al. (2013), Int. J. Refrigeration 36.

W. Chen, M. Liu, D. Chong, J. Yan, A.B. Little, Y. Bartosiewicz, "A 1D model
to predict ejector performance at critical and sub-critical operational
regimes", International Journal of Refrigeration 36 (2013) 1750-1761,
DOI 10.1016/j.ijrefrig.2013.04.009. Source PDF:
sources/ejector/Chen2013_1D_model_critical_subcritical.pdf (gitignored).

Chen's model is Huang 1999's chain (Eqs. 9-18 with phi_m, Eq. 19) run forward
at the experimental mixing area, plus a subcritical branch (the entrained flow
unchokes, P_sy rises until the chain's back pressure equals P_c). Coefficients
are Huang's: eta_p 0.95, eta_py 0.88, eta_s 0.85; phi_m 0.84 for both rigs.

What each figure holds (omega against the back pressure P_c [bar]):

* **Fig. 5a**, the Hemidi et al. (2009) air ejector: d_t 3.3, d_1 4.5,
  d_2 7.6 mm. Chen runs it at P_p = 4/5/6 bar, suction ~1 bar, 20 C. Its
  measured points are Henry 2007 Fig. 5's data (same rig, same lab):
  34 of Henry's 48 points, each within 0.010 in omega and 0.012 in P_c of
  henry2007_fig5.py (a uniform +0.005 omega offset between the two figures).
  Its model points are Chen's OWN evaluation of his model: a fidelity target.
* **Fig. 7a**, Chen's large rectangular air ejector (MEASURED, error bars):
  primary throat height 6.07 mm, nozzle exit 8.35 mm, mixing section 27.06 mm,
  all 50 mm wide (A_t 303.5 mm2, A_2 1353 mm2), suction from ambient,
  P_p = 3.0-5.0 bar. The efficiency-free test (ratios of plateau omega across
  P_p) reads P_p as ABSOLUTE (rms 0.010 against 0.038 for gauge).

Digitization (2026-10-19): 400 dpi render, axes calibrated on every tick
(13+21 ticks in 5a, 11+23 in 7a; held-out tick residual < 0.0006). Points
identified by eye on gridded tiles, then snapped to the scan: model markers by
template, measured points by their error-bar cap pair (centre = mid-caps;
``half_bar`` = half the bar height, the paper's measurement band). Every point
reviewed by the user on a numbered scan overlay. Twelve were flagged as
overlapped (a measured disk under a model circle, or a model marker on an
error bar) and kept after review. The points sit on a P_c grid (0.05 bar
steps plus 1.03); no point exists at 7a 4.0-measured 1.55, 4.5-measured 1.70
or 4.5-model 1.70.

Rows: (id, kind, P_p [bar], P_c [bar], omega, half_bar or None). ``kind`` is
"model" (Chen's model curve) or "measured". The ids match the review page.
These are digitized coordinates (facts); the figures stay in the gitignored PDF.
"""

from __future__ import annotations

#: (A_1/A_t, A_2/A_t) per figure: nozzle exit and mixing section over the throat.
AREA_RATIOS: dict[str, tuple[float, float]] = {
    "5a": ((4.5 / 3.3) ** 2, (7.6 / 3.3) ** 2),
    "7a": (8.35 / 6.07, 27.06 / 6.07),
}

#: Fig. 5a. Sorted model first, then measured; by P_p, then P_c.
FIG5A: list[tuple[str, str, float, float, float, float | None]] = [
    ("5a-M4-1", "model", 4.0, 1.0301, 0.8950, None),
    ("5a-M4-2", "model", 4.0, 1.1273, 0.8950, None),
    ("5a-M4-3", "model", 4.0, 1.2286, 0.8950, None),
    ("5a-M4-4", "model", 4.0, 1.2761, 0.8950, None),
    ("5a-M4-5", "model", 4.0, 1.3226, 0.8584, None),
    ("5a-M4-6", "model", 4.0, 1.3479, 0.7841, None),
    ("5a-M4-7", "model", 4.0, 1.3701, 0.7310, None),
    ("5a-M4-8", "model", 4.0, 1.4165, 0.5929, None),
    ("5a-M4-9", "model", 4.0, 1.4619, 0.4843, None),
    ("5a-M4-10", "model", 4.0, 1.5105, 0.3532, None),
    ("5a-M5-1", "model", 5.0, 1.0296, 0.6808, None),
    ("5a-M5-2", "model", 5.0, 1.1320, 0.6808, None),
    ("5a-M5-3", "model", 5.0, 1.2291, 0.6808, None),
    ("5a-M5-4", "model", 5.0, 1.3284, 0.6808, None),
    ("5a-M5-5", "model", 5.0, 1.4308, 0.6808, None),
    ("5a-M5-6", "model", 5.0, 1.4804, 0.6808, None),
    ("5a-M5-7", "model", 5.0, 1.5247, 0.6253, None),
    ("5a-M5-8", "model", 5.0, 1.5743, 0.5427, None),
    ("5a-M5-9", "model", 5.0, 1.5986, 0.4790, None),
    ("5a-M5-10", "model", 5.0, 1.6377, 0.3987, None),
    ("5a-M5-11", "model", 5.0, 1.6630, 0.3597, None),
    ("5a-M5-12", "model", 5.0, 1.6852, 0.3066, None),
    ("5a-M5-13", "model", 5.0, 1.7126, 0.2724, None),
    ("5a-M5-14", "model", 5.0, 1.7591, 0.1945, None),
    ("5a-M5-15", "model", 5.0, 1.7981, 0.1378, None),
    ("5a-M6-1", "model", 6.0, 1.0412, 0.5427, None),
    ("5a-M6-2", "model", 6.0, 1.1352, 0.5427, None),
    ("5a-M6-3", "model", 6.0, 1.2344, 0.5427, None),
    ("5a-M6-4", "model", 6.0, 1.3252, 0.5427, None),
    ("5a-M6-5", "model", 6.0, 1.4888, 0.5427, None),
    ("5a-M6-6", "model", 6.0, 1.5332, 0.5427, None),
    ("5a-M6-7", "model", 6.0, 1.6430, 0.5427, None),
    ("5a-M6-8", "model", 6.0, 1.7147, 0.5155, None),
    ("5a-M6-9", "model", 6.0, 1.8066, 0.3822, None),
    ("5a-M6-10", "model", 6.0, 1.9005, 0.2535, None),
    ("5a-M6-11", "model", 6.0, 2.0029, 0.1343, None),
    ("5a-E4-1", "measured", 4.0, 1.0301, 0.7581, None),
    ("5a-E4-2", "measured", 4.0, 1.1262, 0.7581, None),
    ("5a-E4-3", "measured", 4.0, 1.2286, 0.7581, None),
    ("5a-E4-4", "measured", 4.0, 1.2761, 0.7581, None),
    ("5a-E4-5", "measured", 4.0, 1.3226, 0.7345, None),
    ("5a-E4-6", "measured", 4.0, 1.3479, 0.7062, None),
    ("5a-E4-7", "measured", 4.0, 1.3701, 0.6601, None),
    ("5a-E4-8", "measured", 4.0, 1.4165, 0.5386, None),
    ("5a-E4-9", "measured", 4.0, 1.4609, 0.4311, None),
    ("5a-E4-10", "measured", 4.0, 1.5105, 0.3143, None),
    ("5a-E5-1", "measured", 5.0, 1.0296, 0.6265, None),
    ("5a-E5-2", "measured", 5.0, 1.1320, 0.6265, None),
    ("5a-E5-3", "measured", 5.0, 1.2291, 0.6265, None),
    ("5a-E5-4", "measured", 5.0, 1.3284, 0.6265, None),
    ("5a-E5-5", "measured", 5.0, 1.4308, 0.6265, None),
    ("5a-E5-6", "measured", 5.0, 1.4804, 0.6265, None),
    ("5a-E5-7", "measured", 5.0, 1.5247, 0.6135, None),
    ("5a-E5-8", "measured", 5.0, 1.5743, 0.5557, None),
    ("5a-E5-9", "measured", 5.0, 1.5986, 0.5120, None),
    ("5a-E5-10", "measured", 5.0, 1.6377, 0.4081, None),
    ("5a-E5-11", "measured", 5.0, 1.6630, 0.3586, None),
    ("5a-E5-12", "measured", 5.0, 1.6852, 0.3196, None),
    ("5a-E5-13", "measured", 5.0, 1.7126, 0.2700, None),
    ("5a-E5-14", "measured", 5.0, 1.7591, 0.2110, None),
    ("5a-E6-1", "measured", 6.0, 1.0412, 0.5096, None),
    ("5a-E6-2", "measured", 6.0, 1.1352, 0.5096, None),
    ("5a-E6-3", "measured", 6.0, 1.2344, 0.5096, None),
    ("5a-E6-4", "measured", 6.0, 1.3252, 0.5132, None),
    ("5a-E6-5", "measured", 6.0, 1.4888, 0.5096, None),
    ("5a-E6-6", "measured", 6.0, 1.5332, 0.5073, None),
    ("5a-E6-7", "measured", 6.0, 1.6430, 0.4730, None),
    ("5a-E6-8", "measured", 6.0, 1.7147, 0.4424, None),
    ("5a-E6-9", "measured", 6.0, 1.8995, 0.3149, None),
    ("5a-E6-10", "measured", 6.0, 2.0029, 0.1591, None),
]

#: Fig. 7a. Sorted model first, then measured; by P_p, then P_c.
FIG7A: list[tuple[str, str, float, float, float, float | None]] = [
    ("7a-M3.0-1", "model", 3.0, 1.0305, 0.9939, None),
    ("7a-M3.0-2", "model", 3.0, 1.0507, 0.9939, None),
    ("7a-M3.0-3", "model", 3.0, 1.0998, 0.9939, None),
    ("7a-M3.0-4", "model", 3.0, 1.1498, 0.9939, None),
    ("7a-M3.0-5", "model", 3.0, 1.1998, 0.9939, None),
    ("7a-M3.0-6", "model", 3.0, 1.2499, 0.8980, None),
    ("7a-M3.0-7", "model", 3.0, 1.2999, 0.7646, None),
    ("7a-M3.0-8", "model", 3.0, 1.3499, 0.6531, None),
    ("7a-M3.0-9", "model", 3.0, 1.3999, 0.5417, None),
    ("7a-M3.5-1", "model", 3.5, 1.1981, 0.8255, None),
    ("7a-M3.5-2", "model", 3.5, 1.2481, 0.8255, None),
    ("7a-M3.5-3", "model", 3.5, 1.2981, 0.8255, None),
    ("7a-M3.5-4", "model", 3.5, 1.3481, 0.7607, None),
    ("7a-M3.5-5", "model", 3.5, 1.3981, 0.6246, None),
    ("7a-M3.5-6", "model", 3.5, 1.4482, 0.4886, None),
    ("7a-M3.5-7", "model", 3.5, 1.4982, 0.3577, None),
    ("7a-M3.5-8", "model", 3.5, 1.5482, 0.2385, None),
    ("7a-M3.5-9", "model", 3.5, 1.5982, 0.1362, None),
    ("7a-M4.0-1", "model", 4.0, 1.2999, 0.7011, None),
    ("7a-M4.0-2", "model", 4.0, 1.3499, 0.7011, None),
    ("7a-M4.0-3", "model", 4.0, 1.3981, 0.7050, None),
    ("7a-M4.0-4", "model", 4.0, 1.4499, 0.6765, None),
    ("7a-M4.0-5", "model", 4.0, 1.4999, 0.5624, None),
    ("7a-M4.0-6", "model", 4.0, 1.5499, 0.4432, None),
    ("7a-M4.0-7", "model", 4.0, 1.5999, 0.3344, None),
    ("7a-M4.0-8", "model", 4.0, 1.6500, 0.2372, None),
    ("7a-M4.0-9", "model", 4.0, 1.6991, 0.1530, None),
    ("7a-M4.5-1", "model", 4.5, 1.3981, 0.5974, None),
    ("7a-M4.5-2", "model", 4.5, 1.4482, 0.5974, None),
    ("7a-M4.5-3", "model", 4.5, 1.4982, 0.5974, None),
    ("7a-M4.5-4", "model", 4.5, 1.5482, 0.5974, None),
    ("7a-M4.5-5", "model", 4.5, 1.6008, 0.5223, None),
    ("7a-M4.5-6", "model", 4.5, 1.6482, 0.4160, None),
    ("7a-M4.5-7", "model", 4.5, 1.7482, 0.2372, None),
    ("7a-M4.5-8", "model", 4.5, 1.8009, 0.1673, None),
    ("7a-M5.0-1", "model", 5.0, 1.4999, 0.5236, None),
    ("7a-M5.0-2", "model", 5.0, 1.5499, 0.5236, None),
    ("7a-M5.0-3", "model", 5.0, 1.5999, 0.5236, None),
    ("7a-M5.0-4", "model", 5.0, 1.6500, 0.5236, None),
    ("7a-M5.0-5", "model", 5.0, 1.7000, 0.4860, None),
    ("7a-M5.0-6", "model", 5.0, 1.7500, 0.4031, None),
    ("7a-M5.0-7", "model", 5.0, 1.8000, 0.3215, None),
    ("7a-M5.0-8", "model", 5.0, 1.8500, 0.2463, None),
    ("7a-M5.0-9", "model", 5.0, 1.9000, 0.1802, None),
    ("7a-E3.0-1", "measured", 3.0, 1.0296, 0.8844, 0.0097),
    ("7a-E3.0-2", "measured", 3.0, 1.0498, 0.8844, 0.0097),
    ("7a-E3.0-3", "measured", 3.0, 1.0998, 0.8838, 0.0104),
    ("7a-E3.0-4", "measured", 3.0, 1.1498, 0.8838, 0.0104),
    ("7a-E3.0-5", "measured", 3.0, 1.1998, 0.8838, 0.0104),
    ("7a-E3.0-6", "measured", 3.0, 1.2499, 0.8624, 0.0097),
    ("7a-E3.0-7", "measured", 3.0, 1.2999, 0.7730, 0.0136),
    ("7a-E3.0-8", "measured", 3.0, 1.3499, 0.6382, 0.0149),
    ("7a-E3.0-9", "measured", 3.0, 1.3999, 0.4471, 0.0052),
    ("7a-E3.5-1", "measured", 3.5, 1.1998, 0.7510, 0.0084),
    ("7a-E3.5-2", "measured", 3.5, 1.2499, 0.7510, 0.0084),
    ("7a-E3.5-3", "measured", 3.5, 1.2999, 0.7510, 0.0084),
    ("7a-E3.5-4", "measured", 3.5, 1.3499, 0.7458, 0.0084),
    ("7a-E3.5-5", "measured", 3.5, 1.3990, 0.6888, 0.0084),
    ("7a-E3.5-6", "measured", 3.5, 1.4499, 0.5812, 0.0071),
    ("7a-E3.5-7", "measured", 3.5, 1.4999, 0.4407, 0.0052),
    ("7a-E3.5-8", "measured", 3.5, 1.5499, 0.2703, 0.0045),
    ("7a-E3.5-9", "measured", 3.5, 1.5999, 0.1323, 0.0052),
    ("7a-E4.0-1", "measured", 4.0, 1.2999, 0.6305, 0.0071),
    ("7a-E4.0-2", "measured", 4.0, 1.3499, 0.6311, 0.0078),
    ("7a-E4.0-3", "measured", 4.0, 1.3999, 0.6305, 0.0071),
    ("7a-E4.0-4", "measured", 4.0, 1.4499, 0.6305, 0.0071),
    ("7a-E4.0-5", "measured", 4.0, 1.4999, 0.6130, 0.0078),
    ("7a-E4.0-6", "measured", 4.0, 1.5999, 0.3953, 0.0052),
    ("7a-E4.0-7", "measured", 4.0, 1.6500, 0.2755, 0.0032),
    ("7a-E4.0-8", "measured", 4.0, 1.7000, 0.1426, 0.0039),
    ("7a-E4.5-1", "measured", 4.5, 1.3999, 0.5404, 0.0065),
    ("7a-E4.5-2", "measured", 4.5, 1.4499, 0.5404, 0.0065),
    ("7a-E4.5-3", "measured", 4.5, 1.4999, 0.5411, 0.0058),
    ("7a-E4.5-4", "measured", 4.5, 1.5499, 0.5385, 0.0084),
    ("7a-E4.5-5", "measured", 4.5, 1.5999, 0.5365, 0.0065),
    ("7a-E4.5-6", "measured", 4.5, 1.6500, 0.4834, 0.0117),
    ("7a-E4.5-7", "measured", 4.5, 1.7500, 0.2709, 0.0039),
    ("7a-E4.5-8", "measured", 4.5, 1.7991, 0.1647, 0.0065),
    ("7a-E5.0-1", "measured", 5.0, 1.4999, 0.4698, 0.0058),
    ("7a-E5.0-2", "measured", 5.0, 1.5499, 0.4698, 0.0058),
    ("7a-E5.0-3", "measured", 5.0, 1.5999, 0.4698, 0.0058),
    ("7a-E5.0-4", "measured", 5.0, 1.6500, 0.4743, 0.0104),
    ("7a-E5.0-5", "measured", 5.0, 1.7000, 0.4698, 0.0058),
    ("7a-E5.0-6", "measured", 5.0, 1.7500, 0.4445, 0.0052),
    ("7a-E5.0-7", "measured", 5.0, 1.8000, 0.3739, 0.0045),
    ("7a-E5.0-8", "measured", 5.0, 1.8500, 0.2748, 0.0039),
    ("7a-E5.0-9", "measured", 5.0, 1.9000, 0.1848, 0.0032),
]
