"""Digitized Figure 5 of Henry et al. (2007), HEFAT2007.

Henry, Leclaire, Hemidi, Seynhaeve, Bartosiewicz (2007), "Analysis of
Supersonic Ejector Operation: Experimental Validation and Two-Phase Aspects",
5th International Conference on Heat Transfer, Fluid Mechanics and
Thermodynamics (HEFAT2007), Sun City, South Africa, paper HF1. Source PDF:
docs/ejector/Bartosiewicz_ejector_2004.pdf (filename predates the confirmed
citation; the paper itself is HEFAT2007, Bartosiewicz corresponding author).

Figure 5, "Ejector characteristics - Experimental results", is the measured
ejector characteristic for an **air** ejector: entrainment ratio omega versus
normalized back pressure Pb / P2^0, for three primary driving pressures
P1^0 = 4, 5, 6 bar. P2^0 is the (roughly constant, near-ambient) secondary
stagnation pressure; the paper reports a whole-range CFD-vs-experiment
deviation below 10% for the k-epsilon model.

Why this fixture matters for the operating-regime extension
-----------------------------------------------------------
This is the off-design **validation** target for the subcritical-droop model
(see validation/ejector/OPERATING_REGIMES_DESIGN.md sec 3.2, 4). Each curve
shows the full regime structure the Tier-1 model must reproduce:

    plateau (critical, omega = omega_crit, flat)   Pb <= P_c*
    -> knee at the critical back pressure P_c*
    -> near-linear droop                           P_c* < Pb < P_b0
    -> breakdown (omega -> 0)                       Pb -> P_b0

and it is **air**, so it is directly comparable against combaero's EOS with no
refrigerant / CoolProp dependency (unlike the R141b Huang critical fixture).
It is used as a *validation* cross-check (trend + shape + anchor pressures),
NOT for parameter fitting -- see the calibration policy in CLAUDE.md and the
design doc.

Trends across the three curves (all physically expected, and reproduced here):
  - plateau omega DECREASES with primary pressure (0.760 -> 0.623 -> 0.504):
    stronger motive flow, lower entrainment ratio.
  - critical (knee) and breakdown back pressures INCREASE with primary
    pressure: a stronger motive stream holds a higher back pressure.

Digitization provenance and cross-check
---------------------------------------
Points were traced from the figure by hand (source: user, 2026-08-29) and
cross-checked against docs/ejector/bartosiewicz_fig5.png. Verdict: faithful.
Plateau values, knee locations, and breakdown endpoints all match the printed
figure to within trace error. The three breakdown endpoints
(1.587, 0.052) / (1.829, 0.073) / (2.047, 0.041) and the three plateaus
(0.760 / 0.623 / 0.504) are the checkable anchors and all agree. The 4 bar
point at Pb/P2 = 1.463 reads slightly low (omega ~ 0.43 vs ~ 0.45 in the
figure) but is well within hand-trace tolerance and does not affect the
plateau / knee / breakdown anchors used for validation.

These are digitized numeric coordinates (facts), stored here in the tracked
validation tree; the copyrighted source figure lives in the gitignored
docs/ejector/ PDF directory.
"""

from __future__ import annotations

# Measured entrainment characteristic omega(Pb / P2^0), air ejector.
# Keyed by primary driving pressure P1^0 [bar]. Each value is a list of
# (pb_over_p2, omega) points sorted by ascending normalized back pressure,
# from the on-design plateau through the subcritical droop to breakdown.
FIG5_OMEGA_VS_PB: dict[int, list[tuple[float, float]]] = {
    4: [
        (1.01784, 0.75929),
        (1.11808, 0.75618),
        (1.22182, 0.75784),
        (1.26374, 0.75364),
        (1.31762, 0.73344),
        (1.34460, 0.70491),
        (1.36490, 0.65866),
        (1.41794, 0.53507),
        (1.46277, 0.42732),
        (1.51233, 0.31245),
        (1.55957, 0.15760),
        (1.58734, 0.05161),
    ],
    5: [
        (1.02835, 0.62211),
        (1.12335, 0.62321),
        (1.22138, 0.62237),
        (1.31811, 0.62299),
        (1.42343, 0.62512),
        (1.47463, 0.62176),
        (1.50424, 0.62022),
        (1.51903, 0.60522),
        (1.54580, 0.57996),
        (1.56880, 0.54929),
        (1.58283, 0.52875),
        (1.59127, 0.51077),
        (1.61670, 0.47146),
        (1.63390, 0.40273),
        (1.64946, 0.37644),
        (1.66242, 0.35099),
        (1.68537, 0.31364),
        (1.71283, 0.26617),
        (1.74240, 0.24149),
        (1.75865, 0.20643),
        (1.79925, 0.13450),
        (1.82893, 0.07313),
    ],
    6: [
        (1.03566, 0.50472),
        (1.13117, 0.50318),
        (1.22718, 0.50365),
        (1.31968, 0.50342),
        (1.41903, 0.50180),
        (1.48292, 0.50107),
        (1.52466, 0.49845),
        (1.60288, 0.48187),
        (1.64185, 0.46718),
        (1.71108, 0.43828),
        (1.80151, 0.37706),
        (1.89686, 0.31074),
        (2.00005, 0.15148),
        (2.04703, 0.04074),
    ],
}

# Approximate anchor pressures read off each curve, in normalized units
# (Pb / P2^0). omega_plateau is the flat critical entrainment ratio;
# pb_knee is the critical back pressure P_c* (plateau end / start of droop);
# pb_breakdown is the last measured point (omega near zero ~ P_b0). These are
# eyeballed from the digitized curve for orientation only -- validation code
# should derive P_c*/P_b0 from the data via a consistent rule (e.g. plateau
# departure threshold and omega -> 0 extrapolation), not hard-code these.
FIG5_ANCHORS: dict[int, dict[str, float]] = {
    4: {"omega_plateau": 0.760, "pb_knee": 1.30, "pb_breakdown": 1.587},
    5: {"omega_plateau": 0.623, "pb_knee": 1.51, "pb_breakdown": 1.829},
    6: {"omega_plateau": 0.504, "pb_knee": 1.55, "pb_breakdown": 2.047},
}

PRIMARY_PRESSURES_BAR: tuple[int, ...] = (4, 5, 6)
