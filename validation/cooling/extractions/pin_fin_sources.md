# Pin-fin array sources (#335)

Every equation below was read from the page image, not the text layer
(scanned NASA reports misread digits, e.g. Faulkner's C2 3.094 OCRs as
5.094). Local copies are in `docs/heat_transfer/pin_fin/` (gitignored).
The scoping history, including dead ends, is on issue #335.

## Canonical basis used by `pin_fin_correlation.h`

| quantity | definition | source of the definition |
|---|---|---|
| Re_D | `mdot D / (mu A_min)`, velocity at the minimum flow area | Armstrong & Winstanley (1988) nomenclature |
| Nu_D | `h D / k`, h on the set's surface (total = pin + endwall, area-weighted) | same |
| f | `dP / (2 rho Vmax^2 N)` per row | same: "f = dP/2 rho V_max^2 N" |

**The code #332 removed** used `dP = N f rho Vmax^2 / 2`, a factor of 4 low
against this definition. It also had:
- a friction exponent of -0.25;
- a Pr^0.4 on Metzger;
- invented `(S/D/2.5)^-0.15 (X/D/2.5)^-0.10 (L/D)^-0.11` terms;
- an "inline 0.092", which matches Metzger 1982a's **staggered** X/D 1.5 coefficient (below; a probable mislabel, not confirmed).

## Geometry converters

Unit cell: one pin per `X S D^2`, in both arrangements.

| quantity | expression (per D) |
|---|---|
| open volume per cell | `H (X S - pi/4)` |
| wetted area per cell | `2 (X S - pi/4) + pi H` |
| `D'/D = 4V/A_t` | `H (4 X S - pi) / (2 X S + pi (H - 1/2))` (A&W Eq. 8, re-derived) |
| `A'/A_min` | `(S - pi/(4 X)) / gap` (Eq. 10 with `gap = S - 1`) |
| `D_h/D` | `4 X H gap / (2 X S + pi (H - 1/2))` (Eq. 13) |
| gap | `min(S - 1, 2 (sqrt((S/2)^2 + X^2) - 1))` for staggered; `S - 1` for inline |

- At VanFossen's H/D 2, S/D 4, X/D 3.464, `D'/D = 3.225`.
- Every source geometry here is transverse-limited.

Fin efficiency, A&W Eqs 14-16: `eta_fin = tanh(mL)/(mL)`, `m = sqrt(4h/(k D))`, `L = H/2`.

## Sets

| set | equation | where printed | geometry / box | notes |
|---|---|---|---|---|
| Metzger, Shepard & Haley 1986 (86-GT-132) | `Nu_D = 0.135 Re_D^0.69 (X/D)^-0.34` | A&W 1988 Eq. 2 (p. 96) | fitted H/D 1, S/D 2.5, 1.5 <= X/D <= 5, 10 rows, 1e3-1e5; A&W limits H/D <= 3, 2 <= S/D <= 4 | A&W: one point of 17 (Arora) outside +/-20%. Han Fig. 4.133 is the same data. |
| Metzger, Fan & Shepard 1982b (Heat Transfer 1982 v3, 137) | `f = 0.317 Re^-0.132` (1e3-1e4); `1.76 Re^-0.318` (1e4-1e5) | A&W Eqs 20-21 (p. 101) | fitted H/D 1, S/D 2.5, 1.5 <= X/D <= 5, +/-15% except X/D 1.79; A&W extend to 0.5 <= H/D <= 6, 2 <= S/D <= 4 on Peng | Branches meet at 1e4 (0.0940 vs 0.0941), slopes differ. Original paper inaccessible (Begell House): fidelity unchecked. |
| VanFossen 1982 (NASA TM-81696) | `Nu_D' = 0.153 Re_D'^0.685` | TM-81696 Eq. 16 (p. 5) | Table I: H/D 0.5 (S 2) and 2 (S 4), equilateral, 4 rows, Re_D' 300-6e4 | Nu from the Eq. 9 fin model with h_pin = h_wall. Wood pins give h_pin/h_wall = 1.345. |
| Damerow, Murtaugh & Burggraf 1972 (NASA CR-120883) | `f = 2.06 (X_T/D)^-1.1 Re^-0.16` | Eq. 18 (p. 42-44 figures); definition Eq. 10 (p. 19) | Fig. 3: X_T/D 4.24 and 7.07, X_L/D 2.12 and 3.54, Z/D 2-4, 10 rows | Eq. 10: `f = dP_T rho / (2 (N-1) G_min^2)`, total pressure first to last row. f rises above inlet M 0.36. |
| Chyu, Hsing, Shih & Natarajan 1998 (98-GT-175) | `Nu/Pr^0.4 = a Re^b` | Han 2012 Table 4.7 (p. 453) | S/D = X/D = 2.5, H/D = 1 (Lyall 2006 Table 2-1); 7 rows | Naphthalene analogy; Han Fig. 4.130 confirms the Total rows. **Re basis (D, Vmax) is assumed**: Lyall labels it Re_d, and Chyu staggered agrees with Metzger 1982a within 5-15% on that basis. Not read from the paper. |

Table 4.7:

| arrangement | surface | a | b |
|---|---|---|---|
| inline | pin | 0.155 | 0.658 |
| inline | endwall | 0.052 | 0.759 |
| inline | total | 0.068 | 0.733 |
| staggered | pin | 0.337 | 0.585 |
| staggered | endwall | 0.315 | 0.582 |
| staggered | total | 0.320 | 0.583 |

**The Chyu inline/staggered modifier** is the ratio of the two Total rows:
`0.2125 Re^0.150`, i.e. 0.76 / 0.85 / 0.94 / 1.00 at Re 5k / 10k / 20k / 30k.
Lyall reports Chyu's combined staggered result as 10-20% above inline.

## Cross-source agreement recorded so far

Read by eye from the figures, not digitised; scoring is PR C.

- **Lawson, Thrift, Thole & Kohli 2011 (Penn State) Fig. 15** vs Metzger 1986: within 5-10%, up to 14% high at Re 3e4 (S2/d 1.73).
- **Uzol & Camci 2005**, 2 rows, endwall: Metzger corrected to 2 rows sits 10% above, with the same slope.
- **Damerow at X_T/D 4.24** vs Metzger 1982b: 2-6% up to Re 1e4, diverging to 22% at 3e4.
- **Han Fig. 4.143 (Metzger 1984, S = X = 2.5)** vs Metzger 1982b: about +/-15%.

## Not carried (declared)

- **Row-count curve:** Metzger 1986 / A&W Fig. 3, to be digitised. The single published point is Nu2/Nu10 = 0.9.
- **Convergence:** Metzger phi = 2.28 Re^-0.096; Brown f multiplier exp(-0.0612 theta); Brigham 1984 attributes the effect to pin height. Unresolved.
- **Long pins:** Faulkner 1971 Eq. D-2, whose pin-height term rests on Kays & London PF-4 only.
- **Pin-endwall fillets:** Chyu 1990.
- **Fillets and inline friction** beyond Chyu (1990)'s single geometry.

## Chyu (1990), J. Heat Transfer 112, 926

Carnegie Mellon, naphthalene. Inline and staggered arrays, H/D 1, S/D = X/D = 2.5, 7 rows, straight and fillet pins.

- **Basis:** `Re = Umax D/nu`, `Umax = Q/A_min` (Eqs 8-9), the canonical basis.
- **Table 2 (printed):** `Nu/Pr^0.4 = A Re^B`, measured on the pins. Table 1: endwall/pin Sh 0.89-1.09.

  | array | pin | A | B |
  |---|---|---|---|
  | inline | straight | 0.463 | 0.537 |
  | inline | fillet | 0.403 | 0.550 |
  | staggered | straight | 0.690 | 0.511 |
  | staggered | fillet | 0.234 | 0.608 |

- **Friction (Fig. 6 only):** `f = 2 dp/(rho Umax^2 N)` (Eq. 16), 4x canonical.
  - Digitised 2026-10-06 into `validation/cooling/data/chyu1990/`; calibration and overlay are recorded in its metadata header.
  - Fits:

    | series | fit | rms |
    |---|---|---|
    | inline straight | constant 0.1693 (slope +0.004) | 1.4% |
    | inline fillet | 0.685 Re^-0.139 | 3 points |
    | staggered straight | 1.616 Re^-0.187 | 2.0% |
    | staggered fillet | 9.387 Re^-0.363 | 1.9% |

- **Cross-check against Gemini's independent digitisation:** where the same marker was read, the two agree within 0.001 in f.
  - The differences, by the calibrated pixel read, are Gemini's: an x offset of about -0.15e4 on the leftmost points, the inline-fillet diamond at 2.21 counted as a circle, and the last staggered-fillet triangle counted as a square. The user's visual check of these points on the scan is pending.
  - A third table supplied as a "digitisation" contains smooth invented values with the inline levels near 0.29; it was discarded.
- **Staggered straight vs Metzger 1982b:** 16-25% lower at the same geometry. Reported, not reconciled.
