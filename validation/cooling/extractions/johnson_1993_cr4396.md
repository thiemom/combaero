# NASA CR-4396 Vol. II (Johnson, Wagner and Steuber, 1993) -- stationary trips

**Status: digitised, user-checked and scored, 2026-10-22 (#452).**
Source: `docs/heat_transfer/CR4396/19940011518.pdf` (gitignored), contract
NAS3-23691, Pratt & Whitney / UTRC. Data: `data/johnson1993/`.

## Why this source matters

Every other rib heat-transfer source in this dataset traces to Texas A&M
(Han's or Lau's laboratory; #403 item 2). This is an industrial rig with
cast-style trips, so it is an **independent-lab accuracy check** on the Han
sets. It is never fidelity, and the harness does not tune to it.

**Circularity gate.** CR-4396 is not among the data any shipped rib set was
fitted on. Its Fig. 10.1 *compares* test 301 with Han et al. (1986). That is
a comparison, not pooling.

## What the report prints

| item | value | where |
|---|---|---|
| passage | square, W/H = 1 (H/w = 1.0) | p. 21 |
| trips | e/d = 0.1, P/e = 10, leading + trailing walls only, **staggered** | p. 9, 21 |
| trip shape | curvilinear (rounded, cast-style) | p. 9, 21 |
| skewed model | 45 deg, both walls parallel in projection (Fig. 3.5) | p. 11, 21 |
| thermal BC | all four walls heated, constant wall temperature, inlet (Tw-Tb)/Tw = 0.13 | p. 18 |
| properties | film temperature (Tw+Tb)/2 | p. 18 |
| area | projected; true area 1.11x projected (normal), 1.15x (skewed) | p. 18 |
| reference | Nu_inf = 0.0176 Re^0.8, Kays and Perkins (1973), constant wall T, "Prandtl number equal to 0.72" | p. 18 |
| uncertainty | Kline-McClintock: +-6% at the inlet to +-30% at the exit (stationary); repeats agree within about half that | p. 18-19 |

**Fig. 3.5's arrangement.** The leading-wall trips (dashed) and the
trailing-wall trips (solid) run the same way in the cross-sectional view,
interleaved. That is Han's "parallel" arrangement with a stagger, not a
"crossed" one.

## The conversion, and why it is not optional

The figures plot Nu/Nu_inf on **Kays and Perkins**. The shipped sets' Nu
ratio is the ribbed-wall Nu over **Dittus-Boelter** 0.023 Re^0.8 Pr^0.4.
Both references scale as Re^0.8, so the conversion is one constant:

    Nu/Nu_DB = (Nu/Nu_KP) * 0.0176 / (0.023 * 0.72^0.4) = 0.873 (Nu/Nu_KP)

The Pr in it is the source's own 0.72, because 0.0176 *is* Kays and Perkins
at 0.72. This is definition, not a fit. It is registered in
`conventions.py` as `("Nu_ratio_ribbed_KP", "Nu_ratio_ribbed_DB")`.

**Skipping it hides the result.** On the raw Kays-Perkins basis the Han sets
read -2.2% (normal) and +3.3% (skewed), which looks like agreement. On a
common basis they over-predict by 12-18% (below). This is the category
error #389 names, and the same size as the ~12% gap #444 criterion 6
warned about.

A second-order assumption sits under every Nu-ratio comparison: the ribbed
and smooth Nu share the same Pr^0.4. Here the ratio is taken at Pr 0.72
against a set evaluated at its own 0.70. If they did not share it, the
difference would be 1.1%.

## Digitisation

400 dpi render of PDF pages 48-49. Scope is the first leg, segments B-D,
both walls, four Re: 48 points.

- **Abscissa.** Each symbol's nominal Re from the legend: square 12,500
  (test 302), circle 25,000 (301), triangle 50,000 (303), diamond 75,000
  (304). The figure's own x is S/d, which identifies the segment.
- **Calibration.** The y axis is log, calibrated per panel on the printed
  ticks. The worst tick residual is 1.5 px, or 0.4% in Nu/Nu_inf.
- **Squares** come from the flat top edge (+19.5 px to centre). **Diamonds**
  come from the bottom tip (-30 px to centre). Both are cross-checked
  against an outline-template match.
- **The diamond offset was measured, not assumed.** One plot diamond has a
  clean upper outline (Fig. 6.2, leading wall, segment R). Its top tip is
  row 1597 and its widest row is 1627, so the centre is tip + 30, and the
  lower half mirrors the upper within 1 px. 1 px is 0.26% in Nu/Nu_inf.
  The legend diamonds are drawn about 15% larger than the plot's and are
  not used.
- **Circles and triangles** come from the template match only. A
  triangle's position is its bounding-box centre. If the plotter centred
  triangles on the centroid, they read about 1.8% high.
- **User review.** Nine points were flagged and checked by the user on a
  scan overlay, all confirmed: two template/anchor disagreements and seven
  weak template matches, mostly Fig. 6.2 triangles inside other outlines.
- **Overplotting.** Four symbols share a cluster about 60 px tall, so the
  score is limited by the read. The edge-anchored half (24 points) is the
  best-located.

## Score (cross-source accuracy, on the Dittus-Boelter basis)

Every point is outside both Han sets' e/D validity (0.1 against
0.047-0.078), so every point is flagged extrapolated.

| figure | set | n | bias | MAE | within +-6% |
|---|---|---|---|---|---|
| 6.1 normal | `han_1988_orthogonal` | 24 | +12.1% | 12.1% | 3/24 |
| 6.2 skewed 45 deg | `han_park_1988_angled` | 24 | +18.4% | 18.4% | 0/24 |

By wall: 6.1 leading +16.4%, trailing +7.7%; 6.2 leading +17.0%, trailing
+19.7%. By Re the bias is flat: 6.1 runs +9.8 to +15.0%, 6.2 runs +15.2 to
+20.3%. So the miss is a level offset, not a Re trend.

**Not scored, recorded:** `rallabandi_2009_high_re` (Texas A&M, high Re).
Its envelope is e/D 0.1-0.18, P/e 5-10, 45 deg, Re >= 30k, so it is the
only shipped set whose validity actually covers the skewed rig. On Fig.
6.2 it gives bias +5.0% and MAE 7.3% over all 24 points. On the 12
in-envelope points (Re 50k and 75k) it gives -0.6% and 3.9%. A series
scores one set in this harness, and the issue scoped the Han R/G sets, so
this is a side measurement rather than a scorecard row.

## Not modelled

Staggered rather than in-line trips, rounded trip profiles, constant wall
temperature, heating on all four walls, and the extra trip area (1.11x and
1.15x the projected area) are all left out. Any of them could carry part of
the offset. None is fitted for.
