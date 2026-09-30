# Extraction: multi-row film superposition (Sellers, and Gao's correction)

**Status: IMPLEMENTED 2026-09-29**, with one deliberate gap: Gao's
correction coefficients ARE published but do not transfer, so `alpha` is a caller input
defaulting to plain Sellers rather than a fitted default.

## Sources, pinned

> Sellers, J.P. (1963). "Gaseous film cooling with multiple slot injection."
> *AIAA Journal* **1**(9), 2154-2156. Not held directly; restated as Eq. (1)
> by Gao et al. below, and independently described in the AFIT thesis
> "Overall Effectiveness Superposition Theory for a Film Cooled Leading
> Edge" (2023), both in `docs/heat_transfer/film/`.

> Gao, Z., Qiu, T., Liu, P., Ding, S., Li, Z., Cheng, R. and Yuan, Q.
> (2025). "A Study on the Film Superposition Method for the Multi-Row Film
> Cooling of the Turbine Outer Ring." *Processes* **13**, 143.
> `docs/heat_transfer/film/processes-13-00143-v2.pdf`. Open access.

## What is implemented

**Sellers, Eq. (1).** The paper prints it as a sum:

    eta = eta_1 + sum_{i=2..n} eta_i prod_{j<i} (1 - eta_j)

which is algebraically `1 - prod_i (1 - eta_i)`. The product form ships,
because it is O(n) and does not accumulate the pairwise rounding of the sum;
a test checks the two against each other on four row counts.

**Gao Eq. (7), the corrected form.**

    eta = sum_i [ eta_i prod_{j=i..n-1} alpha_j prod_{k=i+1..n} (1 - eta_k) ]

`alpha_j` is the correction applied between row `j` and row `j+1`, so there
are `n-1` of them for `n` rows. All ones reproduces Sellers exactly -- a
test asserts that to 1e-14, because getting either product's index range
wrong silently shifts every correction by one row.

**Gao Eq. (5)**, the correction's functional form; **Eq. (9)** equivalent
slot width `s = A_hole/pitch`; **Eq. (10)** equivalent blowing ratio
`M_e = M_0 A_0/A_e`.

Gradients `d eta/d eta_i` for both superposition forms, built from explicit
partial products rather than by dividing the total product -- a fully
effective row (`eta = 1`) would otherwise be 0/0. Tested against finite
differences and specifically at `eta = 1`.

## Why a correction is needed at all

Sellers assumes the rows are independent. Gao measures the cost:

> "the Sellers method accumulates prediction errors as the number of hole
> rows increases, leading to an **overestimation** of the cooling
> efficiency."

That failure is worst exactly where effusion lives -- many closely spaced
rows. A film module built for a few rows and then reused for effusion would
be wrong in the regime it is needed most, which is why #386 and #387 were
designed together rather than in sequence.

The Xia, Chen and Ellis (2024) review reaches the same place from the other
side: Sasaki et al. found superposition from single-hole results works "only
in a far wake region (x > 3D) and only when the holes are laterally widely
spaced". Baldauf already carries the LATERAL interaction as a correlated
parameter, so superposing Baldauf ROWS is on firmer ground than the review's
blanket verdict -- but the STREAMWISE accumulation is still uncorrected, and
that is what Gao's alpha addresses.

## alpha is derived, not invented

Eqs. (3) and (4) are an energy balance on mainstream entrained into the
boundary layer at each injection:

    (T_g - T'_aw)/(T_g - T_aw) = C (m_c/m_g) / (C (m_c/m_g) + 1)

so alpha is the fraction of the film's temperature deficit that survives
mixing on the way to the next row. Three consequences worth keeping:

- **`alpha = 1` is exactly Sellers**, so the identity case is the classical
  model rather than an arbitrary reference point.
- **It is bounded in [0, 1] by construction**, and the implementation
  refuses anything outside -- `alpha > 1` would create coolant.
- **It is per-row**, so it degrades gracefully as rows accumulate, which is
  precisely where the uncorrected model fails.

This replaced an earlier proposal to tune with a coverage exponent
`eta' = 1 - (1 - eta)^k`. That would also have been bounded and monotone,
but it was invented for the purpose; alpha is the source's own parameter
with a derivation behind it.

## What Gao does NOT publish

Eq. (5) gives alpha's form,

    alpha_i = a r / (a r + 1) + b,    r = m_coolant / m_mainstream

**CORRECTION, 2026-09-30: the paper DOES print them.** Section 4.3.1:

> "The empirical coefficients in Equation (5), a and b, were determined to
> be 12 and 0.9465, respectively."

This record previously said they were never printed and that the search had
been exhaustive. It was not. The search looked for `a = 12`; the paper
states the value in prose, in a sentence naming the coefficients rather than
assigning them. A grep over a maths paper's prose is not an exhaustive
search for a constant, and "searched exhaustively" should not have been
written on the strength of one.

They are recorded in code as `film_superposition::gao_a_case1` and
`gao_b_case1`.

**They are still not defaults, and the reason needs no rig arithmetic.**
Since `a r / (a r + 1) >= 0` for `a > 0` and `r >= 0`,

    alpha >= b = 0.9465    for EVERY r

so the published pair can damp a row by at most 5.35% while staying at or
below 1 -- and above 1 the superposition refuses, because `alpha > 1` would
create coolant. Murray's 5.75 D staggered plate needs a per-row alpha of
about 0.85 at low blowing and 0.69 near M = 1. **Both are below the floor,
so no choice of `r` reaches them.**

That argument is deliberately scale-free. A first version of this note
claimed the coefficients give `alpha > 1` above M = 0.20 on Murray's rig,
computed from an assumed `r = M x A_holes_row / A_duct`. It was withdrawn:
Gao's test-section dimensions are not stated in the extractable text, so
`r` cannot be placed on their scale, and a claim resting on that is
unverifiable. Note the asymmetry -- the earlier *fit* was legitimately
scale-invariant, because `r -> k r` is absorbed by `a -> a/k`, but applying
PUBLISHED coefficients is not.

The likely cause is structural: Eq. (5) carries no streamwise-spacing term,
while Murray shows streamwise spacing dominates -- tripling it to 17.25 D
restores plain superposition to under 10%. Gao calibrated at 10.5 d
streamwise spacing; Murray's effective spacing is 2.875 D, 3.7x tighter.
Two plates with the same coolant fraction and different row spacing get the
same alpha from this form.

One thing checked and found NOT to be the issue: Gao's equivalent blowing
ratio `Me = M0 A0/Ae` (Eq. 10) is used to normalise the single-row
efficiency lookup `eta = f(X/(M s))` onto a baseline hole configuration, not
to form alpha's `r`. Alpha's `r` is `m_i/m_g`, the mass-flow ratio from the
Eq. (3) energy balance, and the quoted calibration blowing ratios 0.3 and
1.0 are `M`, not `Me`. combaero needs no `Me` equivalent because Baldauf
takes lateral spacing `s/D` as an explicit argument where Gao's single-row
database is fixed at one spacing.

See `gao_alpha_fit_on_murray.md` for the full comparison, including the
older coupling form that fits Murray better.

## Falsification

Nine perturbations, each applied alone and reverted, all RED:

| # | perturbation | result |
|---|---|---|
| P1 | Sellers sums instead of multiplying | RED |
| P2 | Eq. (7) alpha product off by one row | RED |
| P3 | Eq. (7) dilution starts at `i` not `i+1` | RED |
| P4 | Sellers gradient by division instead of partial products | RED |
| P5 | Eq. (5) offset `b` neutralised | RED |
| P6 | Eq. (9) slot width inverted | RED |
| P7 | Eq. (10) areas swapped | RED |
| P8 | alpha length check weakened | RED |
| P9 | corrected gradient drops its correction sum | RED |

**P8 exposed a real weakness and was fixed rather than accepted.** With the
length guard weakened, the suite did not fail -- it CRASHED, because the
alpha loop then indexed past the end of the vector. gtest reports a dead
process rather than a failed assertion, so a grep for `[  FAILED  ]` finds
nothing and the perturbation looks like a gap. The alpha accesses now use
`.at()`, so a future bypass of the guard raises a diagnosable exception
instead of undefined behaviour, and P8 now produces a clean failure.

## What is not here

- **No fitted `(a, b)`** -- see above. `alpha` defaults to nothing; the
  caller passes 1.0 for plain Sellers.
- **Eq. (11)**, the variable-temperature form with `xi_g,i` and `xi_c,i`,
  is not implemented. It matters when mainstream and coolant temperatures
  vary along the wall, which a network solve can express but the present
  element interface does not yet ask for.
- **No validation against multi-row data yet.** Andrei (2014) and Murray
  (2018) are the targets, and they will also arbitrate Baldauf's Eq. (31)
  discrepancy. Both are held out until the scoring runner exists.
