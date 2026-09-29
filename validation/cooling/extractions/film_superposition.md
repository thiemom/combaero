# Extraction: multi-row film superposition (Sellers, and Gao's correction)

**Status: IMPLEMENTED 2026-09-29**, with one deliberate gap: Gao's
correction coefficients are not published, so `alpha` is a caller input
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

**but the paper never prints the fitted `a` and `b`.** Searched
exhaustively: no coefficient table (Tables 1-5 are existing models, window
positions, plate geometry, infrared calibration and a deviation comparison),
and no inline value anywhere in the text. The shape is sourced; the
constants are not.

So `a` and `b` are REQUIRED arguments rather than defaulted. A made-up
default would acquire an authority it has not earned, and the validation
policy already reserves this slot: matching a specific rig is the user's
job, and the harness never scores a tuner.

Table 5 cannot fill the gap either -- it benchmarks Gao's corrected model
against Zhang et al. [38], not against uncorrected Sellers, so it does not
even quantify what the correction buys.

**The intended path to a default** is our own data: fit `(a, b)` against
Andrei (2014) and Murray (2018), label it tuned, and state its envelope.
That has not been done yet.

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
