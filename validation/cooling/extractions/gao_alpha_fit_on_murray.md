# Fitting Gao's alpha on Murray -- a negative result

Attempted 2026-09-30, after #425 established that `murray2018` is the first
family whose curves all fall the same side of the model, which is the
prerequisite for fitting a correction that only ever reduces.

**Outcome: Gao's published form for alpha cannot be fitted on this data.**
A single constant alpha per blowing ratio works well; the `r` dependence
Gao prescribes does not describe how that constant moves.

## What was being fitted

Gao Eq. (5), the only part of alpha the paper publishes:

    alpha_i = a r / (a r + 1) + b,    r = m_coolant / m_mainstream

`a` and `b` are never printed (searched exhaustively -- see
`film_superposition.md`), which is why they are required arguments in
combaero with no default. This was the attempt to supply them from data.

For one geometry `r` is the same for every row, so alpha collapses to ONE
number per blowing ratio and the fit is a scan over a single parameter per
series -- no optimiser, no local minima, and the objective's shape is
visible rather than assumed.

## A constant alpha works, and nearly reaches the measurement floor

| M | alpha = 1 (Sellers) | best constant alpha | MAE there |
|---|---|---|---|
| 0.19 | 33.1% | **0.850** | 15.8% |
| 0.48 | 43.6% | **0.864** | 21.4% |
| 0.96 | 112.6% | **0.686** | 23.7% |

At M = 0.96 a single knob takes the error from 113% to 24%. The residual is
close to the paper's own stated 15% experimental uncertainty, so a constant
alpha is not far off the noise floor for this rig.

**The minima are real, not plateaus** -- checked, because a flat objective
would mean the "best" alpha is not a measurement:

    M0p19   0.70=25.4%  0.80=17.7%  0.85=15.8%  0.90=17.5%  0.95=23.6%
    M0p48   0.70=32.7%  0.80=23.5%  0.85=21.5%  0.90=24.3%  0.95=33.0%
    M0p96   0.65=24.3%  0.70=23.8%  0.75=27.9%  0.85=53.4%  0.95=89.9%

M = 0.19 and M = 0.48 share a minimum at 0.85 -- their apparent 0.850 vs
0.864 difference sits inside the flat region and is not a real distinction.
M = 0.96 is sharply different at 0.686.

## Why Gao's form cannot carry that

The sequence of best-fit alphas is **0.85, 0.85, 0.69** -- flat, then
falling. Gao's form is strictly monotone in `r`:

    d/dr [a r / (a r + 1)] = a / (a r + 1)^2

whose sign is `a`'s and never changes. A monotone function cannot reproduce
flat-then-falling, and more fundamentally the form is *increasing* for
`a > 0` -- it says MORE coolant needs LESS correction, while Murray needs
more.

Hold-one-out over a grid of `(a, b)` spanning both signs of `a`:

| held out | fitted | predicted alpha | held-out MAE | achievable |
|---|---|---|---|---|
| M = 0.19 | a = -10, b = +0.99 | 0.941 | 22.0% | 15.8% |
| M = 0.48 | a = -5.62, b = +0.86 | 0.789 | 24.2% | 21.4% |
| M = 0.96 | a = +1e4, b = -0.13 | 0.866 | **58.5%** | 23.7% |

Every hold-out is worse than that series' own constant, and the M = 0.96
case -- the one that matters, where the correction earns its keep -- is
2.5x worse. Note also that two of the three fits choose a NEGATIVE `a`,
which is the sign that makes alpha fall with coolant flow: the fit is
fighting the form's intended direction.

## Two things ruled out before concluding

**It is not a matter of getting `r`'s scale right.** Any sensible definition
of `m_coolant/m_mainstream` -- per-row or total, per-hole area or per-row
area -- differs only by a positive constant for fixed geometry, and Gao's
form absorbs a positive scale entirely into `a`, since `r -> k r` is the
same curve with `a -> a/k`. Demonstrated rather than asserted: refitting
with `r` scaled by 0.1, 1, 10 and 100 gives an identical held-out alpha of
0.866 and MAE of 58.5%, with `a` moving 1e5 -> 1e4 -> 1e3 -> 100 exactly as
the algebra predicts.

**It is not the per-row versus cumulative reading.** alpha is per-row, so
`r` might accumulate downstream as `r_i = i r_row`, making alpha rise with
row index. That is a different model, not a rescaling, so it was fitted
separately: held-out MAE 16.3%, 24.8% and **47.6%**. Better than the
constant-`r` reading on the worst case but still twice the achievable
23.7%, and it still selects a negative `a` (-1, -0.316).

## What is NOT being shipped

**No fitted `(a, b)` default.** The form does not describe the data, so any
constants would be a curve forced through points it cannot pass near. `a`
and `b` stay required arguments.

**No invented replacement.** A monotone-decreasing alpha, or a
blowing-ratio-indexed table, would both fit better -- and both would be
made up for the purpose, which is what the tuned-corrections policy
excludes. The shape has to come from a source.

**The per-series constants above are recorded as a REFERENCE for this rig,
not as defaults.** Matching a specific rig is the user's job and the
harness never applies a tuner. Anyone cooling a 5.75D staggered plate at
M below about 0.5 has a starting point of alpha ~ 0.85 here, clearly
labelled as tuned to Murray's geometry and no one else's.

## What would change the verdict

- Data at a **different geometry** at the same blowing ratios. Every point
  here shares one plate, so `r` and `M` move together and the fit cannot
  tell which one alpha actually responds to. That confound is the single
  biggest weakness of this attempt, and no amount of extra blowing ratios
  on the same plate fixes it.
- Gao's own `(a, b)`, if they ever surface. The form may well be right for
  the rig Gao fitted it on.
- Murray's Figure 9a at 17.25D streamwise, where superposition nearly
  works: a fitted alpha should approach 1 there, and a form that cannot
  produce that limit is wrong regardless of what it does at 5.75D.

Related: `film_superposition.md`,
`murray_ireland_2018_effusion_superposition.md`, and #420's finding that
`andrei2014` cannot support this fit at all.
