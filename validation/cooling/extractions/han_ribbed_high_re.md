# Extraction: rib correlations at high Reynolds number

**Status: CONFIRMED -- reviewed 2026-09-14.**

Every item is resolved against the primary source and the reviewer has signed
off. The correlations below come from **the paper itself**, not from the
textbook figure, which misprints one exponent. Released for implementation
under #334.

Changes after this point go in the review log, not silently into the tables.

Tracked by #334, under #339.

---

## For the reviewer

Nothing to judge yet. The blocking need is material, not decisions.

### For the reviewer

All items resolved. Two things to sign off:

| # | what | note |
|---|---|---|
| **33** | the correlations below are taken from **Rallabandi, Yang and Han (2009)** directly, in preference to the textbook figure | ACCEPTED 2026-09-14 |
| **D4** | implement the paper's values; record the textbook discrepancy | ACCEPTED 2026-09-14 |

### D4. Follow the primary source where it disagrees with the spine

**Decision, accepted 2026-09-14.** Implement Rallabandi, Yang and Han (2009)
Eq. (17) and (18) as printed in the paper, not as reprinted in Han, Dutta and
Ekkad Figure 4.193c.

**Why.** Three independent reasons, and the first is sufficient on its own:

1. the paper is the primary source and the textbook is reprinting it
2. the paper's form lands at 0.994 against the digitised data; the textbook's
   at 0.781
3. the paper is a **cleaner print**. The textbook figure's exponents are small
   superscripts on a scanned page, which is where a lost decimal place is most
   likely to originate and hardest to catch

**Scope.** This decision is about this one disagreement, not a general demotion
of the spine source. Where the textbook is the only source, it remains the
source; `han_ribbed.md` is unaffected and was independently corroborated by
this paper's Eq. (12) and (13).

**Recorded for the next reader.** Anyone comparing an implementation against
Figure 4.193c will find a 22% discrepancy in `G` and should not treat it as an
implementation bug.

### The correlations, from the primary source

> Rallabandi, A.P., Yang, H. and Han, J.-C. (2009). Heat Transfer and Pressure
> Drop Correlations for Square Channels With 45 Deg Ribs at High Reynolds
> Numbers. *ASME J. Heat Transfer*, **131**(7), 071703.
> doi:10.1115/1.3090818

```
R = 1.13 (e/D)^-0.17 (p/e)^0.38                          Eq. (17)
G = 1.24 (e/D)^0.014 (p/e)^-0.02  (e+)^0.42              Eq. (18)

30,000 < Re < 400,000,   0.1 < e/D < 0.18,   5 < p/e < 10
```

with `R` and `G` defined exactly as in the confirmed document:

```
R(e+)     = (2/f)^(1/2) + 2.5 ln( (2e/D)(2W/(W+H)) ) + 2.5       Eq. (12)
G(e+,Pr)  = R(e+) + [ f/(2 St_r) - 1 ] / (f/2)^(1/2)             Eq. (13)
```

**Eq. (12) and (13) are character-for-character Han's Eq. 4.15 and 4.16.** The
paper and the textbook agree on the framework and disagree only on one printed
digit -- which is further corroboration of `han_ribbed.md`, from the primary
source this time.

### Item 19/23 closed

The friction correlation absent from Figure 4.193 is Eq. (17). The extended
range now has both `R` and `G`.

---

## Source

> Han, J., Dutta, S. and Ekkad, S. (2012). *Gas Turbine Heat Transfer and
> Cooling Technology*. 2nd Edition. CRC Press.

after:

> Rallabandi, A.P. et al., *ASME J. Heat Transfer*, **131**(7), 071703, 2009a.
>
> Rallabandi, A.P. et al., Heat transfer and pressure drop measurements for a
> square channel with 45 deg round edged ribs at high Reynolds numbers,
> *Proceedings of ASME Turbo-Expo 2009*, Orlando, FL, ASME Paper
> GT2009-59546, June 8-12, 2009b.

**This is a separate regime from `han_ribbed.md`, not a revision of it.**
Different experiments, different parameter range, and -- per the text --
correlations that do not agree with Han's in the extended range. The two must
be selected between explicitly, never blended. The Fig. 4.46 / Fig. 4.47
episode in the companion document is the precedent for why.

---

## Extracted items

| # | item | as extracted | state |
|---|---|---|---|
| 1 | 2009b scope | stationary square channel, 45 deg square/sharp edged ribs, `Re = 30,000-400,000` | confirmed |
| 2 | why that range | "typical of land-based turbines" | confirmed |
| 3 | 2009a scope | adds 45 deg **round-edged** ribs, to account for manufacturing effects | confirmed |
| 4 | both studies | an array of blockage ratios and rib spacing ratios | confirmed |
| 5 | relationship to Han | the correlations of Han (1988), Han et al. (1988) and Han and Park (1988) were **modified** to fit the extended parameter range | confirmed |
| 6 | typical correlation forms | `R = C1 (P/e)^m` and `G = C2 (P/e)^m (e+)^n` | confirmed |
| 7 | what changed | in these studies the blockage ratio `e/D` **had to be explicitly included** in `R` and `G` | confirmed |
| 8 | `e+` range extension | from `e+ = 1,000` (`Re = 70K`, `e/D = 0.078`) to `e+ = 18,000` (`Re = 400K`, `e/D = 0.18`) | confirmed |
| 9 | `e+` definition | `e+ = (e/D)(Re)(f/2)^(1/2)` -- unchanged from Han | confirmed |
| 10 | disagreement | Fig. 4.193c indicates `R` and `G` **do not agree** with the earlier published correlations in the extended range | confirmed |
| 11 | attributed cause | parameter range differences, specifically `e/D` considerably larger than in prior work | confirmed |
| 12 | old correlation, as printed on Fig. 4.193c | `G = 2.24 (W/H)^0.1 (alpha/90)^0.35 (p/e/10)^0.1 (e+)^0.35`, `P/e = 10-20`, `e/D = 0.047-0.078` | confirmed |
| 13 | **new correlation**, as printed on Fig. 4.193c | `G = 1.24 (e/D)^0.14 (p/e)^-0.021 (e+)^0.42`, `P/e = 5-10`, `e/D = 0.1-0.18` | **see item 20** |
| 14 | new correlation `p/e` range | `P/e = 5-10` | confirmed |
| 15 | new correlation `e/D` range | `e/D = 0.1-0.18` | confirmed |
| 16 | data sets plotted | `e/d` = 0.1, 0.15, 0.18 crossed with `p/e` = 5, 7.5, 10, plus 1984 and 1986 reference data | confirmed |
| 17 | plotted axes | `e+` from `10^2` to `2x10^4`; `G` log axis, ~550 px/decade measured | confirmed |
| 18 | curve behaviour | curves 1 and 2 near-coincident at low `e+`; curve **2 (new, solid) rises ABOVE curve 1 (old, dashed)** at high `e+` | confirmed |
| 19 | extended-range friction (`R`) correlation | not present on Fig. 4.193c | **missing** |

### Item 12 independently confirms the CONFIRMED extraction

The "old correlation" printed on Figure 4.193c (p. 513) is character-for-
character Eq. 4.18 from p. 377, with the square-channel exponents substituted:

```
Fig. 4.193c :  G = 2.24 (W/H)^0.1 (alpha/90)^0.35 (p/e/10)^0.1 (e+)^0.35
Eq. 4.18    :  G = 2.24 (W/H)^0.1 (alpha/90)^m    (p/e/10)^n   (e+)^0.35
               square channel: m = 0.35, n = 0.1
```

**136 pages apart, in different chapters, by different routes into this
extraction.** That is the strongest corroboration `han_ribbed.md` has received.

It also settles a reading dispute. A reviewer read the `(p/e/10)` exponent on
Fig. 4.193c as `0.2`; this extraction read `0.1`. Both were flagged as partly
guesswork. The confirmed items on p. 377 -- established from a page image, the
running text, and Figure 4.47 independently -- give `n = 0.1`. The dispute is
resolved by material already in hand, not by re-reading a blurred superscript.
For scale, the difference reaches 7.2% at `p/e = 20`.

### Item 20: the new correlation's `(e/D)` term fails its own figure

| # | item | state |
|---|---|---|
| 20 | `(e/D)^0.14` as literally transcribed | **contradicts Fig. 4.193c** |

Evaluating both printed correlations at the figure's own conditions -- square
channel, 45 deg ribs, `p/e = 10`, `e/D = 0.1`:

| `e+` | old (curve 1) | new, as transcribed | ratio | new, with `(e/D)^-0.14` | ratio |
|---|---|---|---|---|---|
| 100 | 8.81 | 5.92 | 0.67 | 11.28 | 1.28 |
| 1,000 | 19.72 | 15.57 | 0.79 | 29.68 | 1.51 |
| 18,000 | 54.23 | 52.44 | 0.97 | 99.92 | 1.84 |

**As transcribed, curve 2 sits below curve 1 everywhere.** Item 18 records the
opposite: the figure shows them near-coincident at low `e+` with curve 2 rising
above. On the measured 550 px/decade axis the transcribed reading would put a
**95 px gap** at low `e+` -- roughly a quarter of the plot height, impossible to
miss. The curves visibly touch there.

Two readings remove the contradiction, and they are numerically identical
because `0.1^-0.14 = 10^0.14`:

1. the exponent is **negative**: `(e/D)^-0.14`
2. `e/D` is expressed as a **percentage** in this correlation: `(10)^0.14` for
   `e/D = 0.1`

Either reproduces the figure. **Not resolved here.** Choosing between two
hypotheses that fit equally well, on the strength of which looks more natural,
is the move this process exists to prevent -- and the companion document
records what happened the last time a well-fitting inference was preferred to a
printed statement (item 10 there, settled at 0.51% and wrong).

Note the minus sign on `(p/e)^-0.021` **is** legible in the same equation at
the same size, which is evidence against a simply-missed minus on `(e/D)` --
but not conclusive.

### Digitised curves (items 25-27): curve 1 confirmed, curve 2 excludes both readings

Both plotted correlation lines were digitised, five points each, spanning
`e+ = 85` to `31,000`:
`docs/heat_transfer/han/fig_4.193_correlation_{1,2}_data.csv`.

| # | item | state |
|---|---|---|
| 25 | curve 1 fits `G = 1.696 (e+)^0.35` | confirmed |
| 26 | curve 1 matches the OLD correlation at `alpha = 45 deg`, `p/e = 10`, `W/H = 1` -- predicted 1.7575 against 1.6960 measured, **3.5%** | confirmed |
| 27 | curve 2 fits `G = 1.190 (e+)^0.42`, which matches **no** reading of the printed new correlation | **unresolved** |

**Item 26 validates three things at once.** The old correlation as read, the
digitisation itself, and `alpha = 45 deg` -- which is what the text states the
experiments used. A 3.5% agreement across a digitised curve is about as good as
this method gets, and it means the digitisation can be trusted for item 27.

**Item 27 excludes both hypotheses from item 20.** With the exponent fixed at
the stated 0.42, the measured prefactor is 1.190, while:

| reading | predicted prefactor range | measured |
|---|---|---|
| literal, `(e/D)^+0.14` | 0.856 - 0.943 | **1.190 -- outside** |
| alternative, `(e/D)^-0.14` (equivalently `e/D` in percent) | 1.502 - 1.655 | **1.190 -- outside** |

Measured sits *between* them. Working backwards, the base constant each reading
would need is 1.619, 0.952 or 0.850 against the printed 1.24 -- each off by
23-31%, and in different directions.

So the earlier framing was too narrow: this is not a choice between two
readings of the `(e/D)` term. **No sign or units convention on `(e/D)` alone
reconciles the printed equation with its own plotted curve.** Either another
digit is misread -- the `1.24`, or the `0.42`, whose free fit comes out at
0.3996 -- or the curve is drawn for parameters outside the stated ranges, or
there is an error in the figure.

### The unclassified data cloud (items 28-30) narrows it further

38 experimental points digitised without classification, spanning `e+ = 85` to
`25,091`: `fig_4.193_all_data_points_no classifcation.csv`.

| # | item | state |
|---|---|---|
| 28 | free fit to the cloud | `G = 1.3942 (e+)^0.3974` | confirmed |
| 29 | curve 2 **is** the fit to the cloud | curve 2 is `1.3763 (e+)^0.3996` -- within **1.3%** on the constant and **0.5%** on the exponent | confirmed |
| 30 | the old correlation genuinely misses the extended-range data | mean residual **-15.0%**, RMS 16.6% | confirmed |

Item 29 settles where the fault lies. The digitised curve 2 reproduces the
cloud with a mean residual of **+0.6%** and RMS 6.0%, so the curve and the data
agree with each other. **The printed equation is what disagrees with both**:
evaluated anywhere in its stated parameter box it gives only 0.73-0.81 of the
measured data -- 20-27% low.

Item 30 independently confirms the text's own claim that the old correlation
does not extend: -15% mean bias, and one-signed, which is the signature of a
correlation applied outside its range rather than of scatter.

### Item 31 RESOLVED: the textbook misprints one exponent

The primary source gives

```
G = 1.24 (e/D)^0.014 (p/e)^-0.02 (e+)^0.42        Rallabandi et al. Eq. (18)
```

against Figure 4.193c's

```
G = 1.24 (e/D)^0.14  (p/e)^-0.021 (e+)^0.42       Han, Dutta & Ekkad, p. 513
```

**The textbook prints `0.14` where the paper gives `0.014`** -- a lost decimal
place, worth 22% in `G`. Everything else matches: the constant, the `e+`
exponent, and the `p/e` exponent to within the reading precision (`-0.02`
against `-0.021`).

Against the 38 digitised points:

| | mean pred/data | RMS |
|---|---|---|
| **paper Eq. (18)** | **0.994** | **6.9%** |
| textbook Fig. 4.193c | 0.781 | 22.6% |

The paper's form lands within 0.6% of the data it was fitted to, which is what
a correctly transcribed correlation should do.

**The paper's own text agrees with the small exponent**: *"The heat transfer
roughness (G), on the other hand, does not vary much with e/D or p/e."*
Exponents of 0.014 and -0.02 are exactly that; 0.14 would not be. The textbook's
surrounding prose -- that "the blockage ratio (e/D) had to be factored in
explicitly" -- holds for `R`, whose exponent is `-0.17`, and overstates the case
for `G`.

### How the candidate scan performed

The scan against the digitised cloud was run before any paper was opened, and
ranked `1.24, 0.014, -0.021, 0.42` first at 0.992. The paper gives
`1.24, 0.014, -0.02, 0.42`.

That is worth recording in both directions. The measurement identified the
right answer from a figure whose own equation was wrong -- but the document
deliberately declined to adopt it, because at the time a second candidate
(`constant 1.64`) fitted equally well and nothing but the source could separate
them. **The restraint and the measurement were both necessary**; either alone
would have been insufficient.

By contrast, the earlier `(e/D)^-0.14` hypothesis -- proposed because it fixed
the curve ordering -- was excluded by the cloud at 2.0-2.6x, and would have been
wrong.

### The V-rib sibling paper, in hindsight

The 2011 V-shaped-rib paper gives `(e/D)^-0.08`, which corroborated "small
exponent" but not its sign or value; the angled-rib answer is `+0.014`. Useful
as circumstantial evidence, and correctly not substituted for the real thing.

### If the extended range is needed before the papers arrive

The digitised curve is a usable empirical stand-in, and should be labelled as
exactly that:

```
G = 1.394 (e+)^0.397          6% RMS against 38 digitised points,
                              e+ = 85 to 25,000
```

This is **a fit to our own digitisation of a figure**, not a published
correlation. It carries **no `e/D` or `p/e` dependence**, because the cloud was
digitised without symbol classification and those cannot be separated from it.
Given that the text's whole point is that `e/D` had to be included explicitly,
that is a material limitation, not a rounding one.

If it is ever used it belongs in the modelling decisions section with those
caveats attached, never in the extracted items.

### What would separate the exponents

The figure's legend carries **nine symbol classes** -- `e/d` of 0.1, 0.15, 0.18
crossed with `p/e` of 5, 7.5, 10. Digitised **by class**, the data would
separate the `e/D` and `p/e` exponents instead of constraining only their
product, and would show directly whether the printed exponents are consistent
with the plotted points. That is the one measurement left that this figure can
still yield.

**This still needs the Rallabandi papers.** Figure 4.193 is exhausted; the
old correlation is unaffected and now corroborated from three directions.




## Why the correlations diverge, per the text

The log-law velocity profile underlying `R` and `G` holds when the surface
roughness is relatively small. Larger rib thickness invalidates it: form drag
from flow separating and reattaching at the rib grows relative to skin
friction. That makes `R` and `G` depend on `e/D` as well as on `e+`.

This is a **stated limitation of the method**, not a fitting artefact, and it
bounds `han_ribbed.md` from outside: that document's correlations are valid
where the log-law assumption is, and `e/D = 0.18` is 2.3x beyond its validated
band of 0.047-0.078.

## The selection rule (item 13)

| # | item | as extracted | state |
|---|---|---|---|
| 13 | which correlation to use in the extended range | **round-edged ribs: the OLD (Han) correlation. Sharp-edged ribs: the NEW (Rallabandi) correlation.** | confirmed |

Supporting observations from the text: with round-edged ribs the friction was
lower and the pressure drop smaller, with friction performance "quite similar
to smaller ribs", while heat-transfer coefficients were similar to sharp-edged
ribs. So round-edged ribs behave, frictionally, like the smaller ribs Han's
correlation was fitted on.

**The text calls this agreement "coincidently".** That word is worth keeping:
the round-edge agreement is reported as fortuitous rather than principled, so
it should be treated as an empirical observation with the stated range, not as
evidence that Han's correlation extends on physical grounds.

### Consequence for implementation

**Rib edge profile becomes a required model parameter.** Sharp versus round
selects between two different correlations above `Re ~ 60,000`. The
implementation removed in #332 had no such parameter, and neither does the
confirmed `han_ribbed.md` -- which is correct for its own range, where the
distinction does not arise.

---

## Cross-check performed

Item 8 gives two anchors, and Han's confirmed chain should reproduce the first
(inside his validity) while failing on the second (which item 10 says is where
the correlations disagree). Running `han_ribbed.md`'s chain:

| anchor | `f` | `e+` predicted | `e+` stated | difference |
|---|---|---|---|---|
| `Re = 70K`, `e/D = 0.078` | 0.0700 | 1022 | 1,000 | **2.2%** |
| `Re = 400K`, `e/D = 0.18` | 0.1889 | 22,126 | 18,000 | **22.9%** |

Both outcomes are the expected ones. `e/D = 0.078` is the top of Han's
validated band and the chain lands within 2.2% of a figure quoted in a
different chapter of the book. `e/D = 0.18` is 2.3x beyond that band and the
chain misses by 23%, which is what items 10 and 11 predict.

This is an **independent check on `han_ribbed.md`**, using material that
document never saw. It is recorded here rather than there because the anchors
come from this section.

---

## Re-digitisation check (2026-09-15)

Figure 4.193c was digitised a second time, independently. The confirmed
correlations are unaffected.

| | first digitisation | re-digitisation |
|---|---|---|
| curve 1 | `1.5669 (e+)^0.3610` | `1.5856 (e+)^0.3617` |
| curve 2 | `1.3763 (e+)^0.3996` | `1.3558 (e+)^0.4024` |
| cloud | `1.3942 (e+)^0.3974` | `1.3717 (e+)^0.4023` |
| Rallabandi Eq. (18) against the cloud | 0.994, 6.9% RMS | **0.973, 6.9% RMS** |

Agreement between two independent digitisations is within 1.5% on the
constants and 1.3% on the exponents, and Eq. (18) reproduces both at the same
RMS. The textbook's printed form remains the outlier.

## Review log

| date | reviewer | outcome |
|---|---|---|
| 2026-09-10 | extracted by Claude | UNCONFIRMED and INCOMPLETE. Surrounding text extracted; Figure 4.193 needed before anything here can be used |
| 2026-09-13 | reviewer | supplied Figure 4.193 (p. 513) and Figure 4.192. Correlations now extracted as items 12-19 |
| 2026-09-13 | reviewer | read the Fig. 4.193c exponents as `(p/e/10)^0.2` (old) and `(p/e)^-0.031` (new), both flagged as partly guesswork. The old-correlation exponent is settled at `0.1` by the confirmed p. 377 items; the new one is recorded as item 21, bounded at 1.6-2.3% |
| 2026-09-15 | reviewer | re-digitised Fig. 4.193c independently. Both curves and the cloud reproduce within 1.5%; Eq. (18) still fits at 6.9% RMS. **CONFIRMED status unaffected.** |
| 2026-09-14 | reviewer | **CONFIRMED.** D4 accepted: follow the paper over the textbook where they disagree, noting the paper is also the cleaner print -- small superscripts on a scanned textbook page being where a lost decimal is most likely to arise. Released for implementation. |
| 2026-09-14 | reviewer | supplied the primary source, Rallabandi, Yang & Han (2009), *JHT* 131(7), 071703. **Item 31 resolved**: the paper's Eq. (18) gives `(e/D)^0.014`; the textbook figure prints `0.14`, a lost decimal worth 22%. Eq. (17) closes the missing `R` correlation. The paper's Eq. (12) and (13) reproduce Han's Eq. 4.15 and 4.16 exactly, corroborating the confirmed document from the primary source. **No open items.** |
| 2026-09-13 | reviewer | supplied a scanned paper as `Rallabandi_2009.pdf`. **It is a different paper**: Alkhamis, Rallabandi & Han (2011), *JHT* 133, 111901, on **V-shaped** ribs -- not the 45 deg angled-rib papers the book cites for Fig. 4.193. Its Eq. (18) and (19) are recorded, and its `(e/D)^-0.08` corroborates the data-driven estimate, but they are not the Fig. 4.193 correlation. |
| 2026-09-13 | reviewer | confirmed against the page that the constant reads **1.24** and the `(e/D)` exponent reads **0.14**. That eliminates the last surviving misreading candidates, so **item 31 is no longer a transcription question**: the printed equation disagrees with its own figure by 22%. |
| 2026-09-13 | reviewer | supplied the unclassified 38-point data cloud. **Item 29**: curve 2 is the fit to the cloud, agreeing to 1.3%. So the curve and data agree and the printed equation disagrees with both, running 20-27% low. **Item 31**: the `(e/D)^-0.14` hypothesis is excluded at 2.0-2.6x; the fault is one digit, either `(e/D)^0.014` or a constant of 1.64. |
| 2026-09-13 | reviewer | supplied digitised curves for both correlation lines. **Item 26**: curve 1 matches the old correlation to 3.5% at `alpha = 45 deg`, `p/e = 10` -- validating the reading, the digitisation, and the stated rib angle together. **Item 27**: curve 2's prefactor is 1.190 where the two candidate readings predict 0.856-0.943 and 1.502-1.655. Both excluded; needs the Rallabandi papers. |
| 2026-09-13 | extracted by Claude | **Item 20 raised.** The new correlation's `(e/D)^0.14`, as transcribed, puts curve 2 below curve 1 everywhere, where the figure shows the opposite. A negative exponent or a percentage convention both resolve it and are numerically identical. Not resolved here |
