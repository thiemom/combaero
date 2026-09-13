# Extraction: rib correlations at high Reynolds number

**Status: UNCONFIRMED -- Figure 4.193 now extracted, one item unresolved.**

Nothing here may be implemented. The correlations are now present, but the
`(e/D)` term of the new correlation fails a check against the figure it is
printed on -- see item 20. Tracked by #334, under #339.

---

## For the reviewer

Nothing to judge yet. The blocking need is material, not decisions.

### For the reviewer

| # | question | why it matters |
|---|---|---|
| **20** | Does the new correlation's `(e/D)` term carry a **negative** exponent, or is `e/D` expressed as a **percentage** there? As literally transcribed it contradicts the figure it is printed on | a factor of 1.9 in `G` at `e+ = 18,000`, and it inverts which curve lies above the other |
| **21** | Is the `(p/e)` exponent `-0.021` or `-0.031`? | 1.6-2.3% over the stated `p/e = 5-10` |
| **22** | Is the `e+` exponent `0.42`? | read consistently, but the digits are small |
| **23** | Is there a separate **friction** (`R`) correlation for the extended range? | the text says both `R` and `G` were modified, but only `G` appears on Fig. 4.193c |
| **24** | Are there separate coefficient sets for round-edged vs sharp-edged ribs, or does round-edged simply use the old correlation? | the text implies the latter; confirmation would close the selection rule |

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
but not conclusive, and it is the reviewer's call.



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

## Review log

| date | reviewer | outcome |
|---|---|---|
| 2026-09-10 | extracted by Claude | UNCONFIRMED and INCOMPLETE. Surrounding text extracted; Figure 4.193 needed before anything here can be used |
| 2026-09-13 | reviewer | supplied Figure 4.193 (p. 513) and Figure 4.192. Correlations now extracted as items 12-19 |
| 2026-09-13 | reviewer | read the Fig. 4.193c exponents as `(p/e/10)^0.2` (old) and `(p/e)^-0.031` (new), both flagged as partly guesswork. The old-correlation exponent is settled at `0.1` by the confirmed p. 377 items; the new one is recorded as item 21, bounded at 1.6-2.3% |
| 2026-09-13 | extracted by Claude | **Item 20 raised.** The new correlation's `(e/D)^0.14`, as transcribed, puts curve 2 below curve 1 everywhere, where the figure shows the opposite. A negative exponent or a percentage convention both resolve it and are numerically identical. Not resolved here |
