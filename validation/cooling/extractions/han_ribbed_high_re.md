# Extraction: rib correlations at high Reynolds number

**Status: UNCONFIRMED and INCOMPLETE -- the correlations themselves are not
yet extracted.**

Nothing here may be implemented. Unlike `han_ribbed.md`, this document does not
yet contain a usable correlation: the surrounding text has been extracted but
Figure 4.193 has not, and the coefficients live there. Tracked by #334, under
#339.

---

## For the reviewer

Nothing to judge yet. The blocking need is material, not decisions.

### What I still need sent

- **Figure 4.193**, all three panels. (a) and (b) are the rib profiles studied;
  **(c) carries the R and G correlations** and is the one that matters
- the coefficients `C1`, `C2`, `m`, `n` and the explicit `e/D` dependence, if
  they are given in text rather than on the figure
- whether the round-edged and sharp-edged cases have separate coefficient sets

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
| 12 | **the coefficients** | `C1`, `C2`, `m`, `n` and the `e/D` dependence | **missing -- Figure 4.193c** |

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
