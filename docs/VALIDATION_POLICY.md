# Validation policy: fidelity, accuracy, and tuning

How every correlation in combaero is judged, and against what. Applies to the
whole `validation/` tree -- cooling, junction, ejector, orifice -- not just to
the source it was written for.

## The three questions, and why they must stay apart

**1. Fidelity -- does the implementation mirror the paper?**
Judged on **the author's own data**. Nothing else can answer it. A miss here
is *our* bug: a transcribed constant, a dropped term, a convention applied on
the wrong basis. Errata in the source may be fixed, but the fix is stated.

**2. Accuracy -- what will the paper's model do on someone else's rig?**
Judged **cross-source**, and **the result must be labelled cross-source**.
That is not a weaker answer, it is the only honest one: agreement with the
curve a correlation was fitted to is a consistency check, not a validation.
A miss here is the *model's* limitation, not ours.

**3. Matching a particular rig -- the user's job, and out of scope.**
`Nu_multiplier` and `f_multiplier` exist and work. Choosing them for a
specific machine, and validating that choice, is the user's work on the
user's data. **The harness never scores, recommends or applies a tuner.**

They fail differently, so one number covering both is neither. "We match the
paper and the paper disagrees with the data" and "we do not match the paper"
must never look the same on a scorecard.

## Rules

- **A series scores fidelity only if it comes from the correlation's own
  paper.** The `after:` field records the original author for reprinted
  figures; use it, not the containing book.
- **Same lab, different study is ACCURACY, not fidelity.** NASA CR-3837 is
  Han's own laboratory but a different contract, rig and blockage from the
  correlation it scores, so it tests the model rather than the transcription.
- **Cross-source results carry the label in the report**, not only in a
  document beside it.
- **A stated uncertainty band must be a MEASUREMENT band** -- the source's
  own experimental uncertainty, or digitisation precision. Never derived
  from the model's error on that series. Guard it: a declared band that
  correlates with measured model RMS is circular.
- **This applies to the CORRELATION's accuracy too, not only the data's.**
  `RibCorrelationSet::accuracy_R`/`accuracy_G` carry an
  `AccuracyProvenance`, and only `Stated` -- the author's own claim -- may
  be used as a band. Of the three shipped sets exactly one qualifies:
  `han_park_1988_angled`'s 10.5%/8.8% were measured by this project
  *through* `evaluate_rib`, and `rallabandi_2009_high_re`'s 6.9% from the
  printed equation against digitised points. They are reported with their
  provenance and judged against by nothing (#389).
- **An absent figure is not a zero.** `Unstated` carries NaN, so an
  accidental read fails every comparison instead of passing as "perfect".
- **A set is not judged on conditions its authors never claimed.** Scored
  error is partitioned in-domain / out-of-domain against the set's own
  validity box and the two are never pooled. The partition is a REPORT,
  not a filter: `han_1988_orthogonal`'s largest fidelity row scores better
  outside its box than inside it, so dropping extrapolated points would
  have discarded the better half.
- **Two sides must be measuring the same quantity, and say so.** A series
  declares what its y quantity is referenced to and a correlation set
  declares what it produces (`validation/cooling/conventions.py`).
  Scoring matches them, applies a REGISTERED conversion, or refuses. No
  range check can catch this: every coordinate is in range, the quantity
  is simply not the same one.
- **A registered conversion is geometry or algebra, never a fit.** Rohde's
  `Cd` on duct total against McGreehan-Schotsch's on static is
  `sqrt(VHR/(VHR-1))` -- definitional, and registered. Han's
  `G_bar = 1.2 G` looks the same shape and is NOT registered: it is a
  published constant whose measured ratio spans 1.096-1.413, so it is
  applied as the source's own correlation with its cost measured.
- **Fidelity needs completely-sampled data to mean anything.** See below.

## The trap: sampling completeness cuts across both

A figure that overplots several symbol classes yields only its spatially
isolated marks to a digitisation, and those are the ones furthest from the
cluster centre. Such a series' error is an **upper bound**, not an estimate
(see `validation/cooling/recovery.py` and #393).

That interacts with fidelity in a way a single number hides:

| set | judged on | N | MAE | within |
|---|---|---|---|---|
| `han_1988_orthogonal` | FIDELITY / partial | 98 | 5.3% | 65.3% |
| `han_1988_orthogonal` | **FIDELITY / complete** | **8** | 4.4% | 62.5% |
| `han_1988_orthogonal` | ACCURACY / complete | 46 | 4.3% | 78.3% |
| `han_park_1988_angled` | FIDELITY / complete | 154 | 7.4% | 66.2% |
| `rallabandi_2009_high_re` | FIDELITY / complete | 33 | 5.1% | 69.7% |

**`han_1988_orthogonal`'s fidelity rests on eight cleanly-sampled points.**
Its own data survives only as an overplotted scatter, so how faithfully it is
implemented is barely measurable -- and its cross-source accuracy looks
*better* than its fidelity purely because the cross-source data is tabulated.
Report that, do not average it away.

So the honest fidelity number is **own data AND complete sampling**. Where
that set is thin, say so.

## Worked example: what this bought

NASA CR-3837 against `han_park_1988_angled` (#398) is `e/D` = 0.063, a
blockage Han and Park never measured -- accuracy, cross-source, labelled:

- **Eq. 4.17 (`R`) confirmed**: 100% within band at 90 deg, 83-86% at
  75/60/45, inside the source's own 6.6% friction uncertainty.
- **Eq. 4.18 (`G`) reads low at every angle** (-0.8% to -19.9%).

Those are two different verdicts about the same paper, and neither is a
statement about combaero's transcription. Reported as one pooled MAE they
would have been indistinguishable from an implementation bug.

## Reading the scorecard

`validation/cooling/scorecard.py` renders `sampling` (`complete` / `partial` /
`unknown`) per row and segregates the rollup by it. `partial` rows are an
upper bound on the error. `unknown` means the series shares a panel and
nothing independent says how many runs it should have -- neither an estimate
nor a bound.

A row tagged `{out}` was scored **outside the correlation set's own declared
validity box**. Those never pool with in-domain rows. Read them as what
extrapolating that set costs, not as its accuracy:
`baldauf_2002_sellers` scores 138 points of which all 138 are outside its
envelope, so its 48% MAE is entirely an extrapolation figure.

**`Fidelity: checks against quantities the source PRINTS`** answers
question 1 -- does the implementation mirror the paper? A check here
compares against something the source prints: a constant, a table entry,
an equation the author evaluates himself, or an identity between two
printed equations. Scoring against a digitised figure is NOT one: that is
the fidelity/accuracy basis above, and it cannot separate "we transcribed
it wrong" from "the model misses". A printed quantity can -- which is how
Andrews' -13.5% became a model limitation rather than our bug.

**The count is evidence, not a score.** A set with more checks is better
evidenced, not more accurate; `baldauf_2002_sellers` carries one whose
RESULT is that the paper contradicts itself. Sets with no check at all
are named rather than left absent, because "no check" and "nobody looked"
both read as zero.

**`Measurement conventions, series against the set scoring it`** counts
`direct` / `converted` / `undeclared` and lists every conversion applied and
every refusal. A refusal is a category error, not a disagreement -- it
replaces the numbers rather than appearing beside them, because a large
error reported for a mismatched quantity reads as a model limitation.

The report ends with **`Correlation-set accuracy, as each set declares it`**,
which prints each set's `accuracy_R`/`accuracy_G` with its provenance and
says which are usable as bands. A `measured` figure is the error of the model
it describes; nothing is judged against it.
