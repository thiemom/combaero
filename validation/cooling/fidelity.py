"""Question 1, which had no column: does the implementation mirror the paper?

#389 names three questions and the harness answered two of them. Fidelity
was "partly built and entirely invisible" -- the checks existed, scattered
across unit tests and extraction records, and nothing reported "this
implementation was checked against N printed quantities from its source".

WHAT COUNTS AS A FIDELITY CHECK HERE. A comparison against something the
SOURCE PRINTS: a constant, a table entry, an equation the author
evaluates himself, or an algebraic identity between two printed
equations. Scoring against a digitised figure is NOT one of these -- that
is the accuracy/fidelity basis the scorecard already reports, and it
cannot separate "we transcribed it wrong" from "the model misses". A
printed quantity can.

WHY A DECLARED REGISTRY AND NOT DISCOVERY. #389 called this "mostly
aggregation of checks that already run", and the machine-countable part
is thin: `verify.py` finds 9 `printed-curve`/`printed-exponent` findings
across the whole dataset and only one lands on a scored set. The
substantial checks are statements about a paper -- "Eq. (19)'s printed
0.27 equals 0.881 X/(pi L)" -- and no amount of test-name parsing
recovers that. So each is declared with the quantity, the agreement and
the test that pins it, and `test_every_declared_check_exists` keeps the
registry from going stale.

THE COUNT IS NOT A SCORE. A set with many checks is better evidenced,
not more accurate. `baldauf_2002_sellers` carries a check whose RESULT is
that the paper contradicts itself -- that is fidelity evidence of the
most useful kind, and it would be perverse for it to read as a defect.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FidelityCheck:
    """One comparison against a quantity the source prints."""

    #: Correlation set this is evidence about.
    correlation_set: str
    #: What the source prints, specifically enough to find on the page.
    printed: str
    #: Where in the source.
    location: str
    #: What agreement was reached, in the units of the thing compared.
    agreement: str
    #: The test that pins it: a pytest node id, or "ctest:<suite>.<name>".
    pinned_by: str


CHECKS: tuple[FidelityCheck, ...] = (
    # --- Andrews 86-GT-225, the effusion internal correlations -----------
    FidelityCheck(
        "andrews_1986_effusion_internal",
        "Eq. (19)'s printed leading constant 0.27",
        "86-GT-225 Eq. (19), for X = 6.11 mm and L = 6.35 mm",
        "0.26983 against 0.27 -- the only place the paper evaluates "
        "Eq. (18)'s geometry factor numerically",
        "python/tests/test_effusion_internal_runner.py::"
        "test_table_5_geometry_is_internally_consistent",
    ),
    FidelityCheck(
        "andrews_1986_effusion_internal",
        "Eqs. (13) and (14) must agree where they meet",
        "86-GT-225, the L/D = 2 branch junction",
        "2.3600 against 2.3588, i.e. 1.3e-3 -- two independent polynomials "
        "in different variables, which a misread coefficient would not "
        "reproduce",
        "ctest:EffusionInternal.MillsBranchesAgreeWhereTheyMeet",
    ),
    FidelityCheck(
        "andrews_1986_effusion_internal",
        "the author's own evaluation of his own Eq. (19)",
        "86-GT-225 Fig. 10, squares",
        "1.5% mean, 2.5% worst -- digitisation precision. This is what "
        "turned the -13.5% against Fig. 8 from an inference into a "
        "conclusion",
        "python/tests/test_effusion_internal_runner.py::"
        "test_fig10_eq19_series_reproduces_our_summed_correlation",
    ),
    FidelityCheck(
        "andrews_1986_effusion_internal",
        "the author's own evaluation of Eqs. (12)-(14), the throat alone",
        "86-GT-225 Fig. 10, triangles",
        "0.8% mean, 1.6% worst -- the paper effectively plots "
        "mills_entry_length_factor",
        "python/tests/test_effusion_internal_runner.py::"
        "test_fig10_throat_series_reproduces_our_entry_length_factor",
    ),
    FidelityCheck(
        "andrews_1986_effusion_internal",
        "Table 5's printed L/D column, against D recovered from X/D",
        "86-GT-225 Table 5, four plates",
        "1% on every plate, and the four L/D land on Fig. 10's abscissae "
        "to 1% -- a different figure read in a different session",
        "python/tests/test_effusion_internal_runner.py::"
        "test_table_5_geometry_is_internally_consistent",
    ),
    FidelityCheck(
        "andrews_1986_effusion_internal",
        "Table 1's printed A/A_h and hole count N",
        "88-GT-290 Table 1, plates A, B and C",
        "A/A_h within 2% from (X^2 - pi D^2/4)/(pi D L); N = 4306 against "
        "1/X^2 = 4305.8 at X = 0.6 inch",
        "python/tests/test_effusion_plate_element.py::"
        "TestGeometryAgainstAndrewsTable1::test_hole_density_and_area_ratio",
    ),
    FidelityCheck(
        "andrews_1986_effusion_internal",
        "Table 3's printed relative cooling effectiveness",
        "88-GT-290 Table 3, nine cells at G = 0.2, 0.5 and 1.0",
        "0.57 percentage points mean, 2.0 worst -- which confirms the "
        "READING of the table as well as the digitisation, since the "
        "paper never says what the ratio is",
        "python/tests/test_effusion_overall_runner.py::"
        "test_fig10_reproduces_the_papers_own_table_3",
    ),
    # --- Baldauf 2002 -----------------------------------------------------
    FidelityCheck(
        "baldauf_2002_sellers",
        "Table 4's 19 fitted coefficients, from Table 3's inputs",
        "Baldauf et al. (2002), Tables 3 and 4",
        "17 of 19 to better than 2e-6 relative",
        "ctest:BaldaufTable4.SeventeenCoefficientsReproduceThePapersWorkedExample",
    ),
    FidelityCheck(
        "baldauf_2002_sellers",
        "Eq. (31) against the paper's own Table 4 value for b_0",
        "Baldauf et al. (2002), Eq. (31) and Table 4",
        "DISAGREES by 36% (0.83612 against 0.61626) and the equation is "
        "implemented AS PRINTED. A fidelity check whose result is that "
        "the source contradicts itself -- recorded so the next reader "
        "meets it immediately",
        "ctest:BaldaufTable4.EquationThirtyOneContradictsThePapersOwnTable",
    ),
    # --- Han 1988 ---------------------------------------------------------
    FidelityCheck(
        "han_1988_orthogonal",
        "the drawn 90 degree correlation line, stated to BE Han (1988)",
        "Han's textbook figure 4.54, via Han and Zhang (1992)",
        "inside Han's own stated 6% band -- the set meeting a redrawing "
        "of itself in a later paper",
        "python/tests/test_cooling_validation.py::"
        "test_han_1988_reproduces_its_own_figure_454_line",
    ),
)


def machine_counted(dataset) -> dict[str, int]:
    """`printed-curve` and `printed-exponent` findings, per correlation set.

    The part #389 expected to be pure aggregation. It is real but small --
    most such findings sit on `scores: null` correlation curves, which
    belong to no set.
    """
    from collections import Counter

    from validation.cooling import verify

    out: Counter[str] = Counter()
    for series in dataset:
        if not series.scores:
            continue
        for finding in verify.check_series(series):
            if finding.check in ("printed-curve", "printed-exponent"):
                out[series.scores] += 1
    return dict(out)


def count_by_set(dataset=None) -> dict[str, int]:
    """Declared checks per set, plus the machine-counted findings."""
    from collections import Counter

    out: Counter[str] = Counter(c.correlation_set for c in CHECKS)
    if dataset is not None:
        for name, n in machine_counted(dataset).items():
            out[name] += n
    return dict(out)


def sets_without_checks(dataset) -> list[str]:
    """Scoring sets with NO printed-quantity check of any kind.

    Derived, never declared. A hard-coded list said
    `han_park_1988_angled` had none while the machine count found one --
    the report contradicted itself in adjacent lines, which is how a
    stale list fails.
    """
    scoring = {s.scores for s in dataset if s.scores}
    counts = count_by_set(dataset)
    return sorted(name for name in scoring if not counts.get(name))


def render(dataset=None) -> list[str]:
    """The column question 1 never had."""
    counts = count_by_set(dataset)
    if not counts:
        return []
    declared: dict[str, list[FidelityCheck]] = {}
    for c in CHECKS:
        declared.setdefault(c.correlation_set, []).append(c)
    machine = machine_counted(dataset) if dataset is not None else {}

    out = ["Fidelity: checks against quantities the source PRINTS"]
    for name in sorted(counts):
        extra = machine.get(name, 0)
        suffix = f"  ({extra} from figure cards)" if extra else ""
        out.append(f"  {name:<40} {counts[name]:>3}{suffix}")
        for c in declared.get(name, ()):
            out.append(f"      {c.printed}")
            out.append(f"        {c.location} -- {c.agreement}")
    out.append(
        "  A count is EVIDENCE, not a score: a set with more checks is "
        "better evidenced, not more accurate."
    )
    if dataset is not None:
        missing = sets_without_checks(dataset)
        if missing:
            out.append(
                "  No printed-quantity check at all: " + ", ".join(missing)
            )
            out.append(
                "  Named rather than left absent -- 'no check' and 'nobody "
                "looked' both read as zero."
            )
    return out
