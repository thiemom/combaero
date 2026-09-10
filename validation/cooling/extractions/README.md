# Correlation extractions

One document per correlation family, recording what was read out of a source
and how it was checked -- before any of it is implemented.

## The gate

An extraction starts at **Status: UNCONFIRMED** and may not be implemented
while it stays there. A human reviews it against the source, then sets the
status and signs the review log. Only then is it treated as correct.

This exists because the correlations removed in #332 were not wrong through
carelessness so much as through a missing gate: a plausible reading became a
constant, the constant became code, and the tests then measured the code
against itself.

## Rules

- **Source material is not committed.** Page scans live in `docs/sources/`,
  which is gitignored -- they are copyrighted. These documents are the tracked
  artifact and must stand on their own.
- **Three states, never two.** *confirmed*, *suspect* (read correctly, content
  doubtful), *missing*. Plausibility never promotes one to another.
- **Quantify what is suspect.** A flagged item carries the magnitude of what
  turns on it, so a reviewer knows what it costs to leave it open.
- **Every extraction needs a check that could have failed**, and the check
  needs an anchor outside the page -- a theoretical structure, a second
  equation, a limiting case. Two expressions from one reading agreeing proves
  nothing.
- **Two extraction channels.** A visual read and `scripts/ocr_page.swift`
  (macOS Vision). Disagreement is the signal; neither is authoritative. OCR
  garbles typeset mathematics but flags itself with low confidence.
- **Escalate disagreements with a high-resolution visual crop**, not with more
  OCR. Re-running OCR on an upscaled crop was measured to get *worse*.
- **Corrections go in the review log**, not silently into the tables. A
  correction is evidence about the extraction process.

## Separating the two targets

Where a source plots both a fitted correlation and the data behind it, they
answer different questions. The **correlation** is what an implementation must
reproduce, and where it is printed as an equation it needs no digitising. The
**scatter** is how well that correlation matched reality, and it is what a
harness tolerance should be set from. Conflating them is how a tolerance ends
up wide enough to admit a 4-5x error.

## Tracking

Filed under #339. Each correlation family has its own issue: ribbed #334,
pin fin #335, dimples #336, impingement #337, film and effusion #338.
