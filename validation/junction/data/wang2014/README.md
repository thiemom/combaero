# Wang 2014/2015 — Compressible combining flow at 45-deg T-junctions

**Source paper:** Wenhui Wang, Zhenhua Lu, Kangyao Deng, Shuan Qu,
*"An experimental study of compressible combining flow at 45 degree T-junctions"*,
Proc IMechE Part C: J Mech Eng Sci 229(9):1600–1610, 2015 (received 2014).

PDF in `docs/junction/wang-et-al-2014-an-experimental-study-of-compressible-combining-flow-at-45-t-junctions.pdf`
(gitignored — copyright).

## Conventions

- **Combining flow only.** Branches 1 (lateral) and 2 (straight) are inflows;
  branch 3 (common) is the outflow.
- `K_13 = (p0_1 - p0_3) / (p0_3 - p_3)`: lateral inlet to common outlet.
- `K_23 = (p0_2 - p0_3) / (p0_3 - p_3)`: straight inlet to common outlet.
- `q = m_1 / m_3` = lateral / total.
- `M_3` = Mach number in the common branch.
- `a = S_c / S_b` = area ratio (common / lateral). Three values tested:
  a in {1, 1.56, 2.44}. Branch angle fixed at 45 deg throughout.

**Mapping to Bassett naming:**
- `K_13` corresponds to Bassett K12 (joining type 6, branch-to-outlet).
- `K_23` corresponds to Bassett K11 (joining type 6, straight-to-outlet).
- But: Wang reports them as compressible-flow K's at non-zero M_3 (Bassett's
  K11/K12 are incompressible). Direct comparison is only valid at M_3 -> 0
  if the formulas extrapolate cleanly.

## Files

- `wang_fig10_<K_id>_a=1_q=<q>_measured.csv`: q-curves vs M_3 at a=1. K_id
  in {K13, K23}; q in {0, 0.2, 0.5, 0.8, 1}.
- `wang_fig11_<K_id>_a=1.56_q=<q>_measured.csv`: same slice at a=1.56.
- `wang_fig12_<K_id>_a=2.44_q=<q>_measured.csv`: same slice at a=2.44.

Tabulated authoritative values are in `validation/junction/models/wang2014.py`
(Tables 1 and 2 of the paper).

## Quality-check cross-references

- Wang Table 1 (q=0.5, a=1): K_13 / K_23 vs M_3. Cross-check the digitized
  Fig 10 files at q=0.5 against these to catch digitization errors.
- Wang Table 2 (M_3=0.5, a=1): K_13 / K_23 vs q. Same idea on the M_3=0.5
  slice.
- A digitization mistake in an earlier version of Fig 11 K_23 was caught
  this way; the recovered values now agree with Table 2 to <1%.
