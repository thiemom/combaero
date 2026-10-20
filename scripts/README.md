# CombAero Scripts

Utility scripts for code quality, documentation generation, and testing.

## Code Quality Scripts

### check-source-style.sh

Check C++ source code style and quality.

```bash
./scripts/check-source-style.sh
./scripts/check-source-style.sh --verbose
./scripts/check-source-style.sh --strict
```

**Checks:** Non-ASCII characters, block comments, M_PI without `math_constants.h`.
**Scans:** `src/`, `include/`, `examples/`, `tests/`

### check-python-style.sh

Check Python code style and quality.

```bash
./scripts/check-python-style.sh        # Check
./scripts/check-python-style.sh --fix  # Auto-fix formatting
```

**Checks:** Non-ASCII characters, ruff lint, ruff format, optional mypy.
**Scans:** `python/`, `cantera_validation_tests/`, `thermo_data_generator/`

### check-gui-style.sh

Check frontend code style (Biome).

```bash
./scripts/check-gui-style.sh
```

### run-clang-tidy.sh

Run clang-tidy static analysis (requires a CMake build with `compile_commands.json`).

```bash
./scripts/run-clang-tidy.sh
```

## Documentation Scripts

### generate_units_md.py

Regenerate `docs/UNITS.md` from `include/units_data.h`.

```bash
uv run scripts/generate_units_md.py
```

## Stress Testing

### stress/ (run.py, summarize.py, networks.py)

Solve randomly drawn networks COLD, over and past the ranges the test suite
exercises, and record what happened. The suite builds each element near its
validated operating points; these sweeps found what it could not:
- a spurious zero-flow ejector root reported as converged (#619);
- a NaN start reported as converged (#621);
- misfiled residual rows (#614, and a false positive in the row-kind check, #620);
- a choked ejector regime with no root (#624).

```bash
uv run python scripts/stress/run.py --out tmp/stress/after.jsonl                  # all families, seeds 0:40
uv run python scripts/stress/run.py --families ejector,tee3 --seeds 0:100 --jobs 8 --out tmp/stress/x.jsonl
uv run python scripts/stress/run.py --seeds random:40 --out tmp/stress/fresh.jsonl  # 40 fresh seeds; the range is printed
uv run python scripts/stress/run.py --strict ...                                   # raise at the first finding
uv run python scripts/stress/summarize.py tmp/stress/after.jsonl --compare tmp/stress/before.jsonl
```

- **Families:** `networks.FAMILIES` holds 12:
  - layered random networks of the basic elements;
  - wall-coupled hot/cold channels;
  - junctions (MPCE, constant-K, ideal, legacy tee);
  - ejectors;
  - impingement arrays;
  - effusion plates and liners;
  - pin-fin channels;
  - bleed stations;
  - series of vortex/connection elements.
- **Per solve, the runner records:**
  - the outcome;
  - the caught-exception log;
  - anything that escaped `solve()`;
  - the second-law verdicts and ledger;
  - the `check_row_kinds` verdict at the solution.
- **Reproducible seeds:** a (family, seed) pair rebuilds the same network, so a finding is a repro. Keep each generator's draw order when editing it.
- **Fixed vs fresh seeds:** 0:40 is the regression set, for A/B comparison. `random:N` explores: on its first runs it found a generator refusal and an f=0 n-port draw that 0:40 never made (#642). Its records carry their seeds, so a fresh finding is as reproducible as a fixed one. Don't compare a fresh run against a fixed baseline.
- **A failure is a finding to triage, not automatically a bug.** Some draws have no root by construction.
- **A "converged" result is not yet a physical one.** Look at implausible internal states as well.
- **Timeouts are wall-clock,** so for an A/B comparison run both sides with the same `--jobs`.
- **When to run:** after any change to a residual form, a closure's continuation, or the solver's start/retry logic.

## Testing Scripts

### test-pypi-wheel.sh

Build and verify the `combaero-gui` wheel end-to-end before a PyPI release. Builds
the frontend and both wheels, inspects wheel contents, installs in an isolated venv,
and smoke-tests that the server starts and serves JS assets with the correct MIME type.

```bash
./scripts/test-pypi-wheel.sh                        # Full run
./scripts/test-pypi-wheel.sh --skip-frontend        # Reuse existing gui/frontend/dist/
./scripts/test-pypi-wheel.sh --skip-combaero-build  # Use existing dist/combaero-*.whl
```

## Dependency Scripts

### bulk_merge_safe.sh
Merges safe Dependabot PRs (grouped safe-updates, or same-major bumps) one at a
time. `main` requires up-to-date branches and signed commits (ruleset "Iron"),
so each PR is rebased by Dependabot (`@dependabot rebase`, which re-signs;
a GitHub-side `gh pr update-branch` rebase is unsigned and blocks the merge),
waits for CI, is re-checked (all checks passed or skipped, `CLEAN`, pinned
versions agree at its head) and squash-merged at exactly the verified commit.
```bash
./scripts/bulk_merge_safe.sh --dry-run      # report only
./scripts/bulk_merge_safe.sh 510            # one PR
```

### check_version_sync.py
Fails when a tool pinned in several files (biome: package.json, biome.json,
pre-commit, ci.yml) disagrees. `--ref REF` checks a git ref, e.g. a PR head.

## Quick Reference

| Task | Command |
|------|---------|
| Check C++ style | `./scripts/check-source-style.sh` |
| Check Python style | `./scripts/check-python-style.sh` |
| Auto-fix Python formatting | `./scripts/check-python-style.sh --fix` |
| Check GUI style | `./scripts/check-gui-style.sh` |
| Run static analysis | `./scripts/run-clang-tidy.sh` |
| Regenerate units docs | `uv run scripts/generate_units_md.py` |
| Verify GUI wheel before publish | `./scripts/test-pypi-wheel.sh` |

All scripts return `0` on success and `1` on failure.
