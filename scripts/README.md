# CombAero Scripts

Utility scripts for code quality, documentation generation, and testing.

## Source register

### sources.py

The register of every local reference document (`validation/sources.yaml`;
the PDFs under `docs/` are gitignored): citation, content, the questions it
was read for, the verdict, and where the repo uses it.

```bash
uv run scripts/sources.py find impingement effusion   # entries matching ALL terms
uv run scripts/sources.py list --status unscreened
uv run scripts/sources.py list --topic film --long
uv run scripts/sources.py list --wanted [--priority now]   # documents we lack, the question each settles
uv run scripts/sources.py check --docs ../../docs     # schema; unregistered local documents
uv run scripts/sources.py new docs/x/paper.pdf        # stub entry to fill in
```

`python/tests/test_source_register.py` keeps the register valid in CI (the
documents themselves are not there, so CI checks the schema and that every
record and code path it names exists).

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

### gui_result_keys.py

Regenerate `gui/frontend/src/utils/resultKeys.generated.json`: the result keys the GUI backend emits per palette kind, with their SI units from `gui/backend/units.py`. It solves one network per kind through `/solve`, plus the variants that switch on optional diagnostics. The Probe and the Live Telemetry read it, and `quantities.test.ts` checks `QUANTITY_CATALOGUE` against it. Rerun it after adding, renaming or removing a result key. The slow test `test_gui_result_keys_sync.py` fails while the file is stale. Do not edit the JSON by hand.

```bash
uv run python scripts/gui_result_keys.py          # rewrite
uv run python scripts/gui_result_keys.py --check  # exit 1 when stale
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

- **Families:** `networks.FAMILIES` holds 14:
  - layered random networks of the basic elements;
  - wall-coupled hot/cold channels;
  - junctions (MPCE, constant-K, ideal, legacy tee);
  - ejectors;
  - impingement arrays, with and without a bypass crossflow;
  - effusion plates and liners;
  - pin-fin channels;
  - bleed stations;
  - series of vortex/connection elements.
  - crossflow chains fed by inclined jets.
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

### check_ascii.py
The `ascii-only` pre-commit hook (#791): fails on any non-ASCII character in
a C++ or Python source, comments, docstrings and code alike. String literals
are exempt only in the files on its `ALLOWLIST` (display units). Every
tracked source passes; there is no exclude list. ASCII spellings: `--`,
`->`, `~=`, `<=`, `+/-`, `*` for a product or unit dot, `^2`, `degC`,
Greek letters by name (`rho`, `mu`, `Delta_P`), `d` for a partial, `sqrt`.
```bash
uv run scripts/check_ascii.py FILE...          # what the hook runs
uv run scripts/check_ascii.py --all --stats    # counts per category/directory
```

### dependency_floors.py
Prints every runtime requirement of combaero and combaero-gui (core
dependencies, the `examples`/`generator` extras, the GUI). The nightly
`dependency-floors` job installs them with `--resolution lowest-direct` and
runs the fast suite and the examples. Dependabot only ever moves the top of a
range. `--check` fails on a requirement with no lower bound, which
lowest-direct would resolve to the first release ever published. A floor is
the lowest version with no known advisory that passes the suite.
```bash
uv run scripts/dependency_floors.py --check
uv run scripts/dependency_floors.py > floors.in && uv pip install --resolution lowest-direct -r floors.in
```

## Quick Reference

| Task | Command |
|------|---------|
| Check C++ style | `./scripts/check-source-style.sh` |
| Check Python style | `./scripts/check-python-style.sh` |
| Auto-fix Python formatting | `./scripts/check-python-style.sh --fix` |
| Check GUI style | `./scripts/check-gui-style.sh` |
| Run static analysis | `./scripts/run-clang-tidy.sh` |
| Regenerate units docs | `uv run scripts/generate_units_md.py` |
| Regenerate the GUI result-key list | `uv run python scripts/gui_result_keys.py` |
| Verify GUI wheel before publish | `./scripts/test-pypi-wheel.sh` |

All scripts return `0` on success and `1` on failure.
