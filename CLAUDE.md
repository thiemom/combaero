# CombAero — Claude Code Instructions

## Practical Tips
- Avoid shell variable expansions and use literal paths or printenv/env instead. Run shell expansions once per session and use the abs path after. Shell variable expansions require user approval per open bugs (#29616, #18160). Allow list workarounds are not always reliable.
- There are allow rules for sandboxed commands such as uv run, that are best executed in sandbox mode. Disabling sandbox requires user approval and should only be used when absolutely required.
- `Edit(*.pyproject.toml)`-style tools are globally denied for **any** `pyproject.toml` (not just the root one) in Claude Code's own permission settings. Edit these files via a shell command (`perl -i -pe` or similar) instead.
- `git push`/`git tag push`/`git pull` inside the sandbox reliably print `fatal: failed to store: 100001` / `could not lock config file .git/config: Operation not permitted` because writing the upstream-tracking config is a denied path. The underlying git operation still succeeds — verify with `git log`/`gh run list` rather than treating the message as a failure.

## Hard Rules
- **Virtual environment only:** use `uv run` and `uv pip`. Never install into system Python or use global `pip`.
- **Never bypass pre-commit hooks** or CI checks. All code must pass cleanly before pushing.
- **No manual edits to auto-generated files** (`docs/UNITS.md`, `include/thermo_transport_data.h`). Use the scripts below.
- **API sync is mandatory:** adding/removing/changing any function or property requires updating `include/units_data.h` and the relevant API reference ([docs/API_CPP.md](docs/API_CPP.md) or [docs/API_PYTHON.md](docs/API_PYTHON.md)) in the same commit.
- **CHANGELOG is mandatory:** every user-visible change must have an entry added under `[Unreleased]` in `CHANGELOG.md` in the same commit. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) (Added / Changed / Fixed / Removed). At release time the `[Unreleased]` block becomes `[x.y.z] - YYYY-MM-DD`.
- **Documentation Hygiene:** Follow [docs/HYGIENE.md](docs/HYGIENE.md) classification. Archive feature reports to `docs/archive/`.
- **Unit sync test:** every exported Python symbol needs a `units_data.h` entry, or an `IGNORE_LIST` entry in `python/tests/test_units_sync.py`.
- **Solver (f, J) rule:** solver-facing calculations must expose a C++ PyBind11 API returning `std::tuple<double, double>` (value, derivative). Analytical derivatives via chain rule; no finite differences exposed to Python.
- **Define Once:** physics constants go as `constexpr` in the relevant public header namespace — never as magic numbers in `.cpp` files.
- **ASCII only:** no non-ASCII characters in any C++ or Python source file.
- **`python/combaero/` must never import from `validation/`, `cantera_validation_tests/`, `thermo_data_generator/`, or `gui/`.** Those trees are dev-only and excluded from the sdist/wheel (see `pyproject.toml`'s `sdist.exclude`); a production module importing from any of them installs fine from a repo checkout but raises `ModuleNotFoundError` on a real `pip install` (hit in the v0.4.0 release for `MPCEv2Element`/`ConstantKTeeElement` importing `validation.junction.models.mynard2010` — fixed by moving the implementation into `python/combaero/network/_mynard2010.py`). If a model implementation is needed by production code, it belongs in `python/combaero/`, not the validation tree.

## C++ Style (C++17)
- `#pragma once`, sorted minimal includes, `//` comments only (no `/* */`).
- Smart pointers, RAII. No raw `new`/`delete`.
- Use `math_constants.h` for `M_PI` (MSVC compatibility).
- Explicit standard library includes — macOS-only implicit includes break Linux CI.

## Python Style
- `ruff` for lint + format. Run `./scripts/check-python-style.sh --fix` before committing.
- Type annotations on all signatures; `-> None` for void functions; `|` union syntax (Python 3.12).
- No non-ASCII characters, including docstrings and comments.

## Frontend Style (gui/frontend)
- Package manager: `pnpm`. Use `pnpm install`, `pnpm run dev`, etc.
- Linting and formatting: `biome`. Run `./scripts/check-gui-style.sh` before committing.
- Biome config: `biome.json` at repo root.

## Key Scripts
```bash
./scripts/check-source-style.sh             # C++ style
./scripts/check-python-style.sh --fix       # Python lint + format
./scripts/check-gui-style.sh                # Frontend lint + format (biome)
uv pip install -e .                         # rebuild C++ extension
uv run pytest python/tests/                 # Python tests
cd build && ctest                           # C++ tests
uv run scripts/generate_units_md.py         # regenerate docs/UNITS.md
./scripts/test-pypi-wheel.sh                # verify combaero-gui wheel before PyPI publish
```

## Auto-Generated Files
| Target | Source | Regenerate with |
|--------|--------|-----------------|
| `docs/UNITS.md` | `include/units_data.h` | `uv run scripts/generate_units_md.py` |
| `include/thermo_transport_data.h` | `thermo_data_generator/` | see `thermo_data_generator/README.md` |

## Git Workflow
- **Branch protection on `main`:** never push directly to `main`. Always work on a feature/fix branch and open a PR (`feature-branch → main`).
- **Opening PRs:** `git` and `gh` are allowlisted in `.claude/settings.local.json` and run inside the sandbox without approval — including `git push`, `gh pr view/checks/merge`, and (verified) most network operations to `github.com`. `dangerouslyDisableSandbox: true` is only needed when a command writes to a denied path (e.g., `git rebase` that touches `uv.lock`).
- **CI is required to merge:** GitHub Actions checks must pass before a PR can be merged. If main was updated after the PR was opened, re-run checks against the updated base before merging.

## PyPI Release Process
- **Versioning is fully git-tag-driven (`setuptools_scm`) — there is no version string to hand-edit.** Both `pyproject.toml` (combaero) and `gui/pyproject.toml` (combaero-gui) declare `dynamic = ["version"]`; the version comes from the nearest `v*` git tag. Cutting a release is: cut `CHANGELOG.md`'s `[Unreleased]` into a dated `[x.y.z]` block (update the compare-links at the bottom too) → PR → merge to `main` → `git tag -a vx.y.z <commit> && git push origin vx.y.z`.
- **The tag push triggers two independent workflows in parallel:** `publish.yml` (combaero core: wheels for Linux/macOS/Windows via `cibuildwheel` + sdist) and `publish-gui.yml` (combaero-gui: frontend build + wheel + verify smoke test in a clean venv). The core build is slower (3 OSes) — `publish-gui.yml` has to explicitly wait for the matching combaero version to land on PyPI's `/simple/` index (not the JSON API — it updates on a different CDN cache and can falsely report "found" minutes early) before its verify job installs.
- **Whenever a release bumps combaero's minor version, `gui/pyproject.toml`'s `combaero~=x.y` pin must move in lockstep** (edit via shell, see Practical Tips — the Edit tool can't touch it). An unbumped pin resolves a stale combaero from PyPI that may predate classes combaero-gui's backend imports at startup, causing an `ImportError`/`ModuleNotFoundError` the moment the server starts (root cause of the v0.3.1 hotfix and a v0.4.0 near-miss).
- **`combaero-gui`'s publish verify job only installs the built wheel from a clean venv** — it is the one CI path that does *not* run from a repo checkout, so it's the only place that catches a production module importing from an unshipped dev-only tree (see the `validation/` Hard Rule above). `uv run pytest` passing is not sufficient evidence a release is packaging-clean.
- **A pushed git tag should not be moved.** If a release's workflow-only files need a fix after tagging (nothing in the actual package changed), cut the next patch version rather than force-moving/retagging.
- **Yanking a bad PyPI release requires the PyPI web UI** (no API token is available in this environment): project page → the bad version → Options → Yank, with a one-line reason. Yanking doesn't delete anything — pinned installs still work with a warning; unpinned installs skip it.

## Version Sync Checklist
When changing Python version or tool versions, keep in sync:
`.python-version` · `ruff.toml` · `ci.yml` · `validation-tests.yml` · `.pre-commit-config.yaml` · `Makefile`

When changing frontend tool versions, keep in sync:
`gui/frontend/package.json` · `biome.json` (`$schema` version) · `.pre-commit-config.yaml` (biome rev + `@biomejs/biome` pin)
