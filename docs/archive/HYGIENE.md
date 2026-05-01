# Documentation Hygiene & Classification Audit

This document classifies the existing CombAero documentation to distinguish active developer references from historical feature reports and technical debt.

## 1. Active Documentation (High Signal)

These files are the primary sources of truth for developers and agents. They should be maintained in the root `docs/` folder.

| File | Role | Status | Action Needed |
| :--- | :--- | :--- | :--- |
| [API_CPP.md](API_CPP.md) | C++ Interface | **Active** | None |
| [API_PYTHON.md](API_PYTHON.md) | Python Interface | **Active** | None |
| [BUILDING.md](BUILDING.md) | Setup Guide | **Active** | None |
| [DEVELOPMENT.md](DEVELOPMENT.md) | Workflow Guide | **Active** | None |
| [GUI_TECHNICAL.md](GUI_TECHNICAL.md) | Machine-Readable Schemas | **Active** | None |
| [UNITS.md](UNITS.md) | SI Reference | **Active** | Regenerate via script periodically |
| [METHODS_REFERENCE.md](METHODS_REFERENCE.md) | Performance & Design | **Active** | Update benchmarks for Python 3.13+ |

## 2. Roadmap & Technical Debt

These files track the evolution of the project and identified gaps.

| File | Role | Status | Action Needed |
| :--- | :--- | :--- | :--- |
| [NETWORK_ROADMAP.md](NETWORK_ROADMAP.md) | Strategic Vision | **Active** | None (Phases 1–4 archived) |
| [REVIEW_NETWORK_GUI.md](REVIEW_NETWORK_GUI.md) | Technical Debt Audit | **Active** | Tick off items as they are resolved |
| [JACOBIAN_DIFFERENCE_REPORTING.md](JACOBIAN_DIFFERENCE_REPORTING.md) | Dev Tool Guide | **Active / WIP** | Add missing `report_` calls to tests |
| [PACKAGING.md](PACKAGING.md) | Architecture Rationale | **Active** | None (Refactored to Rationale) |

## 3. Historical Feature Reports (Archive)

These documents describe the implementation of specific features. While valuable for context, they clutter the main documentation set.

| File | Role | Status | Recommendation |
| :--- | :--- | :--- | :--- |
| `docs/DEVELOPMENT.md` | **Active** | Core development & environment setup |
| `docs/METHODS_REFERENCE.md` | **Active** | Numerical performance & profiling |
| `docs/JACOBIAN_DIFFERENCE_REPORTING.md` | **Active** | Verification workflow for derivatives |
| `docs/REVIEW_NETWORK_GUI.md` | **Active** | Next-phase design roadmap |
| `gui/frontend/README.md` | **Active** | React/Vite architecture guide |
| `docs/archive/MACH_LIMIT_BENCHMARK_RESULT.md` | **Archive** | Mach limit benchmark snapshot |
| `docs/archive/MACH_LIMIT_SOLVER_LIMITS.md` | **Archive** | Solver convergence limit snapshot |
| `docs/archive/COMPRESSIBLE_FLOW_IMPLEMENTATION.md` | **Archive** | Original design report |

---

## 4. Pending Developer Tasks (Extracted from Docs)

Based on the audit, the following technical tasks are documented but not yet fully implemented:

1. **Jacobian Coverage**: [DONE] Integrated `report_jacobian_difference()` into compressible flow tests in `tests/test_solver_jacobians.cpp`.
2. **Inverse Solvers**: [DONE] Bound `calc_T_from_s_mass`, `calc_T_from_u_mass`, and flash functions in PyBind11 (merged in #108).
3. **Packaging**: [DONE] Populated missing metadata in `pyproject.toml` (authors, keywords, classifiers, URLs).
4. **GUI Refinement**: [DONE] Exposed element-specific friction and Nusselt multipliers across Combustor, Momentum Chamber, and Discrete Loss elements with standardized grid UI.
