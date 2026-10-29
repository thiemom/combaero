"""Solved-network results as one pandas table, one row per case (#784).

The layout::

            inputs                     outputs
            inlet   liner              _solve    inlet   liner
            Tt      L     D     ...    success   Tt      Tt      P     ...
      0     650.0   0.30  0.05         True      650.0   1210.4  1.96e5
      1     650.0   0.30  0.05         True      650.0   1488.9  1.95e5

* One row per case. The label is ``0, 1, 2...`` unless a name is given
  (``append_case(..., name=...)``, or the index of the ``params`` frame passed
  to ``NetworkRunner.sweep``).
* Columns are a 3-level MultiIndex ``(kind, entity, quantity)`` with kind
  ``"inputs"`` or ``"outputs"``. ``df["outputs"]`` is the plain results table
  (``df["outputs"]["liner", "Tt"]``); ``df["inputs"]`` is what defined the
  case: boundary conditions, geometry, element and mixture settings, and the
  network-wide physics settings under the entity ``_network``.
* ``outputs`` carries the solve's status under the entity ``_solve``
  (``success``, ``converged``, ``consistent``, ``outcome``, ``final_norm``).
* Entities are addressed by label, falling back to the id. Two entities with
  the same name are refused, naming both ids.
* Initial guesses, the warm start and the solver options are not inputs: they
  do not change a converged answer. They are kept in ``df.attrs``
  (``initial_guesses``, ``warm_start``, ``solver_options``, each keyed by the
  case label), as are the units: ``df.attrs["units"]`` maps a full column
  ``(kind, entity, quantity)`` to its ASCII unit, so quantity names stay
  unit-free.
* The parameters that varied are the input columns whose values differ
  between rows: ``varied_inputs(df)``.

Every case is read from the long, one-row-per-entity frame of
``NetworkResult.to_dataframe()``, which is also the GUI's CSV export. Besides
the unit-labelled output columns (#774) that frame carries the case's inputs,
so an exported CSV reads back into the full table on its own:

* ``input:<key> [<unit>]`` -- an input of the row's entity (``input:L [m]``,
  ``input:composition.source``). Nested fields are dotted.
* ``initial_guess:<key>`` -- the entity's initial guess.
* ``solver:<key>`` -- on the network row only: the solver options and the
  warm start.
* one row of ``type`` ``network`` (id ``_network``) for the network-wide
  settings, and one of ``type`` ``edge`` for an edge with inputs and no
  result row of its own. A flow edge's inputs go on its auto-linked element
  (``__auto_link__<source>__<target>``).

pandas is a combaero-gui dependency only; core combaero never imports this.
"""

import ast
import copy
import csv
import math
import re
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from .schemas import NetworkGraphSchema
from .units import column_unit, input_unit

if TYPE_CHECKING:
    from .runner import NetworkResult

INPUTS = "inputs"
OUTPUTS = "outputs"
#: Entity holding the network-wide physics settings (inputs); also the id of
#: the export's network row.
NETWORK_ENTITY = "_network"
#: Entity holding the solve's status (outputs).
SOLVE_ENTITY = "_solve"
_RESERVED = frozenset({NETWORK_ENTITY, SOLVE_ENTITY})
_LEVELS = ["kind", "entity", "quantity"]

#: Column prefixes of the long frame / CSV export.
INPUT_PREFIX = "input:"
GUESS_PREFIX = "initial_guess:"
SOLVER_PREFIX = "solver:"

_STATUS = ("success", "converged", "consistent", "outcome", "final_norm")
# Columns of the long export that describe a row rather than a quantity.
_LONG_META = frozenset(
    {"type", "kind", "is_boundary", "combustion_method", "id", "label", *_STATUS}
)
# Node/edge data keys that are not inputs: names, and the initial guess
# (recorded in attrs instead).
_NOT_INPUT = frozenset({"id", "label", "initial_guess"})
# Solver settings that do not define the physical case: how it is solved, and
# the report-only incompressible-validity threshold.
_SOLVER_OPTIONS = ("method", "init_strategy", "timeout", "incompressible_mach_limit")
# A unit-labelled export header, "<key> [<unit>]" ("Y[N2] [-]" included).
_HEADER_RE = re.compile(r"^(.*\S) \[([^\[\]]*)\]$")
# The export prefixes a text cell starting with one of these with "'" so a
# spreadsheet does not run it (main._FORMULA_PREFIXES).
_FORMULA_PREFIXES = ("=", "+", "-", "@", "\t", "\r")
_AUTO_LINK = "__auto_link__"

Column = tuple[str, str, str]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_results(source: "str | Path | pd.DataFrame | NetworkResult") -> pd.DataFrame:
    """Read solved-network results into the one-row-per-case table.

    Args:
        source: a CSV exported from the GUI (or written by
            ``NetworkResult.to_dataframe().to_csv``), a ``to_dataframe()``
            frame, a ``NetworkResult``, a table of this shape saved with
            ``df.to_csv()``, a frame already of this shape (returned as a
            copy), or a 0.7 full-detail ``sweep()`` frame (one case per
            ``_sweep_index``, its override columns as the only inputs).

    Returns:
        DataFrame with ``(kind, entity, quantity)`` columns; see the module
        docstring.
    """
    if isinstance(source, pd.DataFrame) and _is_table(source):
        out = source.copy()
        out.attrs = copy.deepcopy(dict(source.attrs))
        return out
    if isinstance(source, (str, Path)) and _csv_is_table(Path(source)):
        return _read_table_csv(Path(source))
    if isinstance(source, pd.DataFrame) and "_sweep_index" in source.columns:
        return _frame(_legacy_sweep_cases(source))
    return _frame([_case_from(source, 0)])


def append_case(
    df: pd.DataFrame,
    source: "NetworkResult | str | Path | pd.DataFrame",
    name: Any = None,
) -> pd.DataFrame:
    """Return ``df`` with one or more cases appended (``df`` is not changed).

    Args:
        df: a table returned by ``load_results``, ``append_case`` or
            ``NetworkRunner.sweep`` (an empty DataFrame starts a new one).
        source: a ``NetworkResult`` (its inputs are read from the network it
            was solved on), an exported CSV, a ``to_dataframe()`` frame, or a
            table of this shape (all its rows are appended).
        name: the new row's label; the next integer by default. Refused when
            ``df`` already has a row of that name. Not allowed when ``source``
            holds several cases.
    """
    if (isinstance(source, pd.DataFrame) and _is_table(source)) or (
        isinstance(source, (str, Path)) and _csv_is_table(Path(source))
    ):
        new = load_results(source)
        if name is not None:
            if len(new) != 1:
                raise ValueError(f"name= names one case; the source holds {len(new)}.")
            new.index = pd.Index([name])
        cases = _cases_of(new)
    else:
        label = _next_label(df) if name is None else name
        cases = [_case_from(source, label)]
    existing = _cases_of(df) if len(df.columns) else []
    taken = {c.name for c in existing}
    for c in cases:
        if c.name in taken:
            raise ValueError(f"The table already has a case named {c.name!r}.")
        taken.add(c.name)
    return _frame(existing + cases)


def varied_inputs(df: pd.DataFrame) -> pd.DataFrame:
    """The ``inputs`` columns whose values differ between rows: the
    parameters that were varied, as an ``(entity, quantity)`` frame."""
    if not len(df.columns) or INPUTS not in df.columns.get_level_values(0):
        return pd.DataFrame(index=df.index)
    inputs = df[INPUTS]
    keep = [col for col in inputs.columns if len({_hashable(v) for v in inputs[col]}) > 1]
    return inputs.loc[:, keep]


# ---------------------------------------------------------------------------
# Inputs in the long frame (written by NetworkResult.to_dataframe)
# ---------------------------------------------------------------------------


def input_cells(
    schema: NetworkGraphSchema,
    row_ids: Iterable[str],
    solver_options: dict[str, Any],
    warm_start: Any = None,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """The input cells of a network for its long frame.

    Args:
        schema: the network the case was solved on (overrides applied).
        row_ids: ids of the rows the frame already has.
        solver_options: how the case was solved.
        warm_start: what the solve started from, if not cold.

    Returns:
        (cells by row id, as ``{column: value}`` with the ``input:``,
        ``initial_guess:`` and ``solver:`` columns; the meta of the rows to
        add: the network row, and an edge with inputs but no row).
    """
    present = set(row_ids)
    cells: dict[str, dict[str, Any]] = {}
    new_rows: list[dict[str, Any]] = []

    def put(row_id: str, data: dict, extra: dict[str, Any]) -> None:
        row = cells.setdefault(row_id, {})
        for key, v in extra.items():
            row[_input_column(key)] = v
        for key, v in _flatten({k: v for k, v in data.items() if k not in _NOT_INPUT}).items():
            if v is not None:
                row[_input_column(key)] = v
        for key, v in _flatten(dict(data.get("initial_guess") or {})).items():
            row[f"{GUESS_PREFIX}{key}"] = v

    settings = schema.solver_settings.model_dump()
    network = {k: v for k, v in settings.items() if k not in _SOLVER_OPTIONS}
    put(NETWORK_ENTITY, network, {})
    for key, v in {**solver_options, "warm_start": warm_start}.items():
        if v is not None:
            cells[NETWORK_ENTITY][f"{SOLVER_PREFIX}{key}"] = v
    new_rows.append(_meta_row("network", "network", NETWORK_ENTITY, ""))

    for n in schema.nodes:
        put(n.id, n.data, {"type": n.type})
        if n.id not in present:
            new_rows.append(_meta_row("node", n.type, n.id, n.data.get("label") or ""))
    for e in schema.edges:
        data = e.data or {}
        if not any(k not in _NOT_INPUT for k in data) and not data.get("initial_guess"):
            continue  # a plain connection defines nothing
        auto = f"{_AUTO_LINK}{e.source}__{e.target}"
        row_id = auto if auto in present and e.id not in present else e.id
        put(row_id, data, {})
        if row_id not in present:
            new_rows.append(
                _meta_row("edge", data.get("type") or "edge", e.id, data.get("label") or "")
            )
    return cells, new_rows


def _input_column(key: str) -> str:
    unit = input_unit(key)
    return f"{INPUT_PREFIX}{key}" if unit is None else f"{INPUT_PREFIX}{key} [{unit}]"


def _meta_row(type_: str, kind: str, row_id: str, label: str) -> dict[str, Any]:
    return {
        "type": type_,
        "kind": kind,
        "is_boundary": False,
        "combustion_method": "",
        "id": row_id,
        "label": label,
    }


def _schema_frame(
    schema: NetworkGraphSchema, solver_options: dict[str, Any], warm_start: Any = None
) -> pd.DataFrame:
    """The long frame of a network that was not solved: its input rows only."""
    cells, rows = input_cells(schema, (), solver_options, warm_start)
    return pd.DataFrame([{**r, **cells.get(r["id"], {})} for r in rows])


# ---------------------------------------------------------------------------
# One case
# ---------------------------------------------------------------------------


class _Case:
    """One row of the table, before it is framed."""

    def __init__(self, name: Any) -> None:
        self.name = name
        self.values: dict[Column, Any] = {}
        self.units: dict[Column, str] = {}
        self.solver_options: dict[str, Any] = {}
        self.initial_guesses: dict[str, dict] = {}
        self.warm_start: Any = None


def _case_from(
    source: "NetworkResult | str | Path | pd.DataFrame",
    name: Any,
    metrics: list[str] | None = None,
) -> _Case:
    from .runner import NetworkResult

    csv_text = False
    if isinstance(source, NetworkResult):
        long_df = source.to_dataframe()
    elif isinstance(source, (str, Path)):
        long_df = pd.read_csv(source, dtype={"id": str, "label": str})
        csv_text = True
    elif isinstance(source, pd.DataFrame):
        long_df = source
    else:
        raise TypeError(
            "Expected a NetworkResult, a results CSV path or a DataFrame, "
            f"got {type(source).__name__}."
        )
    return _build_case(name, long_df, metrics, csv_text=csv_text)


def _build_case(
    name: Any,
    long_df: pd.DataFrame,
    metrics: list[str] | None = None,
    csv_text: bool = False,
) -> _Case:
    """A case from a long (one row per entity) frame."""
    case = _Case(name)
    parsed = _parse_long(long_df, csv_text)
    names = _entity_names(parsed.ids, parsed.labels)

    for eid in parsed.ids:
        entity = NETWORK_ENTITY if eid == NETWORK_ENTITY else names[eid]
        for q, v in parsed.inputs.get(eid, {}).items():
            col = (INPUTS, entity, q)
            case.values[col] = v
            if q in parsed.input_units:
                case.units[col] = parsed.input_units[q]
        if parsed.guesses.get(eid):
            case.initial_guesses[entity] = parsed.guesses[eid]
    case.solver_options = {
        k: v for k, v in parsed.solver.items() if k != "warm_start" and not _missing(v)
    }
    warm = parsed.solver.get("warm_start")
    case.warm_start = None if _missing(warm) else warm

    if parsed.has_status:
        for q in _STATUS:
            case.values[(OUTPUTS, SOLVE_ENTITY, q)] = parsed.status.get(q)
    wanted = None if metrics is None else set(metrics)
    for eid in parsed.ids:
        if eid == NETWORK_ENTITY:
            continue
        entity = names[eid]
        for q, v in parsed.outputs.get(eid, {}).items():
            if wanted is not None and not ({f"{entity}.{q}", f"{eid}.{q}"} & wanted):
                continue
            col = (OUTPUTS, entity, q)
            case.values[col] = v
            if q in parsed.output_units:
                case.units[col] = parsed.output_units[q]
    if wanted is not None:
        # A metric the solve did not produce is a NaN column, as in 0.7.
        present = {
            f"{key}.{q}"
            for eid, qs in parsed.outputs.items()
            if eid in names
            for key in (eid, names[eid])
            for q in qs
        }
        for m in sorted(wanted - present):
            entity, _, q = m.partition(".")
            case.values.setdefault((OUTPUTS, entity, q), math.nan)
    return case


class _Parsed:
    def __init__(self) -> None:
        self.ids: list[str] = []
        self.labels: dict[str, str] = {}
        self.inputs: dict[str, dict[str, Any]] = {}
        self.outputs: dict[str, dict[str, Any]] = {}
        self.guesses: dict[str, dict[str, Any]] = {}
        self.solver: dict[str, Any] = {}
        self.status: dict[str, Any] = {}
        self.has_status = False
        self.input_units: dict[str, str] = {}
        self.output_units: dict[str, str] = {}


def _parse_long(long_df: pd.DataFrame, csv_text: bool) -> _Parsed:
    """The long, unit-labelled frame (``to_dataframe()`` or its CSV) by row."""
    p = _Parsed()
    if not len(long_df):
        return p
    # column -> (section, key)
    layout: dict[Any, tuple[str, str]] = {}
    for col in long_df.columns:
        m = _HEADER_RE.match(str(col))
        key = m.group(1) if m else str(col)
        unit = m.group(2) if m else None
        if key.startswith(INPUT_PREFIX):
            section, key = "input", key[len(INPUT_PREFIX) :]
            if unit is not None:
                p.input_units[key] = unit
        elif key.startswith(GUESS_PREFIX):
            section, key = "guess", key[len(GUESS_PREFIX) :]
        elif key.startswith(SOLVER_PREFIX):
            section, key = "solver", key[len(SOLVER_PREFIX) :]
        elif key in _LONG_META:
            section = "meta"
        else:
            section = "output"
            if unit is not None:
                p.output_units[key] = unit
        layout[col] = (section, key)
    meta = {key: col for col, (section, key) in layout.items() if section == "meta"}
    if "id" not in meta:
        raise ValueError("Not a results export: it has no 'id' column.")
    p.has_status = any(q in meta for q in _STATUS)

    for _, row in long_df.iterrows():
        eid = str(row[meta["id"]])
        label = row[meta["label"]] if "label" in meta else ""
        if eid not in p.labels:
            p.ids.append(eid)
        p.labels[eid] = _unneutralise(label) if isinstance(label, str) else ""
        if p.has_status and not p.status:
            status = {q: _clean(row[meta[q]], csv_text) for q in _STATUS if q in meta}
            # consistent=None ("nothing checked") reads back from a CSV as NaN
            p.status = {q: None if _missing(v) else v for q, v in status.items()}
        for col, (section, key) in layout.items():
            if section == "meta":
                continue
            v = _clean(row[col], csv_text)
            if _missing(v):
                continue
            if section == "solver":
                p.solver[key] = v
            else:
                target = {"input": p.inputs, "guess": p.guesses, "output": p.outputs}[section]
                target.setdefault(eid, {})[key] = v
    return p


def _entity_names(ids: Iterable[str], labels: dict[str, str]) -> dict[str, str]:
    """Entity name of each id: its label, else its id. Refuses two entities
    sharing a name."""
    names: dict[str, str] = {}
    owner: dict[str, str] = {}
    for eid in ids:
        if eid == NETWORK_ENTITY:
            continue
        name = labels.get(eid) or eid
        if name in _RESERVED:
            raise ValueError(
                f"Entity {eid!r} is named {name!r}, which the results table reserves; relabel it."
            )
        other = owner.get(name)
        if other is not None and other != eid:
            raise ValueError(
                f"Two entities are both named {name!r} (ids {other!r} and {eid!r}); "
                "give them distinct labels to tabulate their results."
            )
        owner[name] = eid
        names[eid] = name
    return names


# ---------------------------------------------------------------------------
# Framing
# ---------------------------------------------------------------------------


def _frame(cases: list[_Case]) -> pd.DataFrame:
    """The table of ``cases``: inputs columns first, then outputs, each in
    first-seen order."""
    seen: dict[Column, None] = {}
    for c in cases:
        for col in c.values:
            seen.setdefault(col, None)
    cols = [c for c in seen if c[0] == INPUTS] + [c for c in seen if c[0] == OUTPUTS]
    data = [[c.values.get(col, math.nan) for col in cols] for c in cases]
    df = pd.DataFrame(
        data,
        index=pd.Index([c.name for c in cases]),
        columns=pd.MultiIndex.from_tuples(cols, names=_LEVELS)
        if cols
        else pd.MultiIndex.from_arrays([[], [], []], names=_LEVELS),
        dtype=object,
    ).infer_objects()
    units: dict[Column, str] = {}
    for c in cases:
        units.update(c.units)
    df.attrs = {
        "units": {col: units[col] for col in cols if col in units},
        "solver_options": {c.name: c.solver_options for c in cases},
        "initial_guesses": {c.name: c.initial_guesses for c in cases},
        "warm_start": {c.name: c.warm_start for c in cases},
    }
    return df


def _cases_of(df: pd.DataFrame) -> list[_Case]:
    """The rows of a table back as cases (to re-frame with new ones)."""
    if not len(df.columns):
        return []
    if not _is_table(df):
        raise ValueError(
            "Not a results table: columns must be (kind, entity, quantity) with kind "
            "'inputs' or 'outputs'. Read it with load_results first."
        )
    attrs = df.attrs or {}
    units = attrs.get("units", {})
    cases = []
    for name, row in df.iterrows():
        c = _Case(name)
        for col, v in row.items():
            c.values[tuple(col)] = v
            if tuple(col) in units:
                c.units[tuple(col)] = units[tuple(col)]
        c.solver_options = dict(attrs.get("solver_options", {}).get(name, {}))
        c.initial_guesses = dict(attrs.get("initial_guesses", {}).get(name, {}))
        c.warm_start = attrs.get("warm_start", {}).get(name)
        cases.append(c)
    return cases


def _next_label(df: pd.DataFrame) -> Any:
    if not len(df.index):
        return 0
    ints = [i for i in df.index if isinstance(i, int) and not isinstance(i, bool)]
    if len(ints) == len(df.index):
        return max(ints) + 1
    n = len(df.index)
    while n in df.index:
        n += 1
    return n


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------


def _is_table(df: pd.DataFrame) -> bool:
    cols = df.columns
    return (
        isinstance(cols, pd.MultiIndex)
        and cols.nlevels == 3
        and set(cols.get_level_values(0)) <= {INPUTS, OUTPUTS}
    )


def _csv_is_table(path: Path) -> bool:
    """True for a table saved with ``df.to_csv()`` (three header rows, the
    first naming the kinds), False for a long export."""
    with path.open(newline="") as fh:
        first = next(csv.reader(fh), [])
    return len(first) > 1 and set(first[1:]) <= {INPUTS, OUTPUTS}


def _read_table_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, header=[0, 1, 2], index_col=0)
    df.index.name = None
    df.columns = df.columns.set_names(_LEVELS)
    for col in df.columns:
        if not pd.api.types.is_numeric_dtype(df[col]):
            df[col] = df[col].map(lambda v: _clean(v, True)).infer_objects()
    # attrs do not travel in this CSV: units are re-derived from the quantity
    # names; solver options, initial guesses and the warm start are lost (the
    # long export keeps them).
    units: dict[Column, str] = {}
    for col in df.columns:
        unit = input_unit(col[2]) if col[0] == INPUTS else column_unit(col[2])
        if unit is not None:
            units[tuple(col)] = unit
    df.attrs = {"units": units, "solver_options": {}, "initial_guesses": {}, "warm_start": {}}
    return df


def _legacy_sweep_cases(df: pd.DataFrame) -> list[_Case]:
    """A 0.7 full-detail sweep frame: one case per ``_sweep_index``; the
    columns between it and ``type`` are the override keys."""
    cols = list(df.columns)
    first = cols.index("type") if "type" in cols else len(cols)
    overrides = [c for c in cols[:first] if c != "_sweep_index"]
    cases = []
    for idx, part in df.groupby("_sweep_index", sort=False):
        case = _build_case(idx, part.drop(columns=["_sweep_index", *overrides]))
        row = part.iloc[0]
        for key in overrides:
            entity, _, q = key.partition(".")
            case.values[(INPUTS, entity, q)] = _clean(row[key], False)
        cases.append(case)
    return cases


# ---------------------------------------------------------------------------
# Values
# ---------------------------------------------------------------------------


def _flatten(d: dict, prefix: str = "") -> dict[str, Any]:
    """Nested dicts, and lists of dicts (a wall's layers), to dotted keys."""
    out: dict[str, Any] = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten(v, f"{key}."))
        elif isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
            for i, x in enumerate(v):
                out.update(_flatten(x, f"{key}.{i}."))
        else:
            out[key] = v
    return out


def _missing(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, float):
        return math.isnan(v)
    if isinstance(v, (list, tuple, dict, str)):
        return False
    try:
        return bool(pd.isna(v))
    except (TypeError, ValueError):
        return False


def _clean(v: Any, csv_text: bool) -> Any:
    """A cell as Python: a numpy scalar to its value; from a CSV, list text
    back to a list, True/False text to a bool and the export's formula guard
    undone."""
    if isinstance(v, str):
        if not csv_text:
            return v
        s = v.strip()
        if s.startswith("[") and s.endswith("]"):
            try:
                return ast.literal_eval(s)
            except (ValueError, SyntaxError):
                return v
        if s in ("True", "False"):
            return s == "True"
        return _unneutralise(v)
    if hasattr(v, "item") and not isinstance(v, (list, tuple, dict)):
        try:
            return v.item()
        except (ValueError, AttributeError):
            return v
    return v


def _unneutralise(text: str) -> str:
    if text.startswith("'") and text[1:].startswith(_FORMULA_PREFIXES):
        return text[1:]
    return text


def _hashable(v: Any) -> Any:
    if _missing(v):
        return ("__missing__",)
    if isinstance(v, (list, tuple)):
        return tuple(_hashable(x) for x in v)
    if isinstance(v, dict):
        return tuple(sorted((k, _hashable(x)) for k, x in v.items()))
    return v
