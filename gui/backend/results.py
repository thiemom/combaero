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
from .units import column_unit

if TYPE_CHECKING:
    from .runner import NetworkResult

INPUTS = "inputs"
OUTPUTS = "outputs"
#: Entity holding the network-wide physics settings (inputs).
NETWORK_ENTITY = "_network"
#: Entity holding the solve's status (outputs).
SOLVE_ENTITY = "_solve"
_RESERVED = frozenset({NETWORK_ENTITY, SOLVE_ENTITY})
_LEVELS = ["kind", "entity", "quantity"]

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


def load_results(
    source: "str | Path | pd.DataFrame | NetworkResult",
    network: "str | Path | dict | NetworkGraphSchema | None" = None,
) -> pd.DataFrame:
    """Read solved-network results into the one-row-per-case table.

    Args:
        source: a CSV exported from the GUI (or written by
            ``NetworkResult.to_dataframe().to_csv``), a table of this shape
            saved with ``df.to_csv()``, a ``to_dataframe()`` frame, a frame
            already of this shape (returned as a copy), a 0.7 full-detail
            ``sweep()`` frame (one case per ``_sweep_index``, its override
            columns as inputs), or a ``NetworkResult``.
        network: the network file (path, dict or schema) a CSV or a
            ``to_dataframe()`` frame was solved from, to fill its inputs. An
            export carries outputs only; without ``network`` the case has no
            inputs columns. Ignored for a ``NetworkResult``, which carries its
            own network.

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
        cases = _legacy_sweep_cases(source, network)
        return _frame(cases)
    return _frame([_case_from(source, 0, network)])


def append_case(
    df: pd.DataFrame,
    source: "NetworkResult | str | Path | pd.DataFrame",
    name: Any = None,
    network: "str | Path | dict | NetworkGraphSchema | None" = None,
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
        network: as for ``load_results``, for a CSV or ``to_dataframe()``
            source.
    """
    if (
        isinstance(source, pd.DataFrame)
        and _is_table(source)
        or (isinstance(source, (str, Path)) and _csv_is_table(Path(source)))
    ):
        new = load_results(source)
        if name is not None:
            if len(new) != 1:
                raise ValueError(f"name= names one case; the source holds {len(new)}.")
            new.index = pd.Index([name])
        cases = _cases_of(new)
    else:
        label = _next_label(df) if name is None else name
        cases = [_case_from(source, label, network)]
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
    network: "str | Path | dict | NetworkGraphSchema | None",
    metrics: list[str] | None = None,
) -> _Case:
    from .runner import NetworkResult

    if isinstance(source, NetworkResult):
        long_df = source.to_dataframe()
        schema = source._schema
        case = _build_case(name, long_df, schema, metrics)
        case.solver_options = dict(getattr(source, "_solver_options", {}))
        case.warm_start = getattr(source, "_warm_start", None)
        return case
    if isinstance(source, (str, Path)):
        long_df = _read_long_csv(Path(source))
    elif isinstance(source, pd.DataFrame):
        long_df = source
    else:
        raise TypeError(
            "Expected a NetworkResult, a results CSV path or a DataFrame, "
            f"got {type(source).__name__}."
        )
    return _build_case(name, long_df, _schema_of(network), metrics)


def _build_case(
    name: Any,
    long_df: pd.DataFrame | None,
    schema: NetworkGraphSchema | None,
    metrics: list[str] | None = None,
) -> _Case:
    """A case from a long (one row per entity) frame and/or the network."""
    case = _Case(name)
    inputs: dict[str, dict[str, Any]] = {}
    labels: dict[str, str] = {}
    aliases: dict[str, str] = {}
    if schema is not None:
        inputs, labels, aliases, guesses, options = _schema_inputs(schema)
        case.initial_guesses = guesses
        case.solver_options = options

    outputs: dict[str, dict[str, Any]] = {}
    status: dict[str, Any] = {}
    units_by_key: dict[str, str] = {}
    if long_df is not None and len(long_df):
        outputs, out_labels, status, units_by_key = _long_outputs(long_df)
        for eid, lbl in out_labels.items():
            labels.setdefault(eid, lbl)

    names = _entity_names([*inputs, *outputs], labels, aliases)

    for eid, quantities in inputs.items():
        entity = NETWORK_ENTITY if eid == NETWORK_ENTITY else names[eid]
        for q, v in quantities.items():
            col = (INPUTS, entity, q)
            case.values[col] = v
            unit = column_unit(q)
            if unit is not None:
                case.units[col] = unit

    if long_df is not None and len(long_df):
        for q in _STATUS:
            case.values[(OUTPUTS, SOLVE_ENTITY, q)] = status.get(q)
    wanted = None if metrics is None else set(metrics)
    for eid, quantities in outputs.items():
        entity = names[eid]
        for q, v in quantities.items():
            if wanted is not None and not ({f"{entity}.{q}", f"{eid}.{q}"} & wanted):
                continue
            col = (OUTPUTS, entity, q)
            case.values[col] = v
            if q in units_by_key:
                case.units[col] = units_by_key[q]
    if wanted is not None:
        # A metric the solve did not produce is a NaN column, as in 0.7.
        present = {
            f"{key}.{q}" for eid, qs in outputs.items() for key in (eid, names[eid]) for q in qs
        }
        for m in sorted(wanted - present):
            entity, _, q = m.partition(".")
            case.values.setdefault((OUTPUTS, entity, q), math.nan)
    return case


def _schema_inputs(
    schema: NetworkGraphSchema,
) -> tuple[
    dict[str, dict[str, Any]], dict[str, str], dict[str, str], dict[str, dict], dict[str, Any]
]:
    """(inputs by id, labels, auto-link aliases, initial guesses, solver
    options) of a network. The inputs are the values the network file gives:
    a field left at its default is absent."""
    inputs: dict[str, dict[str, Any]] = {}
    labels: dict[str, str] = {}
    aliases: dict[str, str] = {}
    guesses: dict[str, dict] = {}
    for n in schema.nodes:
        labels[n.id] = n.data.get("label") or ""
        flat = _flatten({k: v for k, v in n.data.items() if k not in _NOT_INPUT})
        inputs[n.id] = {"type": n.type, **flat}
        if n.data.get("initial_guess"):
            guesses[labels[n.id] or n.id] = dict(n.data["initial_guess"])
    for e in schema.edges:
        data = e.data or {}
        labels[e.id] = data.get("label") or ""
        # A flow edge is solved as an auto-linked element (runner._schema_maps).
        aliases[f"{_AUTO_LINK}{e.source}__{e.target}"] = e.id
        flat = _flatten({k: v for k, v in data.items() if k not in _NOT_INPUT})
        if flat:  # a plain connection defines nothing
            inputs[e.id] = flat
        if data.get("initial_guess"):
            guesses[labels[e.id] or e.id] = dict(data["initial_guess"])
    settings = schema.solver_settings.model_dump()
    network = {k: v for k, v in settings.items() if k not in _SOLVER_OPTIONS and v is not None}
    inputs = {NETWORK_ENTITY: network, **inputs}
    options = {k: settings.get(k) for k in _SOLVER_OPTIONS}
    return inputs, labels, aliases, guesses, options


def _long_outputs(
    long_df: pd.DataFrame,
) -> tuple[dict[str, dict[str, Any]], dict[str, str], dict[str, Any], dict[str, str]]:
    """(quantities by id, labels by id, solve status, unit by quantity) of a
    long, unit-labelled results frame (``to_dataframe()`` or its CSV)."""
    bare: dict[str, str] = {}
    units: dict[str, str] = {}
    for col in long_df.columns:
        m = _HEADER_RE.match(str(col))
        key = m.group(1) if m else str(col)
        bare[col] = key
        if m:
            units[key] = m.group(2)
    if "id" not in bare.values():
        raise ValueError("Not a results export: it has no 'id' column.")
    by_key = {key: col for col, key in bare.items()}
    outputs: dict[str, dict[str, Any]] = {}
    labels: dict[str, str] = {}
    status: dict[str, Any] = {}
    for _, row in long_df.iterrows():
        eid = str(row[by_key["id"]])
        label = row[by_key["label"]] if "label" in by_key else ""
        labels[eid] = _unneutralise(label) if isinstance(label, str) else ""
        if not status:
            # consistent=None ("nothing checked") reads back from a CSV as NaN
            status = {q: _clean(row[by_key[q]]) for q in _STATUS if q in by_key}
            status = {q: None if _missing(v) else v for q, v in status.items()}
        quantities: dict[str, Any] = {}
        for col, key in bare.items():
            if key in _LONG_META:
                continue
            v = _clean(row[col])
            if _missing(v):
                continue
            quantities[key] = v
        outputs[eid] = quantities
    return outputs, labels, status, units


def _entity_names(
    ids: Iterable[str], labels: dict[str, str], aliases: dict[str, str]
) -> dict[str, str]:
    """Entity name of each id: its label, else its id. Refuses two entities
    sharing a name (an auto-linked edge element is its edge)."""
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
        canonical = aliases.get(eid, eid)
        other = owner.get(name)
        if other is not None and other != canonical:
            raise ValueError(
                f"Two entities are both named {name!r} (ids {other!r} and {canonical!r}); "
                "give them distinct labels to tabulate their results."
            )
        owner[name] = canonical
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
            df[col] = df[col].map(_clean).infer_objects()
    # Units do not travel in a CSV: re-derive them from the quantity names.
    df.attrs = {
        "units": {
            tuple(col): unit for col in df.columns if (unit := column_unit(col[2])) is not None
        },
        "solver_options": {},
        "initial_guesses": {},
        "warm_start": {},
    }
    return df


def _read_long_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"id": str, "label": str})


def _legacy_sweep_cases(
    df: pd.DataFrame, network: "str | Path | dict | NetworkGraphSchema | None"
) -> list[_Case]:
    """A 0.7 full-detail sweep frame: one case per ``_sweep_index``; the
    columns between it and ``type`` are the override keys."""
    cols = list(df.columns)
    first = cols.index("type") if "type" in cols else len(cols)
    overrides = [c for c in cols[:first] if c != "_sweep_index"]
    schema = _schema_of(network)
    cases = []
    for idx, part in df.groupby("_sweep_index", sort=False):
        body = part.drop(columns=["_sweep_index", *overrides])
        case = _build_case(idx, body, schema)
        row = part.iloc[0]
        for key in overrides:
            entity, _, q = key.partition(".")
            case.values[(INPUTS, entity, q)] = _clean(row[key])
        cases.append(case)
    return cases


def _schema_of(
    network: "str | Path | dict | NetworkGraphSchema | None",
) -> NetworkGraphSchema | None:
    if network is None or isinstance(network, NetworkGraphSchema):
        return network
    if isinstance(network, (str, Path)):
        import json

        with Path(network).open() as fh:
            network = json.load(fh)
    return NetworkGraphSchema.model_validate(network)


# ---------------------------------------------------------------------------
# Values
# ---------------------------------------------------------------------------


def _flatten(d: dict, prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten(v, f"{key}."))
        else:
            out[key] = v
    return out


def _missing(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, float):
        return math.isnan(v)
    try:
        return bool(pd.isna(v)) if not isinstance(v, (list, tuple, dict, str)) else False
    except (TypeError, ValueError):
        return False


def _clean(v: Any) -> Any:
    """A cell as read from a CSV: list text back to a list, True/False text
    to a bool, NaN to None for the status."""
    if isinstance(v, str):
        s = v.strip()
        if s.startswith("[") and s.endswith("]"):
            try:
                return ast.literal_eval(s)
            except (ValueError, SyntaxError):
                return v
        if s in ("True", "False"):
            return s == "True"
    if hasattr(v, "item") and not isinstance(v, (list, tuple, dict)):
        try:
            return v.item()
        except (ValueError, AttributeError):
            return v
    return v


def _unneutralise(label: str) -> str:
    if label.startswith("'") and label[1:].startswith(_FORMULA_PREFIXES):
        return label[1:]
    return label


def _hashable(v: Any) -> Any:
    if _missing(v):
        return ("__missing__",)
    if isinstance(v, (list, tuple)):
        return tuple(_hashable(x) for x in v)
    if isinstance(v, dict):
        return tuple(sorted((k, _hashable(x)) for k, x in v.items()))
    return v
