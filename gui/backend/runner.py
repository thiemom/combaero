"""Standalone network runner for programmatic and parametric use.

Loads a saved GUI network JSON, applies boundary-condition overrides
keyed by node label, solves, and returns a NetworkResult that exposes
a full-detail DataFrame, raw state access, and a sweep helper for
parametric studies.

Override key format: ``"<label>.<attribute>"``, where ``<label>`` is the
node's GUI label and ``<attribute>`` is the data field to override
(e.g. ``"air_inlet.m_dot"``, ``"fuel_inlet.Tt"``, ``"outlet.Pt"``).

Typical notebook usage::

    from gui.backend.runner import NetworkRunner

    runner = NetworkRunner.from_file("my_model.json")
    result = runner.solve({"air_inlet.m_dot": 1.2})
    print(result.success, result.final_norm)
    df = result.to_dataframe()

    import pandas as pd
    params = pd.DataFrame({"air_inlet.m_dot": [0.8, 1.0, 1.2]})
    sweep_df = runner.sweep(params, metrics=["channel_1.m_dot [kg/s]"])
"""

import json
from pathlib import Path
from typing import Any

import combaero as cb
from combaero.network import NetworkSolver
import pandas as pd

from .graph_builder import build_network_from_schema
from .schemas import ElementResult, NetworkGraphSchema, NodeResult, StateResult
from .units import label_with_unit


# ---------------------------------------------------------------------------
# Shared result-building helper (used by both the HTTP API and NetworkRunner)
# ---------------------------------------------------------------------------


def _schema_maps(
    schema: NetworkGraphSchema,
) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """Build (id_to_label, id_to_kind, id_to_method) dicts from a schema.

    Shared by NetworkRunner.__init__ and the HTTP export endpoint so label
    resolution is consistent in both paths.
    """
    id_to_label: dict[str, str] = {}
    id_to_kind: dict[str, str] = {}
    id_to_method: dict[str, str] = {}
    for n in schema.nodes:
        lbl = n.data.get("label") or ""
        id_to_label[n.id] = lbl
        id_to_kind[n.id] = n.type
        if n.type == "combustor":
            id_to_method[n.id] = n.data.get("method", "complete")
    for e in schema.edges:
        lbl = (e.data or {}).get("label") or ""
        id_to_label[e.id] = lbl
        id_to_kind[e.id] = (e.data or {}).get("type") or "flow"
        auto_id = f"__auto_link__{e.source}__{e.target}"
        id_to_label[auto_id] = lbl
        id_to_kind[auto_id] = "lossless_connection"
    return id_to_label, id_to_kind, id_to_method


def _build_result_objects(
    result: dict,
    net,
) -> tuple[dict[str, NodeResult], dict[str, ElementResult], dict[str, dict]]:
    """Extract typed result objects from the raw solver output dict.

    Factored out of the HTTP solve handler so that both the FastAPI layer
    and the standalone NetworkRunner share identical result construction.
    """
    success = bool(result.get("__success__", False))

    node_results: dict[str, NodeResult] = {}
    for node_id in net.nodes:
        prefix = f"{node_id}."
        node_vals = {k[len(prefix):]: v for k, v in result.items() if k.startswith(prefix)}

        y_vals: list[float] = []
        x_vals: list[float] = []
        i = 0
        while f"Y[{i}]" in node_vals:
            y_vals.append(float(node_vals.pop(f"Y[{i}]")))
            if f"X[{i}]" in node_vals:
                x_vals.append(float(node_vals.pop(f"X[{i}]")))
            i += 1

        state_res = StateResult(
            T=float(node_vals.pop("T", 300.0)),
            P=float(node_vals.pop("P", 101325.0)),
            Pt=node_vals.pop("Pt", None),
            Tt=node_vals.pop("Tt", None),
            Y=y_vals,
            X=x_vals or None,
            **node_vals,
        )
        node_results[node_id] = NodeResult(state=state_res, success=success)

    element_results: dict[str, ElementResult] = {}
    elem_diags = result.get("__element_diag__", {})
    for elem_id in net.elements:
        diag = elem_diags.get(elem_id, {}).copy()
        m_dot = float(result.get(f"{elem_id}.m_dot", diag.get("m_dot_com", 0.0)))
        diag.pop("m_dot", None)
        element_results[elem_id] = ElementResult(m_dot=m_dot, success=success, **diag)

    def _safe_float(v: Any) -> Any:
        try:
            if isinstance(v, (list, tuple, dict)):
                return v
            return float(v)
        except (TypeError, ValueError):
            return v

    edge_results: dict[str, dict] = {}
    for edge_id in net.walls:
        edge_keys = {
            k.split(".", 1)[1]: _safe_float(v)
            for k, v in result.items()
            if k.startswith(f"{edge_id}.")
        }
        if edge_keys:
            edge_results[edge_id] = edge_keys

    return node_results, element_results, edge_results


# ---------------------------------------------------------------------------
# NetworkResult
# ---------------------------------------------------------------------------


class NetworkResult:
    """Result of a single network solve.

    Carries the full solved state for programmatic access, DataFrame export,
    and use as a warm start for follow-on operations such as boundary swaps.

    Attributes:
        success: True when the solver converged within the residual tolerance.
        message: Human-readable solver status string.
        final_norm: Euclidean norm of the residual vector at convergence.
    """

    def __init__(
        self,
        raw: dict,
        node_results: dict[str, NodeResult],
        element_results: dict[str, ElementResult],
        edge_results: dict[str, dict],
        net: Any,
        schema: NetworkGraphSchema,
        id_to_label: dict[str, str],
        id_to_kind: dict[str, str],
        id_to_method: dict[str, str],
    ) -> None:
        self.success: bool = bool(raw.get("__success__", False))
        self.message: str = str(raw.get("__message__", ""))
        self.final_norm: float | None = raw.get("__final_norm__")
        self._raw = raw
        self._node_results = node_results
        self._element_results = element_results
        self._edge_results = edge_results
        self._net = net
        self._schema = schema
        self._id_to_label = id_to_label
        self._id_to_kind = id_to_kind
        self._id_to_method = id_to_method

    def get(self, key: str) -> float:
        """Get a scalar result by raw solver key (e.g. ``'combustor.T'``).

        Keys follow the pattern ``'<element_or_node_id>.<quantity>'``.
        """
        if key not in self._raw:
            raise KeyError(
                f"Result key '{key}' not found. "
                "Keys follow the format '<id>.<quantity>'."
            )
        return float(self._raw[key])

    def node_state(self, label: str) -> dict:
        """Return the full solved thermodynamic state for a node by label or id.

        The returned dict contains T, P, Pt, Tt, rho, mach, Y, X and any
        additional quantities exposed by that node type.
        """
        node_id = self._resolve_label(label)
        res = self._node_results.get(node_id)
        if res is None:
            raise KeyError(f"No node with label or id '{label}' in results.")
        return res.state.model_dump()

    def _resolve_label(self, label: str) -> str:
        """Return the node/element id for a given label, or the label itself if it is an id."""
        for nid, lbl in self._id_to_label.items():
            if lbl == label:
                return nid
        if label in self._node_results or label in self._element_results:
            return label
        raise KeyError(f"No node or element with label or id '{label}'.")

    def to_dataframe(self) -> pd.DataFrame:
        """Full-detail DataFrame equivalent to the GUI CSV export.

        Returns one row per network entity (nodes, flow elements, thermal walls)
        with all available diagnostic quantities.  Column names carry unit
        annotations matching the GUI export (e.g. ``'P [Pa]'``, ``'T [K]'``).
        """
        from combaero.network.components import TeeJunctionElement as _TeeJE

        _BOUNDARY_KINDS = {"mass_boundary", "pressure_boundary"}
        species_labels: list[str] = list(cb.species.names)
        data: list[dict] = []

        for node_id, res in self._node_results.items():
            state_data = res.state.model_dump()
            Y = state_data.pop("Y", [])
            X = state_data.pop("X", None)
            kind = self._id_to_kind.get(node_id, "")
            row: dict = {
                "type": "node",
                "kind": kind,
                "is_boundary": kind in _BOUNDARY_KINDS,
                "combustion_method": self._id_to_method.get(node_id, ""),
                "id": node_id,
                "label": self._id_to_label.get(node_id, ""),
                **state_data,
            }
            for i, y_val in enumerate(Y):
                name = species_labels[i] if i < len(species_labels) else str(i)
                row[f"Y[{name}]"] = y_val
            if X:
                for i, x_val in enumerate(X):
                    name = species_labels[i] if i < len(species_labels) else str(i)
                    row[f"X[{name}]"] = x_val
            data.append(row)

        for elem_id, res in self._element_results.items():
            elem_data = res.model_dump()
            elem_obj = self._net.elements.get(elem_id)
            if isinstance(elem_obj, _TeeJE):
                common_res = self._node_results.get(elem_obj.common_node)
                if common_res:
                    common_state = common_res.state.model_dump()
                    Y = common_state.pop("Y", [])
                    X = common_state.pop("X", None) or []
                    for k, v in common_state.items():
                        if v is not None:
                            elem_data[k] = v
                    for i, y_val in enumerate(Y):
                        name = species_labels[i] if i < len(species_labels) else str(i)
                        elem_data[f"Y[{name}]"] = y_val
                    for i, x_val in enumerate(X):
                        name = species_labels[i] if i < len(species_labels) else str(i)
                        elem_data[f"X[{name}]"] = x_val
            data.append(
                {
                    "type": "element",
                    "kind": self._id_to_kind.get(elem_id, ""),
                    "is_boundary": False,
                    "combustion_method": "",
                    "id": elem_id,
                    "label": self._id_to_label.get(elem_id, ""),
                    **elem_data,
                }
            )

        for wall_id, res in self._edge_results.items():
            data.append(
                {
                    "type": "thermal_wall",
                    "kind": "thermal_wall",
                    "is_boundary": False,
                    "combustion_method": "",
                    "id": wall_id,
                    "label": self._id_to_label.get(wall_id, ""),
                    **res,
                }
            )

        df = pd.DataFrame(data)
        return df.rename(columns={col: label_with_unit(col) for col in df.columns})


# ---------------------------------------------------------------------------
# NetworkRunner
# ---------------------------------------------------------------------------


class NetworkRunner:
    """Programmatic driver for a saved combaero-gui network.

    Load a network JSON file saved from the GUI, apply boundary-condition
    overrides by node label, and run parametric sweeps without starting
    the web server.

    Override key format: ``"<label>.<attribute>"`` where ``<label>`` is the
    node's GUI label and ``<attribute>`` is the data field name
    (e.g. ``"air_inlet.m_dot"``, ``"fuel_inlet.Tt"``, ``"outlet.Pt"``).
    The label falls back to the node id when no label is set.

    Example::

        runner = NetworkRunner.from_file("network.json")

        # Single solve with BC overrides
        result = runner.solve({"air_inlet.m_dot": 1.2, "fuel_inlet.Tt": 850.0})
        print(result.success, result.final_norm)

        # Parametric sweep
        import pandas as pd
        params = pd.DataFrame({
            "air_inlet.m_dot": [0.8, 1.0, 1.2],
            "fuel_inlet.m_dot": [0.03, 0.04, 0.05],
        })
        sweep_df = runner.sweep(params, metrics=["combustor.T", "hot_channel.m_dot"])
    """

    def __init__(self, schema: NetworkGraphSchema) -> None:
        self._schema = schema
        self._id_to_label, self._id_to_kind, self._id_to_method = _schema_maps(schema)
        self._label_to_node_id: dict[str, str] = {
            lbl: nid for nid, lbl in self._id_to_label.items() if lbl
        }

    @classmethod
    def from_file(cls, path: str | Path) -> "NetworkRunner":
        """Load a GUI-saved network JSON file and return a runner."""
        with Path(path).open() as fh:
            return cls.from_dict(json.load(fh))

    @classmethod
    def from_dict(cls, d: dict) -> "NetworkRunner":
        """Construct a runner from an already-loaded dict."""
        return cls(NetworkGraphSchema.model_validate(d))

    # ------------------------------------------------------------------
    # Solve
    # ------------------------------------------------------------------

    def solve(
        self,
        overrides: dict[str, Any] | None = None,
        method: str = "hybr",
        init_strategy: str = "default",
        timeout: float | None = 180.0,
    ) -> NetworkResult:
        """Solve the network, optionally patching boundary conditions.

        Each call is stateless: a deep copy of the schema is made before
        applying overrides, so repeated calls (e.g. inside an optimisation
        loop) do not interfere with each other.

        Args:
            overrides: Dict of ``"<label>.<attribute>": value`` pairs.
                Examples: ``{"air_inlet.m_dot": 1.2}``,
                ``{"outlet.Pt": 200000.0}``.
            method: scipy root-finding method passed to NetworkSolver.
            init_strategy: Solver initialisation strategy
                (``"default"``, ``"incompressible_warmstart"``,
                ``"homotopy"``).
            timeout: Per-solve soft timeout in seconds.

        Returns:
            NetworkResult with solved state and DataFrame export.
        """
        schema = self._schema.model_copy(deep=True)
        if overrides:
            self._apply_overrides(schema, overrides)

        net = build_network_from_schema(schema)
        solver = NetworkSolver(net)
        raw = solver.solve(method=method, init_strategy=init_strategy, timeout=timeout)
        node_results, element_results, edge_results = _build_result_objects(raw, net)

        return NetworkResult(
            raw=raw,
            node_results=node_results,
            element_results=element_results,
            edge_results=edge_results,
            net=net,
            schema=schema,
            id_to_label=dict(self._id_to_label),
            id_to_kind=dict(self._id_to_kind),
            id_to_method=dict(self._id_to_method),
        )

    # ------------------------------------------------------------------
    # Sweep
    # ------------------------------------------------------------------

    def sweep(
        self,
        params: pd.DataFrame,
        metrics: list[str] | None = None,
        method: str = "hybr",
        init_strategy: str = "default",
        timeout: float | None = 180.0,
    ) -> pd.DataFrame:
        """Run a parametric sweep over rows of ``params``.

        Each row of ``params`` is a set of BC overrides (column names are
        override keys in ``<label>.<attribute>`` format).

        Args:
            params: DataFrame whose columns are override keys and whose rows
                are individual parameter combinations.
            metrics: Optional list of raw solver keys to extract per solve
                (e.g. ``["combustor.T", "hot_channel.m_dot"]``).

                * When supplied, returns a compact DataFrame with one row per
                  solve: the param columns, ``success``, ``final_norm``, and
                  the requested metric values.
                * When ``None``, returns the full ``to_dataframe()`` output
                  for each solve, stacked with the param values and a
                  ``_sweep_index`` column prepended.

            method: scipy root-finding method forwarded to the solver.
            init_strategy: Solver initialisation strategy.
            timeout: Per-solve timeout in seconds.

        Returns:
            DataFrame with one or more rows per parameter combination.
        """
        rows: list[Any] = []
        for idx, param_row in params.iterrows():
            overrides = {col: param_row[col] for col in params.columns}
            result = self.solve(
                overrides=overrides,
                method=method,
                init_strategy=init_strategy,
                timeout=timeout,
            )
            if metrics is not None:
                row: dict[str, Any] = dict(param_row)
                row["success"] = result.success
                row["final_norm"] = result.final_norm
                for m in metrics:
                    try:
                        row[m] = result.get(m)
                    except (KeyError, ValueError):
                        row[m] = None
                rows.append(row)
            else:
                detail = result.to_dataframe()
                for col in reversed(list(params.columns)):
                    detail.insert(0, col, param_row[col])
                detail.insert(0, "_sweep_index", idx)
                detail["success"] = result.success
                detail["final_norm"] = result.final_norm
                rows.append(detail)

        if not rows:
            return pd.DataFrame()
        if metrics is not None:
            return pd.DataFrame(rows)
        return pd.concat(rows, ignore_index=True)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _apply_overrides(
        self,
        schema: NetworkGraphSchema,
        overrides: dict[str, Any],
    ) -> None:
        """Patch schema node data in-place with BC override values."""
        for key, value in overrides.items():
            if "." not in key:
                raise ValueError(
                    f"Override key '{key}' must be '<label>.<attribute>' "
                    "(e.g. 'air_inlet.m_dot')."
                )
            label, attr = key.split(".", 1)
            node_id = self._label_to_node_id.get(label)
            if node_id is None:
                if any(n.id == label for n in schema.nodes):
                    node_id = label
                else:
                    raise KeyError(
                        f"No node with label or id '{label}' found in this network."
                    )
            for node in schema.nodes:
                if node.id == node_id:
                    node.data[attr] = value
                    break
