/**
 * Heat source (#792): an EnergyBoundary on the plenum, momentum chamber or
 * combustor it is wired to by ONE edge. The backend (graph_builder) refuses
 * any other wiring by name; these helpers keep the canvas from drawing it.
 */

import type { Connection, Edge, Node } from "reactflow";

export const HEAT_SOURCE_KIND = "heat_source";
/** The heat source's one handle, a connection source. */
export const HEAT_SOURCE_HANDLE = "heat-out";
/** Node kinds a heat source can heat: an energy balance sets their Tt. */
export const HEAT_SOURCE_TARGETS = [
	"plenum",
	"momentum_chamber",
	"combustor",
] as const;

export type HeatSourceMode = "Q" | "fraction" | "Tt_target";

export const HEAT_SOURCE_MODES: { value: HeatSourceMode; label: string }[] = [
	{ value: "Q", label: "Heat flow Q" },
	{ value: "fraction", label: "Fraction of sensible enthalpy" },
	{ value: "Tt_target", label: "Outlet Tt (total)" },
];

/** A dropped heat source: a working example, 1 kW into the node. */
export const HEAT_SOURCE_DEFAULTS = {
	mode: "Q" as HeatSourceMode,
	Q: 1000.0,
	fraction: 0.01,
	Tt_target: 600.0,
};

/**
 * Why a connection touching a heat source cannot be drawn, or null when it
 * can: out of a heat source, into one node it may heat, and only one edge.
 */
export function heatSourceConnectionError(
	c: Connection,
	nodes: Node[],
	edges: Edge[],
): string | null {
	const byId = new Map(nodes.map((n) => [n.id, n]));
	const source = c.source ? byId.get(c.source) : undefined;
	const target = c.target ? byId.get(c.target) : undefined;
	if (target?.type === HEAT_SOURCE_KIND) {
		return "A heat source takes no inflow: draw its edge from it to the node it heats.";
	}
	if (source?.type !== HEAT_SOURCE_KIND) return null;
	if (
		!target ||
		!(HEAT_SOURCE_TARGETS as readonly string[]).includes(target.type ?? "")
	) {
		return "A heat source heats a plenum, momentum chamber or combustor.";
	}
	if (edges.some((e) => e.source === source.id || e.target === source.id)) {
		return "A heat source is wired to one node only.";
	}
	return null;
}

/** The edge drawn from a heat source: a dashed red link, no flow. */
export function heatSourceEdgeProps(): Pick<Edge, "type" | "style" | "data"> {
	return {
		type: "default",
		style: { stroke: "#e11d48", strokeDasharray: "4 3" },
		data: { type: "heat" },
	};
}

/** The set value in one line, for the node card. */
export function heatSourceSummary(data: Record<string, any>): string {
	const mode: HeatSourceMode = data.mode ?? HEAT_SOURCE_DEFAULTS.mode;
	if (mode === "Tt_target") {
		const T = data.Tt_target ?? HEAT_SOURCE_DEFAULTS.Tt_target;
		return `Tt = ${Number(T).toFixed(1)} K (total)`;
	}
	if (mode === "fraction") {
		const f = data.fraction ?? HEAT_SOURCE_DEFAULTS.fraction;
		return `f = ${(Number(f) * 100).toFixed(2)} %`;
	}
	return `Q = ${formatHeat(data.Q ?? HEAT_SOURCE_DEFAULTS.Q)}`;
}

/** A heat flow [W] in W or kW; NaN (undetermined) as "--". */
export function formatHeat(Q: number | null | undefined): string {
	if (Q === null || Q === undefined || !Number.isFinite(Q)) return "--";
	return Math.abs(Q) >= 1e4
		? `${(Q / 1e3).toFixed(2)} kW`
		: `${Q.toFixed(1)} W`;
}

/** The node a heat source is wired to, if any. */
export function heatedNodeId(id: string, edges: Edge[]): string | undefined {
	const e = edges.find((x) => x.source === id || x.target === id);
	return e ? (e.source === id ? e.target : e.source) : undefined;
}

/** A temperature [K] to 0.01 K; missing or NaN as "--". */
export function formatKelvin(T: number | null | undefined): string {
	if (T === null || T === undefined || !Number.isFinite(T)) return "--";
	return `${T.toFixed(2)} K`;
}
