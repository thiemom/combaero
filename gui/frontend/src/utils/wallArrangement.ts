// A thermal wall's arrangement (#627, #629), read from what the model
// declares: each side's run (the element's from_node -> to_node, which the
// backend takes from the flow edges' direction, never from handle sides) and
// the wall's orientation. With the solved flow signs these decide parallel
// or counter flow; nothing here is inferred from where a handle is drawn.

export type WallOrientation = "aligned" | "opposed";

interface GraphNode {
	id: string;
	data?: { label?: string };
}

interface GraphEdge {
	id: string;
	source: string;
	target: string;
	data?: { type?: string };
}

export interface SideRun {
	id: string;
	label: string;
	from: string | null;
	to: string | null;
}

const nameOf = (id: string, nodes: GraphNode[]): string =>
	nodes.find((n) => n.id === id)?.data?.label || id;

// The declared run of one wall side: upstream and downstream along the
// element's flow edges (thermal edges carry no flow and are skipped).
export function sideRun(
	id: string,
	nodes: GraphNode[],
	edges: GraphEdge[],
): SideRun {
	const flow = edges.filter((e) => e.data?.type !== "thermal");
	const inEdge = flow.find((e) => e.target === id);
	const outEdge = flow.find((e) => e.source === id);
	return {
		id,
		label: nameOf(id, nodes),
		from: inEdge ? nameOf(inEdge.source, nodes) : null,
		to: outEdge ? nameOf(outEdge.target, nodes) : null,
	};
}

export const orientationOf = (data?: {
	orientation?: string;
}): WallOrientation =>
	data?.orientation === "opposed" ? "opposed" : "aligned";

// The arrangement the solve reached: the wall reports counterflow as 1/0
// (ThermalWall.is_counterflow). Undefined before a solve.
export function solvedArrangement(result?: {
	counterflow?: number;
}): "parallel" | "counter" | undefined {
	const cf = result?.counterflow;
	if (cf === undefined || cf === null) return undefined;
	return cf >= 0.5 ? "counter" : "parallel";
}

// What forward flow on both sides would run as: the declaration alone.
export const designArrangement = (
	o: WallOrientation,
): "parallel" | "counter" => (o === "opposed" ? "counter" : "parallel");

// Arrow glyphs for the edge marker: aligned runs side by side the same way,
// opposed the opposite way.
export const orientationGlyph = (o: WallOrientation): string =>
	o === "opposed" ? "⇄" : "⇉";
