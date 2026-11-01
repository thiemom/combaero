// Flow edges whose handle ids their nodes lack (#901). The canvas draws an edge
// only between handles its two nodes define; the backend solves an edge by its
// endpoints alone. So a file whose edge names a missing handle (a script that
// wrote "flow-source" for a plenum, whose outlets are flow-right/flow-bottom)
// solves, but shows the node disconnected. On load such an edge is moved to the
// node's first handle on that side when every handle there is positional
// (flow-* or without an id: the same connection, drawn elsewhere), and reported
// otherwise -- a named port (an ejector's secondary inlet, a chamber's side
// inlet) is not guessed.

import table from "./nodeHandles.generated.json";

type HandleId = string | null;

interface NodeHandles {
	inlets: HandleId[];
	outlets: HandleId[];
	thermal: string[];
	dynamic: boolean;
}

const HANDLES = table as Record<string, NodeHandles>;

/** Handles that are interchangeable positions of one connection. */
const positional = (h: HandleId): boolean =>
	h === null || h.startsWith("flow-");

interface NodeLike {
	id: string;
	type?: string;
}

interface EdgeLike {
	id: string;
	source: string;
	target: string;
	sourceHandle?: string | null;
	targetHandle?: string | null;
	type?: string;
	data?: { type?: string } | null;
}

export interface HandleRepair<E> {
	edges: E[];
	/** "<edge>: <node> <old> -> <new>" for each end moved. */
	repaired: string[];
	/** "<edge>: <node> has no handle <old>" for each end left as it was. */
	unresolved: string[];
}

export function repairEdgeHandles<E extends EdgeLike>(
	nodes: NodeLike[],
	edges: E[],
): HandleRepair<E> {
	const kind = new Map(nodes.map((n) => [n.id, n.type ?? ""]));
	const repaired: string[] = [];
	const unresolved: string[] = [];
	const out = edges.map((e) => {
		if (
			e.type === "thermal" ||
			e.type === "probe_link" ||
			e.data?.type === "thermal"
		) {
			return e;
		}
		let next = e;
		for (const end of ["source", "target"] as const) {
			const node = end === "source" ? e.source : e.target;
			const spec = HANDLES[kind.get(node) ?? ""];
			if (!spec || spec.dynamic) continue;
			const key = end === "source" ? "sourceHandle" : "targetHandle";
			const side = end === "source" ? spec.outlets : spec.inlets;
			// No id at all is fine: the canvas then uses the node's first
			// handle on that side. Only an id the node lacks is broken.
			if (e[key] === undefined || e[key] === null) continue;
			const have = e[key] as HandleId;
			if (side.length === 0 || side.includes(have)) continue;
			if (side.every(positional)) {
				next = { ...next, [key]: side[0] };
				repaired.push(
					`${e.id}: ${node} ${have ?? "(none)"} -> ${side[0] ?? "(none)"}`,
				);
			} else {
				unresolved.push(`${e.id}: ${node} has no handle ${have ?? "(none)"}`);
			}
		}
		return next;
	});
	return { edges: out, repaired, unresolved };
}
