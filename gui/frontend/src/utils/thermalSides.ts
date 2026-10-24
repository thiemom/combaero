// Which canvas nodes can be a thermal wall's side (#720).
//
// A wall exchanges heat between two convective surfaces, and the solver now
// refuses a side without one. An orifice or area change is a restriction, not
// a surface: neither has a thermal handle any more. A plenum is a state with
// no flow over a surface. An effusion plate owns its
// wall (its bore is in the internal correlation) and film holes cool a wall
// downstream, so neither ever had one.

interface NodeLike {
	id: string;
	type?: string;
}

interface EdgeLike {
	id: string;
	source: string;
	target: string;
	data?: { type?: string } | null;
}

/** Node types whose element has no convective surface. */
export const NO_THERMAL_SIDE: ReadonlySet<string> = new Set([
	"plenum",
	"orifice",
	"area_change",
	"effusion_plate",
	"film_cooling",
]);

/** Drop the thermal walls a saved network hangs on a node with no surface.
 * Their handles are gone, so the canvas could neither draw nor delete them;
 * the solver would refuse them. Returns the kept edges and the dropped ids. */
export function dropSurfacelessWalls<E extends EdgeLike>(
	nodes: NodeLike[],
	edges: E[],
): { edges: E[]; dropped: string[] } {
	const bare = new Set(
		nodes.filter((n) => n.type && NO_THERMAL_SIDE.has(n.type)).map((n) => n.id),
	);
	const dropped: string[] = [];
	const kept = edges.filter((e) => {
		if (e.data?.type !== "thermal") return true;
		if (bare.has(e.source) || bare.has(e.target)) {
			dropped.push(e.id);
			return false;
		}
		return true;
	});
	return { edges: kept, dropped };
}
