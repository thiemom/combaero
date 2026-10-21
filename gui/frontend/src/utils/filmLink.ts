// A thermal wall's film (#386): which of the wall's two ducts the film from a
// provider (film-cooling holes, an effusion plate) covers by default.
//
// The stored side stays "a" (the wall edge's source) or "b" (its target); the
// GUI names the duct instead, because "A/B" mean nothing on the canvas and
// "hot/cold" follow the solved heat flow, so they swap with its sign.

interface EdgeLike {
	source: string;
	target: string;
	data?: { type?: string } | null;
}

// The film link, tag and handle: orange, a shade off the wall so they read apart.
export const FILM_COLOR = "#f97316";

// The holes' drag handle: dropped on a wall's label box it sets the film; it
// never connects to another handle.
export const FILM_HANDLE = "film-source";

const isFlow = (e: EdgeLike): boolean => e.data?.type !== "thermal";

/** The side whose duct the provider's flow feeds: the duct directly
 * downstream of it, or downstream of the node it discharges into. Falls back
 * to "a" when neither side is (the solver then refuses a wrong pick). */
export function defaultFilmSide(
	wall: { source: string; target: string },
	providerId: string,
	edges: EdgeLike[],
): "a" | "b" {
	const outlets = new Set(
		edges
			.filter((e) => isFlow(e) && e.source === providerId)
			.map((e) => e.target),
	);
	const fed = (duct: string): boolean =>
		outlets.has(duct) ||
		edges.some((e) => isFlow(e) && e.target === duct && outlets.has(e.source));
	if (!fed(wall.source) && fed(wall.target)) return "b";
	return "a";
}

/** The canvas node id of the duct a film side covers. */
export function filmDuctId(
	wall: { source: string; target: string },
	side: string | undefined,
): string {
	return side === "b" ? wall.target : wall.source;
}
