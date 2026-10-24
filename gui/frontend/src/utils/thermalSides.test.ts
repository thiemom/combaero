import { describe, expect, it } from "vitest";
import { dropSurfacelessWalls } from "./thermalSides";

const nodes = [
	{ id: "ch", type: "channel" },
	{ id: "ch2", type: "channel" },
	{ id: "or", type: "orifice" },
	{ id: "ac", type: "area_change" },
	{ id: "w", type: "wall" },
];
const thermal = (id: string, source: string, target: string) => ({
	id,
	source,
	target,
	data: { type: "thermal" },
});

describe("dropSurfacelessWalls", () => {
	it("keeps a wall between two channels and every flow edge", () => {
		const edges = [
			thermal("t1", "ch", "ch2"),
			{ id: "f1", source: "ch", target: "or" },
		];
		const r = dropSurfacelessWalls(nodes, edges);
		expect(r.edges.map((e) => e.id)).toEqual(["t1", "f1"]);
		expect(r.dropped).toEqual([]);
	});

	it("drops a wall with an orifice or area-change side, either end", () => {
		const edges = [
			thermal("t1", "or", "ch"),
			thermal("t2", "ch", "ac"),
			thermal("t3", "w", "or"),
			thermal("t4", "ch", "w"),
			{ id: "f1", source: "or", target: "ac", data: { type: "flow" } },
		];
		const r = dropSurfacelessWalls(nodes, edges);
		expect(r.edges.map((e) => e.id)).toEqual(["t4", "f1"]);
		expect(r.dropped).toEqual(["t1", "t2", "t3"]);
	});
});
