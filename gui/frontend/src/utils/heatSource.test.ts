import type { Edge, Node } from "reactflow";
import { describe, expect, it } from "vitest";
import {
	formatHeat,
	formatKelvin,
	HEAT_SOURCE_HANDLE,
	heatedNodeId,
	heatSourceConnectionError,
	heatSourceSummary,
} from "./heatSource";

const node = (id: string, type: string): Node => ({
	id,
	type,
	position: { x: 0, y: 0 },
	data: {},
});

const nodes = [
	node("hs", "heat_source"),
	node("pl", "plenum"),
	node("mc", "momentum_chamber"),
	node("cb", "combustor"),
	node("ch", "channel"),
	node("pb", "pressure_boundary"),
];

const conn = (source: string, target: string) => ({
	source,
	target,
	sourceHandle: HEAT_SOURCE_HANDLE,
	targetHandle: "flow-left",
});

describe("heatSourceConnectionError", () => {
	it("accepts one edge to a plenum, chamber or combustor", () => {
		for (const t of ["pl", "mc", "cb"]) {
			expect(heatSourceConnectionError(conn("hs", t), nodes, [])).toBeNull();
		}
	});

	it("refuses an element or a boundary", () => {
		expect(heatSourceConnectionError(conn("hs", "ch"), nodes, [])).toMatch(
			/plenum, momentum chamber or combustor/,
		);
		expect(heatSourceConnectionError(conn("hs", "pb"), nodes, [])).toMatch(
			/plenum, momentum chamber or combustor/,
		);
	});

	it("refuses a second edge", () => {
		const edges = [{ id: "e", source: "hs", target: "pl" }] as Edge[];
		expect(heatSourceConnectionError(conn("hs", "mc"), nodes, edges)).toMatch(
			/one node only/,
		);
	});

	it("refuses an edge into a heat source", () => {
		expect(
			heatSourceConnectionError(
				{ ...conn("pl", "hs"), sourceHandle: "flow-right" },
				nodes,
				[],
			),
		).toMatch(/takes no inflow/);
	});

	it("leaves other connections alone", () => {
		expect(
			heatSourceConnectionError(
				{ source: "pb", target: "ch", sourceHandle: null, targetHandle: null },
				nodes,
				[],
			),
		).toBeNull();
	});
});

describe("heat source display", () => {
	it("summarises each mode, the target as a total temperature", () => {
		expect(heatSourceSummary({})).toBe("Q = 1000.0 W");
		expect(heatSourceSummary({ mode: "Q", Q: -25000 })).toBe("Q = -25.00 kW");
		expect(heatSourceSummary({ mode: "fraction", fraction: -0.05 })).toBe(
			"f = -5.00 %",
		);
		expect(heatSourceSummary({ mode: "Tt_target", Tt_target: 650 })).toBe(
			"Tt = 650.0 K (total)",
		);
	});

	it("shows an undetermined heat as --", () => {
		expect(formatHeat(Number.NaN)).toBe("--");
		expect(formatHeat(undefined)).toBe("--");
		expect(formatHeat(null)).toBe("--");
		expect(formatKelvin(Number.NaN)).toBe("--");
		expect(formatKelvin(612.345)).toBe("612.35 K");
	});

	it("finds the heated node", () => {
		const edges = [{ id: "e", source: "hs", target: "pl" }] as Edge[];
		expect(heatedNodeId("hs", edges)).toBe("pl");
		expect(heatedNodeId("hs", [])).toBeUndefined();
	});
});
