import { describe, expect, it } from "vitest";
import {
	designArrangement,
	orientationOf,
	sideRun,
	solvedArrangement,
} from "./wallArrangement";

const nodes = [
	{ id: "hi", data: { label: "Hot In" } },
	{ id: "ho", data: { label: "Hot Out" } },
	{ id: "hch", data: { label: "Hot Duct" } },
	{ id: "ci", data: {} },
	{ id: "co", data: {} },
	{ id: "cch", data: {} },
];

describe("sideRun", () => {
	it("follows the declared edge direction, whatever the handles", () => {
		const edges = [
			{ id: "e1", source: "hi", target: "hch" },
			{ id: "e2", source: "hch", target: "ho" },
			// The cold duct is drawn right to left: the edges still say co -> ci.
			{ id: "e3", source: "co", target: "cch" },
			{ id: "e4", source: "cch", target: "ci" },
			{ id: "w", source: "hch", target: "cch", data: { type: "thermal" } },
		];
		expect(sideRun("hch", nodes, edges)).toEqual({
			id: "hch",
			label: "Hot Duct",
			from: "Hot In",
			to: "Hot Out",
		});
		// The thermal edge is not a flow edge: it is neither upstream nor down.
		expect(sideRun("cch", nodes, edges)).toEqual({
			id: "cch",
			label: "cch",
			from: "co",
			to: "ci",
		});
	});

	it("leaves an unconnected end empty", () => {
		expect(sideRun("hch", nodes, [])).toEqual({
			id: "hch",
			label: "Hot Duct",
			from: null,
			to: null,
		});
	});
});

describe("orientation and arrangement", () => {
	it("defaults a saved file without the field to aligned", () => {
		expect(orientationOf({})).toBe("aligned");
		expect(orientationOf(undefined)).toBe("aligned");
		expect(orientationOf({ orientation: "opposed" })).toBe("opposed");
		expect(orientationOf({ orientation: "bogus" })).toBe("aligned");
	});

	it("reads the solved arrangement from the wall's counterflow flag", () => {
		expect(solvedArrangement(undefined)).toBeUndefined();
		expect(solvedArrangement({})).toBeUndefined();
		expect(solvedArrangement({ counterflow: 0 })).toBe("parallel");
		expect(solvedArrangement({ counterflow: 1 })).toBe("counter");
	});

	it("runs forward flow parallel when aligned, counter when opposed", () => {
		expect(designArrangement("aligned")).toBe("parallel");
		expect(designArrangement("opposed")).toBe("counter");
	});
});
