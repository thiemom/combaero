import { describe, expect, it } from "vitest";
import { cardSize } from "./nodeOrigin";
import {
	resolvedFractions,
	splitterHeight,
	splitterReference,
	validateSplitter,
} from "./splitter";
import { validateNetwork } from "./validation";

describe("flow splitter data (#903)", () => {
	it("resolves the remainder and defaults the reference to it", () => {
		expect(resolvedFractions([0.2, 0.3, null])).toEqual([0.2, 0.3, 0.5]);
		expect(splitterReference({ fractions: [0.2, null, 0.3] })).toBe(1);
		expect(splitterReference({ fractions: [0.2, null], reference: 0 })).toBe(0);
		expect(splitterReference({ fractions: [0.5, 0.5] })).toBeNull();
	});

	it("accepts what the backend accepts", () => {
		expect(validateSplitter({ fractions: [0.5, null] })).toEqual([]);
		// 0.1 + 0.2 + 0.7 is 1 + 2.2e-16: round-off, within the tolerance.
		expect(
			validateSplitter({ fractions: [0.1, 0.2, 0.7], reference: 2 }),
		).toEqual([]);
	});

	it.each([
		[{ fractions: [0.5, 0.4], reference: 0 }, "not 1"],
		[{ fractions: [0.5, 0.5, null] }, "nothing left"],
		[{ fractions: [null, null] }, "at most one"],
		[{ fractions: [0, null] }, "in (0, 1)"],
		[{ fractions: [0.5, 0.5] }, "choose the reference"],
		[{ fractions: [0.5, null], reference: 5 }, "does not exist"],
	])("refuses %j", (data, match) => {
		expect(validateSplitter(data).join("; ")).toContain(match);
	});

	it("is checked before a solve", () => {
		const errors = validateNetwork([
			{
				id: "s",
				type: "flow_splitter",
				position: { x: 0, y: 0 },
				data: { fractions: [0.5, 0.4], reference: 0 },
			},
		]);
		expect(errors.some((e) => e.startsWith("s: Flow splitter"))).toBe(true);
	});

	it("grows with its outlets", () => {
		expect(splitterHeight(2)).toBe(56);
		expect(splitterHeight(6)).toBe(6 * 22 + 12);
		const node = {
			id: "s",
			type: "flow_splitter",
			position: { x: 0, y: 0 },
			data: { fractions: [0.1, 0.1, 0.1, 0.1, null] },
		};
		expect(cardSize(node)).toEqual([140, splitterHeight(5)]);
	});
});
