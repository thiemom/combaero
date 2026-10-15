import { describe, expect, it } from "vitest";
import {
	type IncompressibleFlag,
	summary,
	warningsFor,
} from "./incompressible";

const flags: IncompressibleFlag[] = [
	{ element_id: "ch1", owner_id: "ch1", mach: 0.31, limit: 0.2 },
	{ element_id: "ln__s0", owner_id: "ln", mach: 0.24, limit: 0.2 },
];

describe("incompressible validity (#588)", () => {
	it("attributes each flag to its canvas node", () => {
		expect(warningsFor("ln", flags).map((f) => f.element_id)).toEqual([
			"ln__s0",
		]);
		expect(warningsFor("other", flags)).toEqual([]);
		expect(warningsFor("ln", undefined)).toEqual([]);
	});

	it("summarises with the limit, or says nothing", () => {
		expect(summary(flags)).toBe("2 elements above Ma 0.2 (incompressible)");
		expect(summary([flags[0]])).toBe("1 element above Ma 0.2 (incompressible)");
		expect(summary([])).toBeNull();
		expect(summary(undefined)).toBeNull();
	});
});
