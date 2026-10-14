import { describe, expect, it } from "vitest";
import {
	areaFromDiameter,
	diameterFromArea,
	geometryView,
	storedFrom,
} from "./geometry";

const A01 = areaFromDiameter(0.1); // 0.00785398...

describe("geometryView", () => {
	it("locks an inherited face and shows D, Dh and A from it", () => {
		const v = geometryView({
			stores: "A",
			userValue: null,
			inherited: { value: A01, source: "channel 'ch'" },
		});
		expect(v.locked).toBe(true);
		expect(v.canInherit).toBe(true);
		expect(v.D).toBeCloseTo(0.1, 12);
		expect(v.A).toBeCloseTo(A01, 15);
		expect(v.Dh).toBeCloseTo(0.1, 12); // circular
		expect(v.userSet).toBe(false);
		expect(v.source).toBe("channel 'ch'");
	});

	it("unlocks once the user sets the stored value", () => {
		const v = geometryView({
			stores: "D",
			userValue: 0.2,
			inherited: { value: 0.1, source: "channel 'ch'" },
		});
		expect(v.locked).toBe(false);
		expect(v.userSet).toBe(true);
		expect(v.D).toBe(0.2);
		expect(v.A).toBeCloseTo(areaFromDiameter(0.2), 15);
		expect(v.source).toBe("user");
	});

	it("a Dh override alone also opts the face out", () => {
		const v = geometryView({
			stores: "D",
			userValue: null,
			inherited: { value: 0.1, source: "channel 'ch'" },
			userDh: 0.05,
		});
		expect(v.locked).toBe(false);
		expect(v.D).toBe(0.1); // still the inherited D
		expect(v.Dh).toBe(0.05);
		expect(v.userSetDh).toBe(true);
	});

	it("nothing to inherit: editable, no opt-out box, values unknown", () => {
		const v = geometryView({
			stores: "D",
			userValue: undefined,
			inherited: { value: null, source: "" },
		});
		expect(v.canInherit).toBe(false);
		expect(v.locked).toBe(false);
		expect(v.D).toBeNull();
		expect(v.A).toBeNull();
		expect(v.Dh).toBeNull();
	});

	it("a non-inheritable face (orifice bore) is never locked", () => {
		const v = geometryView({
			stores: "D",
			userValue: null,
			inherited: { value: 0.08, source: "default" },
			inheritable: false,
		});
		expect(v.canInherit).toBe(false);
		expect(v.locked).toBe(false);
	});

	it("reports the solver's own Dh for an element that inherits it", () => {
		// A chamber with a user area still takes Dh from the upstream channel.
		const v = geometryView({
			stores: "A",
			userValue: 0.05,
			inherited: { value: A01, source: "channel 'ch'" },
			userDh: null,
			inheritedDh: { value: 0.1, source: "channel 'ch'" },
		});
		expect(v.D).toBeCloseTo(diameterFromArea(0.05), 12);
		expect(v.Dh).toBe(0.1);
		expect(v.dhSource).toBe("channel 'ch'");
	});
});

describe("storedFrom", () => {
	it("converts an entered D or A to the stored quantity", () => {
		expect(storedFrom("D", "D", 0.1)).toBe(0.1);
		expect(storedFrom("A", "A", A01)).toBe(A01);
		expect(storedFrom("A", "D", 0.1)).toBeCloseTo(A01, 15);
		expect(storedFrom("D", "A", A01)).toBeCloseTo(0.1, 12);
	});
});
