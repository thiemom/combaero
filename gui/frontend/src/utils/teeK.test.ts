import { describe, expect, it } from "vitest";
import { bassettDesignK, designNetLoss } from "./teeK";

// Same golden values as python/tests/test_constant_k_design_K.py: the
// display mirror must show what the solve uses.
const GOLDEN: ["branch" | "merge", number, number, number, [number, number]][] =
	[
		["branch", 0.5, 4.0, 90.0, [0.0, 3.469266]],
		["branch", 0.2, 1.0, 45.0, [-0.06, 0.707412]],
		["merge", 0.5, 4.0, 90.0, [-0.048904, 3.701096]],
		["merge", 0.3, 1.0, 60.0, [0.149539, -0.250461]],
	];

describe("bassettDesignK", () => {
	it.each(GOLDEN)("%s q=%s psi=%s theta=%s", (fd, q, psi, th, want) => {
		const [ks, kb] = bassettDesignK(fd, q, psi, th);
		expect(ks).toBeCloseTo(want[0], 6);
		expect(kb).toBeCloseTo(want[1], 6);
	});
});

describe("designNetLoss", () => {
	it("is the mass-weighted pair", () => {
		expect(designNetLoss(0.3, 0.2, -0.1)).toBeCloseTo(0.11, 12);
	});
	it("is negative for Bassett's joining pair at 30 deg, equal areas (#779)", () => {
		const [ks, kb] = bassettDesignK("merge", 0.5, 1.0, 30.0);
		expect(designNetLoss(0.5, ks, kb)).toBeLessThan(0);
	});
});
