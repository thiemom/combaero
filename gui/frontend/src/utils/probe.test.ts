import { describe, expect, it } from "vitest";
import { isConvergedResult, probeKeyOptions, resolveKey } from "./diagnostics";
import { PORT_ARMS, portArms, withPortArm } from "./portArms";
import { expectedResultKeys } from "./quantities";

describe("probeKeyOptions", () => {
	it("offers every finite number of a node result, state keys prefixed", () => {
		const keys = probeKeyOptions(
			{
				state: { P: 1e5, T: 300, X: [0.79, 0.21], mach: Number.NaN },
				success: true,
			},
			"plenum",
		);
		expect(keys.sort()).toEqual(["state.P", "state.T"]);
	});

	it("offers an element's top-level numbers, not its strings or flags", () => {
		const keys = probeKeyOptions(
			{ m_dot: 0.1, ref_location: "inlet", success: true, rows_Nu: [1, 2] },
			"channel",
		);
		expect(keys).toEqual(["m_dot"]);
	});

	it("before a solve, offers what a solve gives the target's kind", () => {
		expect(probeKeyOptions(undefined, "area_change")).toEqual(
			expectedResultKeys("area_change"),
		);
		expect(probeKeyOptions(undefined, "area_change")).toContain("zeta");
		// The old hand-kept fallback offered P / T / Re for every kind; an
		// area change reports neither a node state nor a node P.
		expect(probeKeyOptions(undefined, "area_change")).not.toContain("P");
		expect(probeKeyOptions(undefined, undefined)).toEqual([]);
	});
});

describe("resolveKey", () => {
	const nodeResult = { state: { P: 2e5 }, P: 1 };
	it("reads state.<k> from the state first", () => {
		expect(resolveKey(nodeResult, "state.P")).toBe(2e5);
		expect(resolveKey(nodeResult, "P")).toBe(1);
	});
	it("returns undefined for a string or non-finite value", () => {
		expect(resolveKey({ tap: "Flange" }, "tap")).toBeUndefined();
		expect(resolveKey({ zeta: Number.NaN }, "zeta")).toBeUndefined();
		expect(resolveKey(undefined, "zeta")).toBeUndefined();
	});
});

describe("isConvergedResult", () => {
	it("is false only for success: false", () => {
		expect(isConvergedResult({ success: false })).toBe(false);
		expect(isConvergedResult({ success: true })).toBe(true);
		expect(isConvergedResult({ Q: 1 })).toBe(true);
		expect(isConvergedResult(undefined)).toBe(false);
	});
});

describe("portArms", () => {
	// The port order is pinned against a solve in
	// python/tests/test_gui_port_arms.py.
	it("numbers a tee's ports inlets first", () => {
		expect(portArms("mpce_tee", { flow_direction: "branch" })).toEqual([
			"C",
			"S",
			"B",
		]);
		expect(portArms("mpce_tee", { flow_direction: "merge" })).toEqual([
			"S",
			"B",
			"C",
		]);
		expect(portArms("mpce_tee", {})).toEqual(PORT_ARMS.mpce_tee_branch);
		expect(portArms("ejector", {})).toEqual(["P", "S", "O"]);
		expect(portArms("channel", {})).toBeNull();
	});

	it("labels port keys with the arm, and nothing else", () => {
		const merge = { flow_direction: "merge" };
		expect(withPortArm("p0", "port_0_P", "mpce_tee", merge)).toBe("p0 (S)");
		expect(withPortArm("p2", "port_2_m_dot", "mpce_tee", merge)).toBe("p2 (C)");
		expect(withPortArm("p1", "port_1_T", "ejector", {})).toBe("p1 (S)");
		expect(withPortArm("K", "K_branch", "mpce_tee", merge)).toBe("K");
		expect(withPortArm("p0", "port_0_P", "channel", {})).toBe("p0");
		expect(withPortArm("p9", "port_9_P", "mpce_tee", merge)).toBe("p9");
	});
});
