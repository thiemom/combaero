import type { Node } from "reactflow";
import { describe, expect, it } from "vitest";
import { validateNetwork } from "./validation";

// The tee branch angle must satisfy 0 < |theta| <= 90 deg, the same rule as
// the backend's MPCETeeData validator. The closures read it only through
// cos/sin, so the check here is the only thing that refuses 900 (which would
// solve as a 180 deg tee) or 0 (the branch on the straight run).

const tee = (theta_deg: number): Node => ({
	id: "tee",
	type: "mpce_tee",
	position: { x: 0, y: 0 },
	data: { theta_deg },
});

const angleErrors = (theta_deg: number) =>
	validateNetwork([tee(theta_deg)]).filter((e) => e.includes("branch angle"));

describe("tee branch angle", () => {
	it.each([90, 45, -90, 1e-3])("accepts %s deg", (theta) => {
		expect(angleErrors(theta)).toHaveLength(0);
	});

	it.each([0, 900, 180, 91, -91])("refuses %s deg", (theta) => {
		expect(angleErrors(theta)).toHaveLength(1);
	});
});
