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

describe("tee loss factor", () => {
	const ideal = (loss_factor: number): Node => ({
		id: "tee",
		type: "mpce_tee",
		position: { x: 0, y: 0 },
		data: { junction_model: "ideal_momentum", loss_factor },
	});
	const errs = (f: number) =>
		validateNetwork([ideal(f)]).filter((e) => e.includes("loss factor"));

	it.each([0.25, 1.0, 3.0])("accepts %s", (f) => {
		expect(errs(f)).toHaveLength(0);
	});

	it.each([0.2, 0, -1])("refuses %s", (f) => {
		expect(errs(f)).toHaveLength(1);
	});
});

// Ejector recovery multiplies the lossless mixed stagnation pressure: above 1
// is a pressure gain beyond lossless mixing, the same bound as the backend's
// EjectorData field and EjectorElement (#769).
describe("ejector recovery efficiency", () => {
	const ejector = (recovery_efficiency: number): Node => ({
		id: "ej",
		type: "ejector",
		position: { x: 0, y: 0 },
		data: { recovery_efficiency },
	});
	const errs = (r: number) =>
		validateNetwork([ejector(r)]).filter((e) =>
			e.includes("recovery_efficiency"),
		);

	it.each([1.0, 0.9, 1e-3])("accepts %s", (r) => {
		expect(errs(r)).toHaveLength(0);
	});

	it.each([1.0000001, 2, 0, -0.5, Number.NaN])("refuses %s", (r) => {
		expect(errs(r)).toHaveLength(1);
	});
});
