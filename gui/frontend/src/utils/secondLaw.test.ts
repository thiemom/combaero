import { describe, expect, it } from "vitest";
import {
	balanceIsOpen,
	describeElement,
	displayName,
	failuresFor,
	type Ledger,
} from "./secondLaw";

const flag = (element_id: string, owner_id: string) => ({
	element_id,
	owner_id,
	sigma_norm: -1e-3,
	sigma_norm_high: -1e-3,
});

const ledger = (mass: number, energy: number): Ledger => ({
	streams: [],
	m_dot_in: 0.1,
	m_dot_out: 0.1,
	mass_closure_rel: mass,
	H_in: 3e4,
	H_out: 3e4,
	Q_given: 0,
	Q_wall_out: 0,
	Q_withheld: 0,
	energy_closure: energy * 3e4,
	energy_closure_rel: energy,
});

describe("second law (#545)", () => {
	it("flags a node for the failures it owns, hidden ones included", () => {
		const fs = [flag("__tee_step__t1_straight", "t1"), flag("ch", "ch")];
		expect(failuresFor("t1", fs).map((f) => f.element_id)).toEqual([
			"__tee_step__t1_straight",
		]);
		expect(failuresFor("other", fs)).toEqual([]);
		expect(failuresFor("t1", undefined)).toEqual([]);
	});
	it("names a hidden port step by its port", () => {
		expect(describeElement(flag("__tee_step__t1_branch", "t1"))).toBe(
			"the area step at its branch port",
		);
		expect(describeElement(flag("ch", "ch"))).toBe("this element");
		expect(describeElement(flag("ln__s0", "ln"))).toBe("ln__s0");
	});
	it("calls the balance open past 1e-6 in mass or energy", () => {
		expect(balanceIsOpen(ledger(1e-9, 1e-8))).toBe(false);
		expect(balanceIsOpen(ledger(0, 2e-6))).toBe(true);
		expect(balanceIsOpen(ledger(-2e-6, 0))).toBe(true);
		expect(balanceIsOpen(undefined)).toBe(false);
	});
	it("names a node by its label, else its kind and the end of its id", () => {
		const nodes = [
			{ id: "node_1791967900397", type: "mass_boundary", data: {} },
			{ id: "node_2", type: "plenum", data: { label: "Supply" } },
		];
		expect(displayName("node_1791967900397", nodes)).toBe("mass boundary 0397");
		expect(displayName("node_2", nodes)).toBe("Supply");
		expect(displayName("gone", nodes)).toBe("gone");
	});
});
