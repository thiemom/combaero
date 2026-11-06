import { describe, expect, it } from "vitest";
import { FLOW_REVERSAL_M_SCALE, solvedFlowAsMref } from "./fractionLoss";

describe("solvedFlowAsMref (#1004)", () => {
	it("is the solved |m_dot| after a converged solve", () => {
		expect(solvedFlowAsMref({ result: { success: true, m_dot: 0.25 } })).toBe(
			0.25,
		);
		expect(solvedFlowAsMref({ result: { success: true, m_dot: -0.25 } })).toBe(
			0.25,
		);
	});

	it("is off with no result or a failed solve", () => {
		expect(solvedFlowAsMref({})).toBeNull();
		expect(solvedFlowAsMref({ result: null })).toBeNull();
		expect(
			solvedFlowAsMref({ result: { success: false, m_dot: 0.25 } }),
		).toBeNull();
	});

	it("is off where the element carries no flow (its zero-flow blend)", () => {
		expect(
			solvedFlowAsMref({
				result: { success: true, m_dot: 0.5 * FLOW_REVERSAL_M_SCALE },
			}),
		).toBeNull();
		expect(
			solvedFlowAsMref({ result: { success: true, m_dot: null } }),
		).toBeNull();
	});

	it("is off below a set m_ref, where the loss was faded", () => {
		expect(
			solvedFlowAsMref({ m_ref: 0.5, result: { success: true, m_dot: 0.25 } }),
		).toBeNull();
		expect(
			solvedFlowAsMref({ m_ref: 0.1, result: { success: true, m_dot: 0.25 } }),
		).toBe(0.25);
	});
});
