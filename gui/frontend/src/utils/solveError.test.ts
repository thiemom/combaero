import axios from "axios";
import { describe, expect, it } from "vitest";
import { BACKEND_UNREACHABLE_MESSAGE, solveErrorMessage } from "./solveError";

const axiosError = (response?: unknown) =>
	new axios.AxiosError(
		"request failed",
		"ERR",
		undefined,
		undefined,
		response as never,
	);

describe("solveErrorMessage", () => {
	it("a non-axios error falls back to the generic message", () => {
		const info = solveErrorMessage(new Error("boom"), true);
		expect(info.message).toContain("Check console for details");
		expect(info.isTopologyMismatch).toBeUndefined();
	});

	it("an unreachable backend says how to start it, not a solve error", () => {
		// No HTTP response (connection refused) -> health probe also fails.
		const info = solveErrorMessage(axiosError(undefined), false);
		expect(info.message).toBe(BACKEND_UNREACHABLE_MESSAGE);
	});

	it("unreachable takes precedence even if a stale response detail exists", () => {
		const info = solveErrorMessage(
			axiosError({ data: { detail: "something" } }),
			false,
		);
		expect(info.message).toBe(BACKEND_UNREACHABLE_MESSAGE);
	});

	it("a reachable backend's error detail is surfaced", () => {
		const info = solveErrorMessage(
			axiosError({ data: { detail: "The iteration is not making progress" } }),
			true,
		);
		expect(info.message).toBe(
			"Solve failed: The iteration is not making progress",
		);
		expect(info.isTopologyMismatch).toBe(false);
	});

	it("a topology-mismatch detail sets the flag", () => {
		const info = solveErrorMessage(
			axiosError({ data: { detail: "Topology mismatch: unknowns differ" } }),
			true,
		);
		expect(info.isTopologyMismatch).toBe(true);
	});

	it("a reachable backend error without a detail stays generic", () => {
		const info = solveErrorMessage(axiosError({ data: {} }), true);
		expect(info.message).toContain("Check console for details");
	});

	it("names the export when the export failed", () => {
		const info = solveErrorMessage(
			axiosError({ data: { detail: "Element 'ch' is not fully connected." } }),
			true,
			"Export",
		);
		expect(info.message).toBe(
			"Export failed: Element 'ch' is not fully connected.",
		);
		expect(solveErrorMessage(new Error("x"), true, "Export").message).toBe(
			"Export failed. Check console for details.",
		);
	});
});
