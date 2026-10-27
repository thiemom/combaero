import axios from "axios";

export interface SolveErrorInfo {
	message: string;
	isTopologyMismatch?: boolean;
}

// Shown when the solve request never reached the backend (no HTTP response, or
// the dev proxy could not connect) -- almost always the backend not running,
// not a problem with the network itself.
export const BACKEND_UNREACHABLE_MESSAGE =
	"Cannot reach the solver backend at localhost:8000. Is the server running? " +
	"Start it with `./scripts/run-gui.sh` (or `combaero-gui`), then solve again.";

/**
 * Map a failed /solve (or /export: `action` "Export") request to a
 * user-facing message.
 *
 * `backendReachable` is the result of a /health probe taken after the failure:
 * it separates "the backend is down" (show how to start it) from "the backend
 * answered with an error" (show that error's detail).
 */
export const solveErrorMessage = (
	err: unknown,
	backendReachable: boolean,
	action = "Solve",
): SolveErrorInfo => {
	const generic = `${action} failed. Check console for details.`;
	if (!axios.isAxiosError(err)) {
		return { message: generic };
	}
	if (!backendReachable) {
		return { message: BACKEND_UNREACHABLE_MESSAGE };
	}
	const detail = err.response?.data?.detail as string | undefined;
	if (!detail) {
		return { message: generic };
	}
	return {
		message: `${action} failed: ${detail}`,
		isTopologyMismatch:
			detail.includes("Topology mismatch") ||
			detail.includes("Initial guess mismatch") ||
			detail.includes("x0 length does not match"),
	};
};
