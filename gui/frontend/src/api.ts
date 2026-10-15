import axios from "axios";
import type { IncompressibleFlag } from "./utils/incompressible";
import type { Ledger, SecondLawFlag } from "./utils/secondLaw";

const API_BASE_URL = "";

export interface NetworkGraph {
	nodes: any[];
	edges: any[];
	solver_settings?: Record<string, any>;
}

export interface NodeResult {
	T: number;
	P: number;
	Y: number[];
	success: boolean;
}

export interface ConvergencePoint {
	eval: number;
	t: number;
	norm: number;
}

export interface WorstResidual {
	name: string;
	residual: number;
}

export interface SolveResponse {
	success: boolean;
	message: string;
	final_norm?: number;
	node_results: Record<string, NodeResult>;
	element_results?: Record<string, any>;
	convergence_history?: ConvergencePoint[];
	worst_residuals?: WorstResidual[];
	solver_settings_used?: Record<string, unknown>;
	lm_started_at_eval?: number | null;
	second_law_failures?: SecondLawFlag[];
	incompressible_warnings?: IncompressibleFlag[];
	ledger?: Ledger | Record<string, never>;
}

export const solveNetwork = async (
	graph: NetworkGraph,
): Promise<SolveResponse> => {
	const response = await axios.post(`${API_BASE_URL}/solve`, graph);
	return response.data;
};

// --- geometry-inheritance preview (resolve without solving, #559) ----------

export interface PreviewField {
	value: number | null;
	// "user" | "inherited from <id>" | "channel '<id>'" | "default" | ""
	source: string;
}

export interface PreviewItem {
	status: "resolved" | "unresolved";
	reason: string;
	fields: Record<string, PreviewField>;
}

export interface ResolvePreviewResponse {
	items: Record<string, PreviewItem>;
}

export const resolvePreview = async (
	graph: NetworkGraph,
): Promise<ResolvePreviewResponse> => {
	const response = await axios.post(`${API_BASE_URL}/resolve_preview`, graph);
	return response.data;
};

export const exportResults = async (graph: NetworkGraph): Promise<Blob> => {
	const response = await axios.post(`${API_BASE_URL}/export`, graph, {
		responseType: "blob",
	});
	return response.data;
};

export const fetchMaterials = async (): Promise<string[]> => {
	const response = await axios.get(`${API_BASE_URL}/metadata/materials`);
	return response.data.names || [];
};

export const checkContinuationAvailable = async (): Promise<boolean> => {
	const response = await axios.get(
		`${API_BASE_URL}/solver/continuation_available`,
	);
	return response.data.available;
};

// Probe used after a failed request to tell "backend is down" apart from
// "backend answered with an error". Returns false on any failure (no response,
// connection refused, timeout).
export const isBackendReachable = async (): Promise<boolean> => {
	try {
		await axios.get(`${API_BASE_URL}/health`, { timeout: 3000 });
		return true;
	} catch {
		return false;
	}
};
