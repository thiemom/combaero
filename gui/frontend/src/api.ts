import axios from "axios";

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
}

export const solveNetwork = async (
	graph: NetworkGraph,
): Promise<SolveResponse> => {
	const response = await axios.post(`${API_BASE_URL}/solve`, graph);
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
