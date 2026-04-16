import axios from "axios";

const API_BASE_URL = "http://localhost:8000";

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

export interface SolveResponse {
	success: boolean;
	message: string;
	node_results: Record<string, NodeResult>;
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
