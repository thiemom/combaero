import axios from "axios";

const API_BASE_URL = "http://localhost:8000";

export interface NetworkGraph {
	nodes: any[];
	edges: any[];
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

export const exportResults = async (graph: NetworkGraph): Promise<void> => {
	const response = await axios.post(`${API_BASE_URL}/export`, graph, {
		responseType: "blob",
	});

	// Create a download link
	const url = window.URL.createObjectURL(new Blob([response.data]));
	const link = document.createElement("a");
	link.href = url;
	link.setAttribute("download", "combaero_results.csv");
	document.body.appendChild(link);
	link.click();
	link.remove();
};
