import {
	addEdge,
	applyEdgeChanges,
	applyNodeChanges,
	type Connection,
	type Edge,
	type EdgeChange,
	type Node,
	type NodeChange,
} from "reactflow";
import { create } from "zustand";

export interface RFState {
	nodes: Node[];
	edges: Edge[];
	solveResults: any | null;
	onNodesChange: (changes: NodeChange[]) => void;
	onEdgesChange: (changes: EdgeChange[]) => void;
	onConnect: (connection: Connection) => void;
	setNodes: (nodes: Node[]) => void;
	setEdges: (edges: Edge[]) => void;
	updateNodeData: (nodeId: string, data: any) => void;
	updateEdgeData: (edgeId: string, data: any) => void;
	setSolveResults: (results: any) => void;
	speciesMetadata: { names: string[]; molar_masses: number[] } | null;
	fetchSpeciesMetadata: () => Promise<void>;
	displaySettings: string[];
	setDisplaySettings: (settings: string[]) => void;
	solverSettings: {
		global_regime: "incompressible" | "compressible";
		init_strategy: "default" | "incompressible_warmstart" | "homotopy";
	};
	updateSolverSettings: (settings: Partial<RFState["solverSettings"]>) => void;
	unitPreferences: { pressure: "Pa" | "kPa" | "MPa" };
	setPressureUnit: (unit: "Pa" | "kPa" | "MPa") => void;
	saveNetwork: (filename?: string) => void;
	loadNetwork: (data: any) => void;
	dismissSolveStatus: () => void;
}

const useStore = create<RFState>((set, get) => ({
	nodes: [],
	edges: [],
	solveResults: null,

	onNodesChange: (changes: NodeChange[]) => {
		set({
			nodes: applyNodeChanges(changes, get().nodes),
		});
	},

	onEdgesChange: (changes: EdgeChange[]) => {
		set({
			edges: applyEdgeChanges(changes, get().edges),
		});
	},

	onConnect: (connection: Connection) => {
		const isThermal =
			connection.sourceHandle?.includes("thermal") ||
			connection.targetHandle?.includes("thermal");

		const edge: Edge = {
			...connection,
			id: `edge_${Date.now()}`,
			type: isThermal ? "thermal" : "default",
			data: isThermal
				? { type: "thermal", thickness: 0.003, conductivity: 20.0, area: 0.05 }
				: null,
		} as Edge;

		set({
			edges: addEdge(edge, get().edges),
		});
	},

	setNodes: (nodes: Node[]) => {
		set({ nodes });
	},

	setEdges: (edges: Edge[]) => {
		set({ edges });
	},

	updateNodeData: (nodeId: string, data: any) => {
		set({
			nodes: get().nodes.map((node) => {
				if (node.id === nodeId) {
					return { ...node, data: { ...node.data, ...data } };
				}
				return node;
			}),
		});
	},
	updateEdgeData: (edgeId: string, data: any) => {
		set({
			edges: get().edges.map((edge) => {
				if (edge.id === edgeId) {
					return { ...edge, data: { ...edge.data, ...data } };
				}
				return edge;
			}),
		});
	},

	setSolveResults: (results: any) => {
		set({ solveResults: results });

		// Also inject results into node data for reactive UI
		if (results) {
			set({
				nodes: get().nodes.map((node) => {
					// Check for node results
					if (results.node_results?.[node.id]) {
						return {
							...node,
							data: { ...node.data, result: results.node_results[node.id] },
						};
					}
					// Check for element results (since elements are nodes now)
					if (results.element_results?.[node.id]) {
						return {
							...node,
							data: { ...node.data, result: results.element_results[node.id] },
						};
					}
					return node;
				}),
				edges: get().edges.map((edge) => {
					if (results.edge_results?.[edge.id]) {
						return {
							...edge,
							data: { ...edge.data, result: results.edge_results[edge.id] },
						};
					}
					return edge;
				}),
			});
		}
	},

	speciesMetadata: null,
	fetchSpeciesMetadata: async () => {
		try {
			const res = await fetch("http://localhost:8000/metadata/species");
			const data = await res.json();
			set({ speciesMetadata: data });
		} catch (error) {
			console.error("Failed to fetch species metadata:", error);
		}
	},

	displaySettings: ["T", "P", "m_dot", "rho"],
	setDisplaySettings: (settings: string[]) => {
		set({ displaySettings: settings });
	},
	solverSettings: {
		global_regime: "compressible",
		init_strategy: "default",
	},
	updateSolverSettings: (settings: Partial<RFState["solverSettings"]>) => {
		set({
			solverSettings: { ...get().solverSettings, ...settings },
		});
	},

	unitPreferences: { pressure: "kPa" },
	setPressureUnit: (unit: "Pa" | "kPa" | "MPa") => {
		set({ unitPreferences: { pressure: unit } });
	},

	saveNetwork: (filename?: string) => {
		const { nodes, edges, solverSettings } = get();
		const data = { nodes, edges, solverSettings, version: "1.1" };
		const blob = new Blob([JSON.stringify(data, null, 2)], {
			type: "application/json",
		});
		const url = URL.createObjectURL(blob);
		const link = document.createElement("a");
		link.href = url;
		link.download = filename
			? `${filename}.json`
			: `combaero_network_${new Date().toISOString().split("T")[0]}.json`;
		document.body.appendChild(link);
		link.click();
		document.body.removeChild(link);
		URL.revokeObjectURL(url);
	},

	loadNetwork: (data: any) => {
		if (data.nodes && data.edges) {
			set({
				nodes: data.nodes.map((n: any) => ({
					...n,
					data: { ...n.data, result: undefined },
				})),
				edges: data.edges,
				solverSettings: data.solverSettings || get().solverSettings,
				solveResults: null,
			});
		}
	},

	dismissSolveStatus: () => {
		set({ solveResults: null });
	},
}));

export default useStore;
