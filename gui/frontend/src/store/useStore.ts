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
	setSolveResults: (results: any) => void;
	speciesMetadata: { names: string[]; molar_masses: number[] } | null;
	fetchSpeciesMetadata: () => Promise<void>;
	displaySettings: string[];
	setDisplaySettings: (settings: string[]) => void;
	unitPreferences: { pressure: "Pa" | "kPa" | "MPa" };
	setPressureUnit: (unit: "Pa" | "kPa" | "MPa") => void;
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
		set({
			edges: addEdge(connection, get().edges),
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

	unitPreferences: { pressure: "kPa" },
	setPressureUnit: (unit: "Pa" | "kPa" | "MPa") => {
		set({ unitPreferences: { pressure: unit } });
	},
}));

export default useStore;
