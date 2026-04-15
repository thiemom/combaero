import {
	addEdge,
	applyEdgeChanges,
	applyNodeChanges,
	type Connection,
	type Edge,
	type EdgeChange,
	MarkerType,
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
	unitPreferences: {
		pressure: "Pa" | "kPa" | "MPa";
		length: "m" | "mm";
		area: "m²" | "mm²";
	};
	setPressureUnit: (unit: "Pa" | "kPa" | "MPa") => void;
	setLengthUnit: (unit: "m" | "mm") => void;
	setAreaUnit: (unit: "m²" | "mm²") => void;
	saveNetwork: (filename?: string) => void | Promise<void>;
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
			solveResults: null,
		});
	},

	onEdgesChange: (changes: EdgeChange[]) => {
		set({
			edges: applyEdgeChanges(changes, get().edges),
			solveResults: null,
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
			markerEnd: isThermal
				? {
						type: MarkerType.ArrowClosed,
						color: "rgba(255, 152, 0, 0.45)",
					}
				: undefined,
			data: isThermal
				? { type: "thermal", thickness: 0.003, conductivity: 20.0, area: 0.05 }
				: null,
		} as Edge;

		set({
			edges: addEdge(edge, get().edges),
			solveResults: null,
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
			solveResults: null,
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
			solveResults: null,
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

	unitPreferences: { pressure: "kPa", length: "m", area: "m²" },
	setPressureUnit: (unit: "Pa" | "kPa" | "MPa") => {
		set((state) => ({
			unitPreferences: { ...state.unitPreferences, pressure: unit },
		}));
	},
	setLengthUnit: (unit: "m" | "mm") => {
		set((state) => ({
			unitPreferences: { ...state.unitPreferences, length: unit },
		}));
	},
	setAreaUnit: (unit: "m²" | "mm²") => {
		set((state) => ({
			unitPreferences: { ...state.unitPreferences, area: unit },
		}));
	},

	saveNetwork: async (filename?: string) => {
		const { nodes, edges, solverSettings } = get();
		const data = { nodes, edges, solverSettings, version: "1.1" };
		const content = JSON.stringify(data, null, 2);

		// Try File System Access API (Native "Save As" Dialog)
		if ("showSaveFilePicker" in window) {
			try {
				const handle = await (window as any).showSaveFilePicker({
					suggestedName: filename
						? `${filename.trim()}.json`
						: `combaero_network_${new Date().toISOString().split("T")[0]}.json`,
					types: [
						{
							description: "CombAero Network JSON",
							accept: { "application/json": [".json"] },
						},
					],
				});
				const writable = await handle.createWritable();
				await writable.write(content);
				await writable.close();
				return;
			} catch (err) {
				// AbortError is perfectly normal (user clicked Cancel)
				if ((err as Error).name === "AbortError") return;
				console.warn(
					"Native Save Picker failed, falling back to legacy download",
					err,
				);
			}
		} else {
			console.info(
				"File System Access API not supported in this browser. Using legacy download mode.",
			);
		}

		// Fallback: the classic <a download> method
		const blob = new Blob([content], {
			type: "application/json",
		});
		const url = URL.createObjectURL(blob);
		const link = document.createElement("a");
		link.href = url;

		// Smarter filename for sorting in legacy mode (adds timestamp if collision likely)
		const timestamp = new Date()
			.toLocaleTimeString("en-GB", {
				hour: "2-digit",
				minute: "2-digit",
				second: "2-digit",
			})
			.replace(/:/g, "-");

		link.download = filename
			? `${filename.trim()}.json`
			: `combaero_network_${new Date().toISOString().split("T")[0]}_${timestamp}.json`;
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
