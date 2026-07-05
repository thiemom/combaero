import ReactFlow, {
	Background,
	Controls,
	MiniMap,
	useReactFlow,
} from "reactflow";
import "reactflow/dist/style.css";
import { useCallback } from "react";
import { useCopyPaste } from "../hooks/useCopyPaste";
import useStore, { stripResults } from "../store/useStore";
import { edgeTypes, nodeTypes } from "./flowTypes";
import SolverStatusBadge from "./SolverStatusBadge";

const NetworkCanvas = () => {
	const { nodes, edges, setNodes, onNodesChange, onEdgesChange, onConnect } =
		useStore();
	useCopyPaste();
	const { screenToFlowPosition } = useReactFlow();

	const onDragOver = useCallback((event: React.DragEvent) => {
		event.preventDefault();
		event.dataTransfer.dropEffect = "move";
	}, []);

	const onDrop = useCallback(
		(event: React.DragEvent) => {
			event.preventDefault();

			const type = event.dataTransfer.getData("application/reactflow");

			if (typeof type === "undefined" || !type) {
				return;
			}

			// Map client drop position to flow coordinates based on current zoom/pan panning
			const position = screenToFlowPosition({
				x: event.clientX,
				y: event.clientY,
			});

			const id = `node_${Date.now()}`;
			let data = { id };

			if (type === "mass_boundary") {
				data = {
					...data,
					m_dot: 1.0,
					Tt: 300,
					composition: { source: "humid_air", mode: "mole" },
				} as any;
			} else if (type === "pressure_boundary") {
				data = {
					...data,
					Pt: 101325,
					Tt: 300,
					composition: { source: "humid_air", mode: "mole" },
				} as any;
			} else if (type === "channel") {
				data = {
					...data,
					L: 1.0,
					D: 0.1,
					roughness: 1e-5,
					Nu_multiplier: 1.0,
					f_multiplier: 1.0,
				} as any;
			} else if (type === "orifice") {
				data = { ...data, diameter: 0.08, Cd: 0.6 } as any;
			} else if (type === "discrete_loss") {
				data = {
					...data,
					correlation_type: "constant_fraction",
					xi: 0.03,
					zeta: 1.0,
					k: 0.001,
					xi0: 0.02,
					zeta0: 1.0,
					area: null,
					theta_source: null,
				} as any;
			} else if (type === "combustor") {
				data = {
					...data,
					method: "complete",
					area: 0.1,
					Nu_multiplier: 1.0,
				} as any;
			} else if (type === "momentum_chamber") {
				data = {
					...data,
					area: 0.1,
					Nu_multiplier: 1.0,
					f_multiplier: 1.0,
				} as any;
			} else if (type === "mpce_tee") {
				// F_C / F_branch / psi intentionally omitted so the backend
				// inherits geometry from connected channels (per schema default
				// `None = inherit`).
				data = {
					...data,
					flow_direction: "branch",
					theta_deg: 90,
				} as any;
			} else if (type === "vortex") {
				data = {
					...data,
					r_c: 0.02,
					r_out: 0.1,
					r_in: 0.0,
					omega_rpm: null,
					n: 2.0,
				} as any;
			}

			const newNode = {
				id,
				type,
				position,
				data,
			};

			// Adding a node is a structural edit: stored results solved a
			// different network and must not warm-start the next solve.
			setNodes(stripResults(nodes).concat(newNode));
		},
		[nodes, setNodes, screenToFlowPosition],
	);

	return (
		<main
			className="flex-grow h-full bg-stone-50 relative"
			onDragOver={onDragOver}
			onDrop={onDrop}
		>
			<SolverStatusBadge />
			<ReactFlow
				nodes={nodes}
				edges={edges}
				onNodesChange={onNodesChange}
				onEdgesChange={onEdgesChange}
				onConnect={onConnect}
				nodeTypes={nodeTypes}
				edgeTypes={edgeTypes}
				fitView
				fitViewOptions={{ padding: 0.5, maxZoom: 1.0 }}
			>
				<Background />
				<Controls />
				<MiniMap />
			</ReactFlow>
		</main>
	);
};

export default NetworkCanvas;
