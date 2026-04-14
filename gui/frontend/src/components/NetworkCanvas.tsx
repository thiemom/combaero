import ReactFlow, {
	Background,
	Controls,
	MiniMap,
	useReactFlow,
} from "reactflow";
import "reactflow/dist/style.css";
import { useCallback } from "react";
import useStore from "../store/useStore";
import { edgeTypes, nodeTypes } from "./flowTypes";
import SolverStatusBadge from "./SolverStatusBadge";

const NetworkCanvas = () => {
	const { nodes, edges, setNodes, onNodesChange, onEdgesChange, onConnect } =
		useStore();
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
					T_total: 300,
					composition: {
						source: "dry_air",
						mode: "mole",
						relative_humidity: 0.0,
					},
				} as any;
			} else if (type === "pressure_boundary") {
				data = {
					...data,
					P_total: 101325,
					T_total: 300,
					composition: {
						source: "dry_air",
						mode: "mole",
						relative_humidity: 0.0,
					},
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
			} else if (type === "combustor") {
				data = { ...data, method: "complete" } as any;
			} else if (type === "momentum_chamber") {
				data = {
					...data,
					area: 0.1,
					Nu_multiplier: 1.0,
					f_multiplier: 1.0,
				} as any;
			}

			const newNode = {
				id,
				type,
				position,
				data,
			};

			setNodes(nodes.concat(newNode));
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
