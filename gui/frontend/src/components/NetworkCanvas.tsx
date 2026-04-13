import ReactFlow, {
	Background,
	Controls,
	MiniMap,
	useReactFlow,
} from "reactflow";
import "reactflow/dist/style.css";
import { useCallback } from "react";
import useStore from "../store/useStore";
import CombustorNode from "./nodes/CombustorNode.tsx";
import LosslessNode from "./nodes/LosslessNode";
import MassBoundaryNode from "./nodes/MassBoundaryNode";
import MomentumChamberNode from "./nodes/MomentumChamberNode.tsx";
import OrificeNode from "./nodes/OrificeNode";
import PipeNode from "./nodes/PipeNode";
import PlenumNode from "./nodes/PlenumNode";
import PressureBoundaryNode from "./nodes/PressureBoundaryNode";
import ProbeNode from "./nodes/ProbeNode";
import ThermalEdge from "./ThermalEdge";
import SolverStatusBadge from "./SolverStatusBadge";

const NODE_TYPES = Object.freeze({
	plenum: PlenumNode,
	mass_boundary: MassBoundaryNode,
	pressure_boundary: PressureBoundaryNode,
	pipe: PipeNode,
	orifice: OrificeNode,
	combustor: CombustorNode,
	momentum_chamber: MomentumChamberNode,
	lossless_connection: LosslessNode,
	probe: ProbeNode,
});

const EDGE_TYPES = Object.freeze({
	thermal: ThermalEdge,
});

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
			} else if (type === "pipe") {
				data = { ...data, L: 1.0, D: 0.1, roughness: 1e-5 } as any;
			} else if (type === "orifice") {
				data = { ...data, area: 0.01, Cd: 0.6 } as any;
			} else if (type === "combustor") {
				data = { ...data, method: "complete" } as any;
			} else if (type === "momentum_chamber") {
				data = { ...data, area: 0.1 } as any;
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
				nodeTypes={NODE_TYPES}
				edgeTypes={EDGE_TYPES}
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
