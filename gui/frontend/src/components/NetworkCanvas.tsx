import ReactFlow, { Background, Controls, MiniMap } from "reactflow";
import "reactflow/dist/style.css";
import useStore from "../store/useStore";
import MassBoundaryNode from "./nodes/MassBoundaryNode";
import OrificeNode from "./nodes/OrificeNode";
import PipeNode from "./nodes/PipeNode";
import PlenumNode from "./nodes/PlenumNode";
import PressureBoundaryNode from "./nodes/PressureBoundaryNode";

const NODE_TYPES = Object.freeze({
	plenum: PlenumNode,
	mass_boundary: MassBoundaryNode,
	pressure_boundary: PressureBoundaryNode,
	pipe: PipeNode,
	orifice: OrificeNode,
});

const NetworkCanvas = () => {
	const { nodes, edges, onNodesChange, onEdgesChange, onConnect } = useStore();

	return (
		<div className="flex-grow h-full bg-stone-50">
			<ReactFlow
				nodes={nodes}
				edges={edges}
				onNodesChange={onNodesChange}
				onEdgesChange={onEdgesChange}
				onConnect={onConnect}
				nodeTypes={NODE_TYPES}
				fitView
			>
				<Background />
				<Controls />
				<MiniMap />
			</ReactFlow>
		</div>
	);
};

export default NetworkCanvas;
