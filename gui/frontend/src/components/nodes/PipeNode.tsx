import { Settings2 } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";

const PipeNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();

	// Text counter-rotation: L→R for horizontal nodes, bottom→top for vertical
	// Node 0° → text 0° (L→R), Node 90° → text 180° (net 270° = B→T)
	// Node 180° → text 180° (net 0° = L→R), Node 270° → text 0° (net 270° = B→T)
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	return (
		<div
			className={`shadow-md rounded border-2 bg-white flex items-center gap-2 px-3 py-1 ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-green-400"
						: "border-stone-300"
			}`}
			style={{
				width: 140,
				height: 56,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center p-1 bg-stone-100 rounded shrink-0">
				<Settings2 size={16} className="text-stone-600" />
			</div>

			<div
				className="flex flex-col items-start"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<span className="text-[10px] font-bold text-gray-400 uppercase leading-none whitespace-nowrap">
					{data.label ? data.label : "Pipe"}
				</span>
				<span className="text-[10px] font-mono text-gray-500 whitespace-nowrap">
					L:{data.L}m D:{data.D}m
				</span>
			</div>

			{/* Handles: FIXED positions. CSS rotation moves them visually. */}
			<Handle type="target" position={Position.Left} id="flow-target" />
			<Handle type="source" position={Position.Right} id="flow-source" />
			<Handle type="target" position={Position.Top} id="thermal-target" />
			<Handle type="source" position={Position.Bottom} id="thermal-source" />
		</div>
	);
};

export default PipeNode;
