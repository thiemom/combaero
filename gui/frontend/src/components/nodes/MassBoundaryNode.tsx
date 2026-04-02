import { ArrowRightToLine } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";

const MassBoundaryNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	return (
		<div
			className={`shadow-md rounded border-2 bg-stone-50 flex items-center gap-2 px-4 py-2 ${
				selected ? "border-blue-500 shadow-blue-100" : "border-stone-400"
			}`}
			style={{
				width: 144,
				height: 56,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<ArrowRightToLine size={20} className="text-blue-600 shrink-0" />

			<div
				className="flex flex-col items-start"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-xs font-bold whitespace-nowrap">
					{data.label ? data.label : "Mass Flow"}
				</div>
				<div className="text-[10px] text-gray-500 font-mono whitespace-nowrap">
					{data.m_dot} kg/s
				</div>
			</div>

			<Handle type="source" position={Position.Right} id="flow-source" />
		</div>
	);
};

export default MassBoundaryNode;
