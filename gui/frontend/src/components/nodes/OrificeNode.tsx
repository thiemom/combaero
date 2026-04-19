import { ChevronRight } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";

const OrificeNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	return (
		<div
			className={`shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 px-3 py-1 overflow-hidden ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-green-400"
						: "border-stone-300"
			}`}
			style={{
				width: 120,
				height: 48,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center p-1 bg-white rounded border border-stone-200 shrink-0">
				<ChevronRight size={12} className="text-orange-400" />
			</div>

			<div
				className="flex flex-col items-start flex-1 min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap">
					{data.label ? data.label : "Orifice"}
				</div>
				<div className="text-[9px] font-mono text-gray-500 whitespace-nowrap">
					Cd: {(data.result?.Cd ?? data.Cd ?? 0.6).toFixed(3)}
				</div>
			</div>

			<Handle type="target" position={Position.Left} id="flow-target" />
			<Handle type="source" position={Position.Right} id="flow-source" />
			<Handle type="target" position={Position.Top} id="thermal-target" />
			<Handle type="source" position={Position.Bottom} id="thermal-source" />
		</div>
	);
};

export default OrificeNode;
