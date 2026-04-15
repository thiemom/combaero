import { Flame } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";

const CombustorNode = ({ id, data, selected }: NodeProps) => {
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
			className={`shadow-md rounded border-2 bg-white flex items-center gap-2 px-3 py-2 overflow-hidden ${
				selected
					? "border-red-500 shadow-red-100"
					: isSolved
						? "border-red-500"
						: "border-red-300"
			}`}
			style={{
				width: 160,
				height: 72,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center rounded-full w-8 h-8 bg-red-100 p-1.5 shrink-0">
				<Flame size={16} className="text-red-600" />
			</div>

			<div
				className="flex flex-col items-center flex-1 min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-xs font-bold whitespace-nowrap">
					{data.label ? data.label : "Combustor"}
				</div>
				{isSolved && data.result.state && (
					<>
						<div className="text-[9px] text-gray-500 font-mono whitespace-nowrap">
							{data.result.state.T.toFixed(0)} K
						</div>
						<div className="text-[9px] text-gray-500 font-mono whitespace-nowrap">
							{(data.result.state.P / 1e5).toFixed(2)} bar
						</div>
					</>
				)}
			</div>

			{/* Flow handles */}
			<Handle type="target" position={Position.Left} id="flow-target" />
			<Handle type="source" position={Position.Right} id="flow-source" />
			<Handle type="target" position={Position.Top} id="thermal-target" />
			<Handle type="source" position={Position.Bottom} id="thermal-source" />
		</div>
	);
};

export default CombustorNode;
