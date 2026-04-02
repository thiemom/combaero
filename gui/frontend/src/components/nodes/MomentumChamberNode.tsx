import { CircleDot } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";

const MomentumChamberNode = ({ id, data, selected }: NodeProps) => {
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
			className={`shadow-md rounded border-2 bg-white flex items-center justify-center gap-2 px-4 py-2 ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-blue-500"
						: "border-stone-300"
			}`}
			style={{
				width: 140,
				height: 80,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center rounded-full w-10 h-10 bg-blue-50 p-2 shrink-0">
				<CircleDot size={20} className="text-blue-600" />
			</div>

			<div
				className="flex flex-col items-start"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-sm font-bold whitespace-nowrap">
					{data.label ? data.label : "Mom. Chamber"}
				</div>
				{isSolved && data.result.state && (
					<div className="text-[10px] text-gray-500 font-mono whitespace-nowrap">
						{data.result.state.T.toFixed(1)} K |{" "}
						{(data.result.state.P / 1e5).toFixed(2)} bar
					</div>
				)}
			</div>

			<Handle type="target" position={Position.Left} id="flow-target" />
			<Handle type="source" position={Position.Right} id="flow-source" />
		</div>
	);
};

export default MomentumChamberNode;
