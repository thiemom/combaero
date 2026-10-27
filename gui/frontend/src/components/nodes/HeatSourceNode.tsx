import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import {
	formatHeat,
	HEAT_SOURCE_HANDLE,
	heatSourceSummary,
} from "../../utils/heatSource";
import { handleStyle, rotPos } from "../../utils/nodeUtils";
import NodeIcon from "../NodeIcon";

/**
 * A heat source/sink (#792): an EnergyBoundary on the plenum, chamber or
 * combustor its one edge points at. The card shows what is set and, once
 * solved, the heat applied (the solved Q in an outlet-Tt mode).
 */
const HeatSourceNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;
	const result = data.result;

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	return (
		<div
			className={`shadow-md rounded border-2 bg-rose-50 flex items-center gap-2 px-3 py-2 ${
				selected
					? "border-blue-500 shadow-blue-100"
					: result
						? "border-rose-400"
						: "border-rose-200"
			}`}
			style={{
				width: 150,
				height: 60,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<NodeIcon type="heat_source" size={20} className="shrink-0" />
			<div
				className="flex flex-col items-start min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<span className="text-[10px] font-bold text-gray-600 uppercase whitespace-nowrap truncate max-w-full">
					{data.label || "Heat Source"}
				</span>
				<span className="text-[10px] text-stone-500 font-mono whitespace-nowrap truncate max-w-full">
					{heatSourceSummary(data)}
				</span>
				{result && (
					<span className="text-[10px] text-rose-700 font-mono font-bold whitespace-nowrap truncate max-w-full">
						Q = {formatHeat(result.Q)}
					</span>
				)}
			</div>
			<Handle
				type="source"
				position={rotPos(Position.Right, rotation)}
				style={handleStyle(Position.Right, rotation)}
				id={HEAT_SOURCE_HANDLE}
			/>
		</div>
	);
};

export default HeatSourceNode;
