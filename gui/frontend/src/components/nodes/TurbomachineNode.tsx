import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { handleStyle, rotPos } from "../../utils/nodeUtils";
import NodeIcon from "../NodeIcon";

/**
 * A compressor or turbine at a fixed total-to-total pressure ratio (#893).
 * The flow is set by the network; the card shows the pressure ratio and the
 * efficiency, and the shaft power once a solve has run (power > 0 for both:
 * absorbed by a compressor, delivered by a turbine). No surface: no thermal
 * handles.
 */
const TurbomachineNode = ({ id, type, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;
	const kind = type === "turbine" ? "turbine" : "compressor";

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	const power = data.result?.power;
	const eta = data.efficiency ?? (kind === "turbine" ? 0.9 : 0.88);
	const etaName = data.efficiency_type === "isentropic" ? "ηs" : "ηp";
	return (
		<div
			className={`shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 px-3 py-1 ${
				selected
					? "border-blue-500 shadow-blue-100"
					: data.result
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
			<div className="flex items-center justify-center p-1 bg-white rounded border border-stone-200 shrink-0">
				<NodeIcon type={kind} size={14} />
			</div>
			<div
				className="flex flex-col items-start flex-1 min-w-0 overflow-hidden"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap truncate max-w-full">
					{data.label
						? data.label
						: kind === "turbine"
							? "Turbine"
							: "Compressor"}
				</div>
				<div className="text-[9px] font-mono text-gray-500 whitespace-nowrap truncate max-w-full">
					{`PR ${Number(data.pressure_ratio ?? 4).toPrecision(3)} ${etaName} ${Number(eta).toFixed(3)}`}
				</div>
				{typeof power === "number" && (
					<div className="text-[9px] font-mono text-gray-500 whitespace-nowrap truncate max-w-full">
						{`${(power / 1000).toPrecision(4)} kW`}
					</div>
				)}
			</div>
			<Handle
				type="target"
				position={rotPos(Position.Left, rotation)}
				style={handleStyle(Position.Left, rotation)}
				id="flow-target"
			/>
			<Handle
				type="source"
				position={rotPos(Position.Right, rotation)}
				style={handleStyle(Position.Right, rotation)}
				id="flow-source"
			/>
		</div>
	);
};

export default TurbomachineNode;
