import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import useStore from "../../store/useStore";
import { handleStyle, rotPos } from "../../utils/nodeUtils";
import { FUEL_LINK_COLOR } from "../../utils/probeLink";
import NodeIcon from "../NodeIcon";

/**
 * A fuel controller (#892): supplies the fuel flow that puts the combustor
 * named by data.combustor_id at Tt_target. The flow is solved, not given;
 * the card shows it once a solve has run. The dashed link to the combustor
 * is derived from combustor_id (utils/probeLink.ts), like a probe's.
 */
const FuelBoundaryNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;
	const combustor = useStore((s) =>
		data.combustor_id
			? s.nodes.find((n) => n.id === data.combustor_id)
			: undefined,
	);

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	const m = data.result?.state?.m_dot;
	const branch = data.branch ?? "lean";
	return (
		<div
			className={`shadow-md rounded border-2 bg-stone-50 flex items-center gap-2 px-4 py-2 ${
				selected ? "border-blue-500 shadow-blue-100" : "border-stone-400"
			}`}
			style={{
				width: 150,
				height: 64,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<NodeIcon type="fuel_boundary" size={20} className="shrink-0" />

			<div
				className="flex flex-col items-start min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-xs font-bold whitespace-nowrap">
					{data.label ? data.label : "Fuel"}
				</div>
				<div className="text-[10px] text-gray-500 font-mono whitespace-nowrap">
					{`${data.Tt_target ?? 1500} K ${branch}`}
				</div>
				<div className="text-[10px] text-gray-500 font-mono whitespace-nowrap truncate max-w-full">
					{combustor
						? typeof m === "number"
							? `${m.toPrecision(4)} kg/s`
							: `sets ${combustor.data?.label || combustor.id}`
						: "no combustor"}
				</div>
			</div>

			<Handle
				type="source"
				position={rotPos(Position.Right, rotation)}
				style={handleStyle(Position.Right, rotation)}
				id="flow-source"
			/>
			{/* The link to the combustor it controls (utils/probeLink.ts): the
			    derived edge loops probe-out -> probe-in on this card. Drag from
			    probe-out onto a combustor to choose it. */}
			<Handle
				type="target"
				position={Position.Left}
				id="probe-in"
				style={{ opacity: 0, width: 6, height: 6 }}
			/>
			<Handle
				type="source"
				position={Position.Top}
				id="probe-out"
				title="Drag onto the combustor whose exit Tt this fuel sets"
				style={{
					width: 8,
					height: 8,
					background: "#fff",
					border: `1.5px dashed ${FUEL_LINK_COLOR}`,
				}}
			/>
		</div>
	);
};

export default FuelBoundaryNode;
