import { Activity } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { NodeDiagRows } from "../NodeDiagRows";

const DiscreteLossNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	const formulaByType: Record<string, string> = {
		constant_fraction: "ΔP = ξ·P",
		constant_head: "ΔP = ζ·½ρv²",
		linear_theta_fraction: "ξ = k·Θ + ξ₀",
		linear_theta_head: "ζ = k·Θ + ζ₀",
	};
	const formula = formulaByType[data.correlation_type] ?? "ΔP = ξ·P";

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	return (
		<div
			className={`shadow-sm rounded bg-white border-2 flex items-center gap-2 px-3 py-1 overflow-hidden ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-green-400"
						: "border-stone-300"
			}`}
			style={{
				width: 144,
				height: 58,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center p-1 bg-stone-100 rounded shrink-0">
				<Activity
					size={14}
					className="text-stone-600"
					style={{ transform: "rotate(90deg)" }}
				/>
			</div>

			<div
				className="flex flex-col items-start flex-1 min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap">
					{data.label ? data.label : "Discrete Loss"}
				</div>
				{isSolved ? (
					<NodeDiagRows result={data.result} maxRows={1} />
				) : (
					<div className="text-[9px] font-mono text-gray-500 whitespace-nowrap">
						{formula}
					</div>
				)}
			</div>

			<Handle type="target" position={Position.Left} id="flow-target" />
			<Handle type="source" position={Position.Right} id="flow-source" />
			<Handle type="target" position={Position.Top} id="thermal-target" />
			<Handle type="source" position={Position.Bottom} id="thermal-source" />
		</div>
	);
};

export default DiscreteLossNode;
