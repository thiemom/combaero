import { GitMerge } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { handleStyle, rotPos } from "../../utils/nodeUtils";
import { NodeDiagRows } from "../NodeDiagRows";

// Port colours: blue = input arm, green = output arm (matches handle triangle colours).
const portColor = (isInput: boolean) => (isInput ? "#3b82f6" : "#22c55e");

const MPCETeeNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;
	const isMerge = (data.flow_direction ?? "branch") === "merge";

	// Merge:  S(left, in) + B(bottom, in)  → C(right, out)
	// Branch: C(left, in) → S(right, out) + B(bottom, out)
	const straightBase = isMerge ? Position.Left : Position.Right;
	const commonBase = isMerge ? Position.Right : Position.Left;
	const branchBase = Position.Bottom;
	const straightPos = rotPos(straightBase, rotation);
	const commonPos = rotPos(commonBase, rotation);
	const branchPos = rotPos(branchBase, rotation);

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation and flow_direction both trigger handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, data.flow_direction, updateNodeInternals]);

	return (
		<div
			className={`shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 px-3 py-1 ${
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
			<div
				className={`flex items-center justify-center p-1 rounded border shrink-0 ${
					isMerge
						? "bg-indigo-50 border-indigo-200 text-indigo-500"
						: "bg-rose-50 border-rose-200 text-rose-500"
				}`}
			>
				<GitMerge size={16} />
			</div>

			<div
				className="flex flex-col items-start flex-1 min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap">
					{data.label ? data.label : "MPCE Tee"}
				</div>
				<div
					className={`text-[9px] font-mono whitespace-nowrap ${isMerge ? "text-indigo-500" : "text-rose-500"}`}
				>
					{isMerge ? "merge" : "branch"}{" "}
					{Math.round(Math.abs(data.theta_deg ?? 90))}&deg;
				</div>
				{isSolved && <NodeDiagRows result={data.result} maxRows={1} />}
			</div>

			{/* Port labels: S/C swap left/right per flow_direction; B stays at bottom. */}
			<div
				className={`absolute ${isMerge ? "-left-3" : "-right-3"} top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5`}
				style={{
					color: portColor(isMerge),
					transform: `rotate(${-rotation}deg)`,
				}}
			>
				S
			</div>
			<div
				className={`absolute ${isMerge ? "-right-3" : "-left-3"} top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5`}
				style={{
					color: portColor(!isMerge),
					transform: `rotate(${-rotation}deg)`,
				}}
			>
				C
			</div>
			<div
				className="absolute -bottom-3 left-1/2 -translate-x-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5"
				style={{
					color: portColor(isMerge),
					transform: `rotate(${-rotation}deg)`,
				}}
			>
				B
			</div>

			{/* Straight arm — left for merge, right for branch */}
			<Handle
				type="target"
				position={straightPos}
				style={handleStyle(straightBase, rotation)}
				id="port-straight-target"
			/>
			<Handle
				type="source"
				position={straightPos}
				style={handleStyle(straightBase, rotation)}
				id="port-straight-source"
			/>

			{/* Common arm — right for merge, left for branch */}
			<Handle
				type="target"
				position={commonPos}
				style={handleStyle(commonBase, rotation)}
				id="port-common-target"
			/>
			<Handle
				type="source"
				position={commonPos}
				style={handleStyle(commonBase, rotation)}
				id="port-common-source"
			/>

			{/* Branch arm — bottom at 0°, rotates with node */}
			<Handle
				type="target"
				position={branchPos}
				style={handleStyle(branchBase, rotation)}
				id="port-branch-target"
			/>
			<Handle
				type="source"
				position={branchPos}
				style={handleStyle(branchBase, rotation)}
				id="port-branch-source"
			/>
		</div>
	);
};

export default MPCETeeNode;
