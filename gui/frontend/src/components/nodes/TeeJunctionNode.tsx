import { Split } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { HANDLE_CSS, rotPos } from "../../utils/nodeUtils";
import { NodeDiagRows } from "../NodeDiagRows";

// Port colours: blue = input arm, green = output arm (matches handle triangle colours).
const portColor = (isInput: boolean) => (isInput ? "#3b82f6" : "#22c55e");

const TeeJunctionNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;
	const isMerging = (data.tee_type ?? "merging") === "merging";
	const isExtrapolated =
		isSolved && (data.result?.correlation_extrapolated ?? 0) > 0;

	// Merging:  S(left, in) + B(bottom, in)  → C(right, out)
	// Branching: C(left, in) → S(right, out) + B(bottom, out)
	// Swap the left/right positions of S and C when the type changes.
	const straightBase = isMerging ? Position.Left : Position.Right;
	const commonBase = isMerging ? Position.Right : Position.Left;
	const branchBase = Position.Bottom;
	const straightPos = rotPos(straightBase, rotation);
	const commonPos = rotPos(commonBase, rotation);
	const branchPos = rotPos(branchBase, rotation);

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation and tee_type both trigger handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, data.tee_type, updateNodeInternals]);

	return (
		<div
			className={`shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 px-3 py-1 ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isExtrapolated
						? "border-amber-400"
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
					isMerging
						? "bg-blue-50 border-blue-200 text-blue-500"
						: "bg-orange-50 border-orange-200 text-orange-500"
				}`}
			>
				<Split size={16} />
			</div>

			<div
				className="flex flex-col items-start flex-1 min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap">
					{data.label ? data.label : "Tee"}
				</div>
				<div
					className={`text-[9px] font-mono whitespace-nowrap ${isMerging ? "text-blue-500" : "text-orange-500"}`}
				>
					{isMerging ? "merging" : "branching"}{" "}
					{Math.round(Math.abs(data.theta_deg ?? 90))}&deg;
				</div>
				{isSolved && <NodeDiagRows result={data.result} maxRows={1} />}
				{isExtrapolated && (
					<div className="text-[8px] font-bold text-amber-600 uppercase leading-none mt-0.5 whitespace-nowrap">
						extrapolated
					</div>
				)}
			</div>

			{/* Port labels sit just outside the visual bounding box.
			    S and C swap left/right when tee_type changes; B stays at bottom. */}
			<div
				className={`absolute ${isMerging ? "-left-3" : "-right-3"} top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5`}
				style={{
					color: portColor(isMerging),
					transform: `rotate(${-rotation}deg)`,
				}}
			>
				S
			</div>
			<div
				className={`absolute ${isMerging ? "-right-3" : "-left-3"} top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5`}
				style={{
					color: portColor(!isMerging),
					transform: `rotate(${-rotation}deg)`,
				}}
			>
				C
			</div>
			<div
				className="absolute -bottom-3 left-1/2 -translate-x-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5"
				style={{
					color: portColor(isMerging),
					transform: `rotate(${-rotation}deg)`,
				}}
			>
				B
			</div>

			{/* Straight arm — left for merging, right for branching */}
			<Handle
				type="target"
				position={straightPos}
				style={HANDLE_CSS[straightBase]}
				id="port-straight-target"
			/>
			<Handle
				type="source"
				position={straightPos}
				style={HANDLE_CSS[straightBase]}
				id="port-straight-source"
			/>

			{/* Common arm — right for merging, left for branching */}
			<Handle
				type="target"
				position={commonPos}
				style={HANDLE_CSS[commonBase]}
				id="port-common-target"
			/>
			<Handle
				type="source"
				position={commonPos}
				style={HANDLE_CSS[commonBase]}
				id="port-common-source"
			/>

			{/* Branch arm — bottom at 0°, rotates with node */}
			<Handle
				type="target"
				position={branchPos}
				style={HANDLE_CSS[branchBase]}
				id="port-branch-target"
			/>
			<Handle
				type="source"
				position={branchPos}
				style={HANDLE_CSS[branchBase]}
				id="port-branch-source"
			/>
		</div>
	);
};

export default TeeJunctionNode;
