import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { NodeDiagRows } from "../NodeDiagRows";

const TeeIcon = () => (
	<svg
		width="20"
		height="20"
		viewBox="0 0 20 20"
		fill="none"
		stroke="currentColor"
		strokeWidth="2"
		strokeLinecap="round"
		aria-hidden="true"
	>
		<title>Tee junction</title>
		<line x1="1" y1="8" x2="19" y2="8" />
		<line x1="10" y1="8" x2="10" y2="19" />
		<circle cx="10" cy="8" r="1.5" fill="currentColor" stroke="none" />
	</svg>
);

const TeeJunctionNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;
	const isMerging = (data.tee_type ?? "merging") === "merging";
	const isExtrapolated =
		isSolved && (data.result?.correlation_extrapolated ?? 0) > 0;

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	return (
		<div
			className={`shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 px-3 py-1 overflow-hidden ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isExtrapolated
						? "border-amber-400"
						: isSolved
							? "border-green-400"
							: "border-stone-300"
			}`}
			style={{
				width: 128,
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
				<TeeIcon />
			</div>

			<div
				className="flex flex-col items-start flex-1 min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap">
					{data.label ? data.label : "Tee Junction"}
				</div>
				<div
					className={`text-[9px] font-mono whitespace-nowrap ${isMerging ? "text-blue-500" : "text-orange-500"}`}
				>
					{isMerging ? "merging" : "branching"}
				</div>
				{isSolved && <NodeDiagRows result={data.result} maxRows={1} />}
				{isExtrapolated && (
					<div className="text-[8px] font-bold text-amber-600 uppercase leading-none mt-0.5">
						EXT
					</div>
				)}
			</div>

			{/* Straight arm -- left */}
			<Handle
				type="target"
				position={Position.Left}
				id="port-straight-target"
			/>
			<Handle
				type="source"
				position={Position.Left}
				id="port-straight-source"
			/>

			{/* Common arm -- right */}
			<Handle type="target" position={Position.Right} id="port-common-target" />
			<Handle type="source" position={Position.Right} id="port-common-source" />

			{/* Branch arm -- bottom */}
			<Handle
				type="target"
				position={Position.Bottom}
				id="port-branch-target"
			/>
			<Handle
				type="source"
				position={Position.Bottom}
				id="port-branch-source"
			/>
		</div>
	);
};

export default TeeJunctionNode;
