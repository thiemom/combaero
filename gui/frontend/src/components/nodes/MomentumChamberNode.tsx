import { CircleDot } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { NodeDiagRows } from "../NodeDiagRows";

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
			className={`shadow-md rounded border-2 bg-white flex items-center gap-2 px-3 py-2 ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-blue-500"
						: "border-stone-300"
			}`}
			style={{
				width: 160,
				height: 72,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center rounded-full w-8 h-8 bg-blue-50 p-1.5 shrink-0">
				<CircleDot size={16} className="text-blue-600" />
			</div>

			<div
				className="flex flex-col items-center flex-1 min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-xs font-bold whitespace-nowrap">
					{data.label ? data.label : "Mom. Chamber"}
				</div>
				{isSolved && <NodeDiagRows result={data.result} />}
			</div>

			<div
				className="absolute -left-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-blue-500"
				style={{ transform: `rotate(${-rotation}deg)` }}
			>
				i
			</div>
			<div
				className="absolute -right-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-green-500"
				style={{ transform: `rotate(${-rotation}deg)` }}
			>
				o
			</div>
			<div
				className="absolute -top-3 left-1/2 -translate-x-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-orange-400"
				style={{ transform: `rotate(${-rotation}deg)` }}
			>
				q
			</div>
			<div
				className="absolute -bottom-3 left-1/2 -translate-x-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-orange-400"
				style={{ transform: `rotate(${-rotation}deg)` }}
			>
				q
			</div>

			{/* Flow handles */}
			<Handle type="target" position={Position.Left} id="flow-target" />
			<Handle type="source" position={Position.Right} id="flow-source" />
			<Handle type="target" position={Position.Top} id="thermal-target" />
			<Handle type="source" position={Position.Bottom} id="thermal-source" />
		</div>
	);
};

export default MomentumChamberNode;
