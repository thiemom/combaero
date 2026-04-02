import { Settings2 } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { getHandlePosition, getWrapperSize } from "../../utils/nodeUtils";

// Base dimensions of the Pipe node (unrotated)
const NODE_W = 140;
const NODE_H = 56;
const WRAPPER = getWrapperSize(NODE_W, NODE_H);

const PipeNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();

	// Re-measure handles whenever rotation changes
	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	return (
		<div
			className="relative"
			style={{
				width: WRAPPER,
				height: WRAPPER,
				display: "flex",
				alignItems: "center",
				justifyContent: "center",
			}}
		>
			{/* Inner visible box — rotates as a unit */}
			<div
				className={`shadow-md rounded border-2 bg-white flex items-center gap-2 px-3 py-1 ${
					selected
						? "border-blue-500 shadow-blue-100"
						: isSolved
							? "border-green-400"
							: "border-stone-300"
				}`}
				style={{
					width: NODE_W,
					height: NODE_H,
					transform: `rotate(${rotation}deg)`,
					transformOrigin: "center center",
				}}
			>
				{/* Icon */}
				<div className="flex items-center justify-center p-1 bg-stone-100 rounded shrink-0">
					<Settings2 size={16} className="text-stone-600" />
				</div>

				{/* Text: counter-rotate to stay horizontal */}
				<div
					className="flex flex-col items-start"
					style={{ transform: `rotate(${-rotation}deg)` }}
				>
					<span className="text-[10px] font-bold text-gray-400 uppercase leading-none whitespace-nowrap">
						{data.label ? data.label : "Pipe"}
					</span>
					<span className="text-[10px] font-mono text-gray-500 whitespace-nowrap">
						L:{data.L}m D:{data.D}m
					</span>
				</div>
			</div>

			{/* Handles — on OUTER wrapper, unaffected by inner rotation */}
			<Handle
				type="target"
				position={getHandlePosition(Position.Left, rotation)}
				id="flow-target"
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Right, rotation)}
				id="flow-source"
			/>
			<Handle
				type="target"
				position={getHandlePosition(Position.Top, rotation)}
				id="thermal-target"
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Bottom, rotation)}
				id="thermal-source"
			/>
		</div>
	);
};

export default PipeNode;
