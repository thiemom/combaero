import { ArrowRightToLine } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { getHandlePosition, getWrapperSize } from "../../utils/nodeUtils";

const NODE_W = 128;
const NODE_H = 56;
const WRAPPER = getWrapperSize(NODE_W, NODE_H);

const PressureBoundaryNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const updateNodeInternals = useUpdateNodeInternals();

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
			<div
				className={`shadow-md rounded border-2 bg-stone-50 flex items-center gap-2 px-4 py-2 ${
					selected ? "border-blue-500 shadow-blue-100" : "border-stone-400"
				}`}
				style={{
					width: NODE_W,
					height: NODE_H,
					transform: `rotate(${rotation}deg)`,
					transformOrigin: "center center",
				}}
			>
				<ArrowRightToLine size={20} className="text-blue-600 shrink-0" />

				<div
					className="flex flex-col items-start"
					style={{ transform: `rotate(${-rotation}deg)` }}
				>
					{data.label && (
						<span className="text-[10px] font-bold text-gray-500 uppercase whitespace-nowrap">
							{data.label}
						</span>
					)}
					<div className="text-[10px] text-gray-500 font-mono whitespace-nowrap">
						{(data.P_total / 1e5).toFixed(2)} bar
					</div>
				</div>
			</div>

			<Handle
				type="source"
				position={getHandlePosition(Position.Right, rotation)}
				id="flow-source"
			/>
			<Handle
				type="target"
				position={getHandlePosition(Position.Left, rotation)}
				id="flow-target"
			/>
		</div>
	);
};

export default PressureBoundaryNode;
