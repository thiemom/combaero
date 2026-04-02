import { Database } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { getHandlePosition, getWrapperSize } from "../../utils/nodeUtils";

const NODE_W = 140;
const NODE_H = 80;
const WRAPPER = getWrapperSize(NODE_W, NODE_H);

const PlenumNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
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
				className={`shadow-md rounded border-2 bg-stone-50 flex items-center justify-center gap-2 px-4 py-2 ${
					selected
						? "border-blue-500 shadow-blue-100"
						: isSolved
							? "border-amber-400"
							: "border-stone-300"
				}`}
				style={{
					width: NODE_W,
					height: NODE_H,
					transform: `rotate(${rotation}deg)`,
					transformOrigin: "center center",
				}}
			>
				<div className="flex items-center justify-center rounded-full w-10 h-10 bg-amber-100 p-2 shrink-0">
					<Database size={20} className="text-amber-600" />
				</div>

				<div
					className="flex flex-col items-start"
					style={{ transform: `rotate(${-rotation}deg)` }}
				>
					<div className="text-sm font-bold whitespace-nowrap">
						{data.label ? data.label : "Plenum"}
					</div>
					{isSolved && data.result.state && (
						<div className="text-[10px] text-gray-500 font-mono whitespace-nowrap">
							{data.result.state.T.toFixed(1)} K |{" "}
							{(data.result.state.P / 1e5).toFixed(2)} bar
						</div>
					)}
				</div>
			</div>

			{/* Flow Handles: L/T = target (inlet), R/B = source (outlet) */}
			<Handle
				type="target"
				position={getHandlePosition(Position.Left, rotation)}
				id="flow-left"
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Right, rotation)}
				id="flow-right"
			/>
			<Handle
				type="target"
				position={getHandlePosition(Position.Top, rotation)}
				id="flow-top"
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Bottom, rotation)}
				id="flow-bottom"
			/>

			{/* Thermal Handles */}
			<Handle
				type="target"
				position={getHandlePosition(Position.Top, rotation)}
				id="thermal-target"
				style={{ left: "25%" }}
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Bottom, rotation)}
				id="thermal-source"
				style={{ left: "25%" }}
			/>
		</div>
	);
};

export default PlenumNode;
