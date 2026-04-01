import { ArrowRightToLine } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const PressureBoundaryNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isVertical = rotation === 90 || rotation === 270;

	// Orientation logic for text only
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	return (
		<div
			className={`w-32 px-4 py-2 shadow-md rounded border-2 bg-stone-50 transition-colors flex items-center justify-center ${
				selected ? "border-blue-500 shadow-blue-100" : "border-stone-400"
			}`}
			style={{
				flexDirection: isVertical ? "row" : "column",
				minHeight: isVertical ? 32 : 80,
				minWidth: isVertical ? 80 : 128,
			}}
		>
			<div
				className="flex items-center justify-center pointer-events-none"
				style={{
					transform: `rotate(${rotation}deg)`,
				}}
			>
				<ArrowRightToLine size={24} className="text-blue-600" />
			</div>

			<div
				className="flex flex-col items-center ml-2"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				{data.label && (
					<span className="text-[10px] font-bold text-gray-500 uppercase">
						{data.label}
					</span>
				)}
				<div className="text-[10px] text-gray-500 font-mono">
					{(data.P_total / 1e5).toFixed(2)} bar
				</div>
			</div>

			{/* Use Source for Inlet flow, Target for Outlet sink - Having One Source is common for Boundaries */}
			<Handle
				type="source"
				position={getHandlePosition(Position.Right, rotation)}
				id="flow-source"
				className="w-2 h-2 !bg-blue-600"
			/>
			<Handle
				type="target"
				position={getHandlePosition(Position.Left, rotation)}
				id="flow-target"
				className="w-2 h-2 !bg-blue-300"
			/>
		</div>
	);
};

export default PressureBoundaryNode;
