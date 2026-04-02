import { CircleDot } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const MomentumChamberNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const isVertical = rotation === 90 || rotation === 270;

	// Orientation logic: Vertical text stays Bottom-to-Top
	const textRotation = isVertical ? 270 : 0;

	// Layout sizing and flex direction
	const containerStyle = {
		display: "flex",
		alignItems: "center",
		justifyContent: "center",
		gap: "8px",
		minHeight: isVertical ? 120 : 80,
		minWidth: isVertical ? 80 : 140,
		flexDirection:
			rotation === 90
				? ("column" as const)
				: rotation === 180
					? ("row-reverse" as const)
					: rotation === 270
						? ("column-reverse" as const)
						: ("row" as const),
	};

	return (
		<div
			className={`px-4 py-2 shadow-md rounded border-2 bg-white transition-all ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-blue-500"
						: "border-stone-300"
			}`}
			style={containerStyle}
		>
			{/* Icon: Rotates with node */}
			<div
				className="flex items-center justify-center rounded-full w-10 h-10 bg-blue-50 p-2 pointer-events-none"
				style={{ transform: `rotate(${rotation}deg)` }}
			>
				<CircleDot size={20} className="text-blue-600" />
			</div>

			{/* Text: Stays upright or vertical B-T */}
			<div
				className="flex flex-col items-center pointer-events-none"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-sm font-bold">
					{data.label ? data.label : "Mom. Chamber"}
				</div>
				{isSolved && data.result.state && (
					<div className="text-[10px] text-gray-500 font-mono whitespace-nowrap">
						{data.result.state.T.toFixed(1)} K |{" "}
						{(data.result.state.P / 1e5).toFixed(2)} bar
					</div>
				)}
			</div>

			{/* Flow Handles */}
			<Handle
				type="target"
				position={getHandlePosition(Position.Left, rotation)}
				id="flow-target"
				className="w-2 h-2 !bg-blue-500"
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Right, rotation)}
				id="flow-source"
				className="w-2 h-2 !bg-blue-500"
			/>
		</div>
	);
};

export default MomentumChamberNode;
