import { ArrowRightToLine } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const MassBoundaryNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
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
		minWidth: isVertical ? 80 : 144,
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
			className={`px-4 py-2 shadow-md rounded border-2 bg-stone-50 transition-all ${
				selected ? "border-blue-500 shadow-blue-100" : "border-stone-400"
			}`}
			style={containerStyle}
		>
			{/* Icon: Rotates with node */}
			<div
				className="flex items-center justify-center pointer-events-none"
				style={{ transform: `rotate(${rotation}deg)` }}
			>
				<ArrowRightToLine size={24} className="text-blue-600" />
			</div>

			{/* Text: Stays upright or vertical B-T */}
			<div
				className="flex flex-col items-center pointer-events-none"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-xs font-bold text-center">
					{data.label ? data.label : "Mass Flow"}
				</div>
				<div className="text-[10px] text-gray-500 font-mono whitespace-nowrap">
					{data.m_dot} kg/s
				</div>
			</div>

			<Handle
				type="source"
				position={getHandlePosition(Position.Right, rotation)}
				id="flow-source"
				className="w-2 h-2 !bg-blue-600"
			/>
		</div>
	);
};

export default MassBoundaryNode;
