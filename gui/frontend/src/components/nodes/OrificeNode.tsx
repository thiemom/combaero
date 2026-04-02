import { ChevronRight } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const OrificeNode = ({ data, selected }: NodeProps) => {
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
		gap: "4px",
		minHeight: isVertical ? 120 : 48,
		minWidth: isVertical ? 48 : 120,
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
			className={`px-3 py-1 shadow-sm rounded bg-stone-50 border-2 transition-all ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-green-400"
						: "border-stone-300"
			}`}
			style={containerStyle}
		>
			{/* Icon: Rotates with node */}
			<div
				className="flex items-center justify-center p-1 bg-white rounded border border-stone-200 pointer-events-none"
				style={{ transform: `rotate(${rotation}deg)` }}
			>
				<ChevronRight size={14} className="text-orange-400" />
			</div>

			{/* Text: Stays upright or vertical B-T */}
			<div
				className="flex flex-col items-center pointer-events-none"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold text-gray-400 uppercase leading-none">
					{data.label ? data.label : "Orifice"}
				</div>
				<div className="text-xs font-semibold whitespace-nowrap">
					Cd: {data.Cd}
				</div>
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

			{/* Thermal Handles */}
			<Handle
				type="target"
				position={getHandlePosition(Position.Top, rotation)}
				id="thermal-target"
				style={{ background: "#ff9800" }}
				className="w-2 h-2 border-none"
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Bottom, rotation)}
				id="thermal-source"
				style={{ background: "#ff9800" }}
				className="w-2 h-2 border-none"
			/>
		</div>
	);
};

export default OrificeNode;
