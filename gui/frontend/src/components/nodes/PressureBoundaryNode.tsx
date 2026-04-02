import { ArrowRightToLine } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const PressureBoundaryNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isVertical = rotation === 90 || rotation === 270;

	// Inner shape dimensions
	const innerWidth = isVertical ? 80 : 128;
	const innerHeight = isVertical ? 128 : 80;

	// Orientation logic: Vertical text stays Bottom-to-Top
	const textRotation = isVertical ? 270 : 0;

	return (
		<div className="relative flex items-center justify-center w-[128px] h-[128px]">
			<div
				className={`relative shadow-md rounded border-2 bg-stone-50 transition-all flex items-center justify-center ${
					selected ? "border-blue-500 shadow-blue-100" : "border-stone-400"
				}`}
				style={{
					width: innerWidth,
					height: innerHeight,
					gap: "8px",
					flexDirection:
						rotation === 90
							? ("column" as const)
							: rotation === 180
								? ("row-reverse" as const)
								: rotation === 270
									? ("column-reverse" as const)
									: ("row" as const),
				}}
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
					{data.label && (
						<span className="text-[10px] font-bold text-gray-500 uppercase">
							{data.label}
						</span>
					)}
					<div className="text-[10px] text-gray-500 font-mono whitespace-nowrap">
						{(data.P_total / 1e5).toFixed(2)} bar
					</div>
				</div>

				{/* Flow Handles */}
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
		</div>
	);
};

export default PressureBoundaryNode;
