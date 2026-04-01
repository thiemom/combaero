import { ArrowRightToLine } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const PressureBoundaryNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isVertical = rotation === 90 || rotation === 270;

	// Rotation 0 -> world 0 (H L-R)
	// Rotation 90 -> world 270 (V B-T)
	// Rotation 180 -> world 360/0 (H L-R)
	// Rotation 270 -> world 270 (V B-T)
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	return (
		<div
			className={`w-32 px-4 py-2 shadow-md rounded border-2 bg-stone-50 transition-colors ${
				selected ? "border-blue-500 shadow-blue-100" : "border-stone-400"
			}`}
			style={{
				display: "flex",
				flexDirection: isVertical ? "row" : "column",
				alignItems: "center",
				justifyContent: "space-between",
				minHeight: isVertical ? 32 : data.label ? 64 : 48,
				minWidth: isVertical ? (data.label ? 64 : 48) : data.label ? 128 : 96,
			}}
		>
			<Handle
				type="target"
				position={getHandlePosition(Position.Left, rotation)}
				style={{ background: "#555" }}
			/>

			<div
				className="flex items-center justify-center pointer-events-none"
				style={{
					transform: `rotate(${rotation}deg)`,
				}}
			>
				<div
					className="flex flex-col items-center"
					style={{ transform: `rotate(${textRotation}deg)` }}
				>
					<ArrowRightToLine size={24} className="text-blue-600 mb-1" />
					{data.label && (
						<span className="text-[10px] font-bold text-gray-500 uppercase">
							{data.label}
						</span>
					)}
				</div>
			</div>
		</div>
	);
};

export default PressureBoundaryNode;
