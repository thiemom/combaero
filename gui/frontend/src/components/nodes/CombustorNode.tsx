import { Flame } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const CombustorNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const isVertical = rotation === 90 || rotation === 270;

	// Inner shape dimensions
	const innerWidth = isVertical ? 80 : 140;
	const innerHeight = isVertical ? 140 : 80;

	// Orientation logic: Vertical text stays Bottom-to-Top
	const textRotation = isVertical ? 270 : 0;

	return (
		<div className="relative flex items-center justify-center w-[140px] h-[140px]">
			<div
				className={`relative shadow-md rounded border-2 bg-white transition-all flex items-center justify-center ${
					selected
						? "border-red-500 shadow-red-100"
						: isSolved
							? "border-red-500"
							: "border-red-300"
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
					className="flex items-center justify-center rounded-full w-10 h-10 bg-red-100 p-2 pointer-events-none"
					style={{ transform: `rotate(${rotation}deg)` }}
				>
					<Flame size={20} className="text-red-600" />
				</div>

				{/* Text: Stays upright or vertical B-T */}
				<div
					className="flex flex-col items-center pointer-events-none"
					style={{ transform: `rotate(${textRotation}deg)` }}
				>
					<div className="text-sm font-bold">
						{data.label ? data.label : "Combustor"}
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
		</div>
	);
};

export default CombustorNode;
