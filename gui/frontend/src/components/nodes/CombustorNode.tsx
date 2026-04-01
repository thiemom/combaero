import { Flame } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const CombustorNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const isVertical = rotation === 90 || rotation === 270;

	// Orientation logic
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	return (
		<div
			className={`px-4 py-2 shadow-md rounded border-2 bg-white transition-colors flex items-center justify-center ${
				selected
					? "border-red-500 shadow-red-100"
					: isSolved
						? "border-red-500"
						: "border-red-300"
			}`}
			style={{
				flexDirection: isVertical ? "row" : "column",
				minHeight: isVertical ? 64 : 80,
				minWidth: isVertical ? 80 : 140,
			}}
		>
			<div
				className="flex items-center gap-2 pointer-events-none"
				style={{
					transform: `rotate(${rotation}deg)`,
				}}
			>
				<div
					className="flex flex-col items-center"
					style={{ transform: `rotate(${textRotation}deg)` }}
				>
					<div className="flex items-center gap-3">
						<div className="rounded-full w-10 h-10 flex justify-center items-center bg-red-100 p-2">
							<Flame size={20} className="text-red-600" />
						</div>
						<div className="flex flex-col">
							<div className="text-sm font-bold">
								{data.label ? data.label : "Combustor"}
							</div>
							{isSolved && data.result.state && (
								<div className="text-[10px] text-gray-500 font-mono">
									{data.result.state.T.toFixed(1)} K |{" "}
									{(data.result.state.P / 1e5).toFixed(2)} bar
								</div>
							)}
						</div>
					</div>
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

export default CombustorNode;
