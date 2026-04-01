import { Settings2 } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const PipeNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const isVertical = rotation === 90 || rotation === 270;

	// Orientation logic: keeps text readable
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	return (
		<div
			className={`px-4 py-2 shadow-md rounded border-2 bg-white transition-colors flex items-center justify-center ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-green-400"
						: "border-stone-300"
			}`}
			style={{
				flexDirection: isVertical ? "row" : "column",
				minHeight: isVertical ? 32 : data.label ? 64 : 48,
				minWidth: isVertical ? (data.label ? 64 : 48) : data.label ? 128 : 96,
			}}
		>
			<div
				className="flex items-center gap-3 pointer-events-none"
				style={{
					transform: `rotate(${rotation}deg)`,
				}}
			>
				<div
					className="flex flex-col items-center"
					style={{ transform: `rotate(${textRotation}deg)` }}
				>
					<div className="flex items-center gap-2">
						<div className="p-1 bg-stone-100 rounded">
							<Settings2 size={16} className="text-stone-600" />
						</div>
						<div className="flex flex-col">
							<span className="text-[10px] font-bold text-gray-400 uppercase leading-tight">
								{data.label ? data.label : "Pipe"}
							</span>
							<span className="text-xs font-bold">
								L: {data.L}m | D: {data.D}m
							</span>
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

export default PipeNode;
