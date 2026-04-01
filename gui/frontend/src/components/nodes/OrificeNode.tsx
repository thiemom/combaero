import { ChevronRight } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const OrificeNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;

	// Orientation logic
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	return (
		<div
			className={`px-3 py-1 shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 transition-colors ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-green-400"
						: "border-stone-300"
			}`}
			style={{ minWidth: "120px" }}
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
					<div className="flex items-center gap-2">
						<div className="flex flex-col items-center justify-center p-1 bg-white rounded border border-stone-200">
							<ChevronRight size={14} className="text-orange-400" />
						</div>
						<div className="flex flex-col">
							<div className="text-[10px] font-bold text-gray-400 uppercase leading-none">
								{data.label ? data.label : "Orifice"}
							</div>
							<div className="text-xs font-semibold">Cd: {data.Cd}</div>
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

export default OrificeNode;
