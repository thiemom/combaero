import { Database } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const PlenumNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;

	// Orientation logic
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	return (
		<div
			className={`px-4 py-2 shadow-md rounded border-2 bg-stone-50 transition-colors ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-amber-400"
						: "border-stone-300"
			}`}
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
						<div className="rounded-full w-10 h-10 flex justify-center items-center bg-amber-100 p-2">
							<Database size={20} className="text-amber-600" />
						</div>
						<div className="flex flex-col">
							<div className="text-sm font-bold">
								{data.label ? data.label : "Plenum"}
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

			{/* Flow Handles (Multiple) */}
			<Handle
				type="target"
				position={getHandlePosition(Position.Left, rotation)}
				id="flow-left"
				className="w-2 h-2 !bg-blue-500"
			/>
			<Handle
				type="target"
				position={getHandlePosition(Position.Right, rotation)}
				id="flow-right"
				className="w-2 h-2 !bg-blue-500"
			/>
			<Handle
				type="target"
				position={getHandlePosition(Position.Top, rotation)}
				id="flow-top"
				className="w-2 h-2 !bg-blue-500"
			/>
			<Handle
				type="target"
				position={getHandlePosition(Position.Bottom, rotation)}
				id="flow-bottom"
				className="w-2 h-2 !bg-blue-500"
			/>

			{/* Thermal Handles */}
			<Handle
				type="target"
				position={getHandlePosition(Position.Top, rotation)}
				id="thermal-target"
				style={{ background: "#ff9800", left: "20%" }}
				className="w-2 h-2 border-none"
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Bottom, rotation)}
				id="thermal-source"
				style={{ background: "#ff9800", left: "20%" }}
				className="w-2 h-2 border-none"
			/>
		</div>
	);
};

export default PlenumNode;
