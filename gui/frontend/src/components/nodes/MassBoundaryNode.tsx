import { ArrowRightToLine } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { getHandlePosition } from "../../utils/nodeUtils";

const MassBoundaryNode = ({ data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;

	// Orientation logic
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	return (
		<div
			className={`w-36 px-4 py-2 shadow-md rounded border-2 bg-stone-50 transition-colors flex flex-col items-center justify-center ${
				selected ? "border-blue-500 shadow-blue-100" : "border-stone-400"
			}`}
			style={{ minHeight: "80px" }}
		>
			<div
				className="flex items-center justify-center pointer-events-none mb-2"
				style={{
					transform: `rotate(${rotation}deg)`,
				}}
			>
				<div
					className="flex flex-col items-center"
					style={{ transform: `rotate(${textRotation}deg)` }}
				>
					<ArrowRightToLine size={24} className="text-blue-600 mb-1" />
					<div className="text-xs font-bold text-center">
						{data.label ? data.label : "Mass Flow"}
					</div>
					<div className="text-[10px] text-gray-500 font-mono">
						{data.m_dot} kg/s
					</div>
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
