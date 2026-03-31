import { Wind } from "lucide-react";
import { Handle, Position } from "reactflow";

const MomentumChamberNode = ({ data }: { data: any }) => {
	const isSolved = !!data.result;

	return (
		<div
			className={`px-4 py-2 shadow-md rounded-md bg-white border-2 ${isSolved ? "border-indigo-500" : "border-indigo-300"}`}
		>
			<div className="flex items-center">
				<div className="rounded-full w-10 h-10 flex justify-center items-center bg-indigo-100">
					<Wind size={20} className="text-indigo-600" />
				</div>
				<div className="ml-2">
					<div className="text-sm font-bold">
						{data.label ? `${data.label} (Momentum)` : "Momentum"}
					</div>
					<div className="text-xs text-gray-500">
						Area: {data.area?.toFixed(4)} m²
					</div>
				</div>
			</div>

			{isSolved && data.result.state && (
				<div className="mt-2 text-xs font-mono bg-indigo-50 p-1 rounded border border-indigo-200">
					<div className="flex justify-between">
						<span className="text-indigo-600">P:</span>
						<span>{(data.result.state.P / 1e5).toFixed(2)} bar</span>
					</div>
					<div className="flex justify-between">
						<span className="text-indigo-600">Mach:</span>
						<span>{data.result.state.mach?.toFixed(3) || "N/A"}</span>
					</div>
				</div>
			)}

			<Handle
				type="target"
				position={Position.Left}
				id="t-left"
				className="w-2 h-2 !bg-indigo-500"
			/>
			<Handle
				type="source"
				position={Position.Right}
				id="s-right"
				className="w-2 h-2 !bg-indigo-500"
			/>

			{/* Thermal handles */}
			<Handle
				type="target"
				position={Position.Top}
				id="t-top"
				className="w-2 h-2 !bg-orange-500 rounded-none transform rotate-45"
			/>
			<Handle
				type="source"
				position={Position.Bottom}
				id="s-bottom"
				className="w-2 h-2 !bg-orange-500 rounded-none transform rotate-45"
			/>
		</div>
	);
};

export default MomentumChamberNode;
