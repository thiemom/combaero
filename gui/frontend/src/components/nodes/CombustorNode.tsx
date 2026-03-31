import { Flame } from "lucide-react";
import { Handle, Position } from "reactflow";

const CombustorNode = ({ data }: { data: any }) => {
	const isSolved = !!data.result;

	return (
		<div
			className={`px-4 py-2 shadow-md rounded-md bg-white border-2 ${isSolved ? "border-red-500" : "border-red-300"}`}
		>
			<div className="flex items-center">
				<div className="rounded-full w-10 h-10 flex justify-center items-center bg-red-100">
					<Flame size={20} className="text-red-600" />
				</div>
				<div className="ml-2">
					<div className="text-sm font-bold">
						{data.label ? `${data.label} (Combustor)` : "Combustor"}
					</div>
					<div className="text-xs text-gray-500">
						Method: {data.method || "complete"}
					</div>
				</div>
			</div>

			{isSolved && data.result.state && (
				<div className="mt-2 text-xs font-mono bg-red-50 p-1 rounded border border-red-200">
					<div className="flex justify-between">
						<span className="text-red-600">P:</span>
						<span>{(data.result.state.P / 1e5).toFixed(2)} bar</span>
					</div>
					<div className="flex justify-between">
						<span className="text-red-600">T:</span>
						<span>{data.result.state.T.toFixed(1)} K</span>
					</div>
				</div>
			)}

			<Handle
				type="target"
				position={Position.Left}
				id="t-left"
				className="w-2 h-2 !bg-red-500"
			/>
			<Handle
				type="source"
				position={Position.Right}
				id="s-right"
				className="w-2 h-2 !bg-red-500"
			/>
		</div>
	);
};

export default CombustorNode;
