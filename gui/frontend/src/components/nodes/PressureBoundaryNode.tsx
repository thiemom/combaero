import { LogOut } from "lucide-react";
import { Handle, Position } from "reactflow";

const PressureBoundaryNode = ({ data }: { data: any }) => {
	return (
		<div
			className={`px-4 py-2 shadow-md rounded-md bg-white border-2 border-blue-400`}
		>
			<div className="flex items-center">
				<div className="rounded-full w-10 h-10 flex justify-center items-center bg-blue-100">
					<LogOut size={20} className="text-blue-600" />
				</div>
				<div className="ml-2">
					<div className="text-sm font-bold">Pressure Boundary</div>
					<div className="text-xs text-gray-400">
						P_total: {(data.P_total / 1e5).toFixed(2)} bar
					</div>
				</div>
			</div>

			{data.result?.state && (
				<div className="mt-2 text-xs font-mono bg-green-50 p-1 rounded border border-green-200">
					<div className="flex justify-between">
						<span className="text-green-600">T:</span>
						<span>{data.result.state.T.toFixed(1)} K</span>
					</div>
				</div>
			)}

			{/* Cardinal Target Handles */}
			<Handle
				type="target"
				position={Position.Top}
				id="t-top"
				className="w-2 h-2 !bg-blue-500"
			/>
			<Handle
				type="target"
				position={Position.Bottom}
				id="t-bottom"
				className="w-2 h-2 !bg-blue-500"
			/>
			<Handle
				type="target"
				position={Position.Left}
				id="t-left"
				className="w-2 h-2 !bg-blue-500"
			/>
			<Handle
				type="target"
				position={Position.Right}
				id="t-right"
				className="w-2 h-2 !bg-blue-500"
			/>
		</div>
	);
};

export default PressureBoundaryNode;
