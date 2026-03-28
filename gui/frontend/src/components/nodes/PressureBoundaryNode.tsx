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

			<Handle
				type="target"
				position={Position.Top}
				className="w-2 h-2 !bg-blue-500"
			/>
		</div>
	);
};

export default PressureBoundaryNode;
