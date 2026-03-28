import { LogIn } from "lucide-react";
import { Handle, Position } from "reactflow";

const MassBoundaryNode = ({ data }: { data: any }) => {
	return (
		<div
			className={`px-4 py-2 shadow-md rounded-md bg-white border-2 border-orange-400`}
		>
			<div className="flex items-center">
				<div className="rounded-full w-10 h-10 flex justify-center items-center bg-orange-100">
					<LogIn size={20} className="text-orange-600" />
				</div>
				<div className="ml-2">
					<div className="text-sm font-bold">Mass Boundary</div>
					<div className="text-xs text-gray-400">
						Fixed Mdot: {data.m_dot} kg/s
					</div>
				</div>
			</div>

			<Handle
				type="source"
				position={Position.Bottom}
				className="w-2 h-2 !bg-orange-500"
			/>
		</div>
	);
};

export default MassBoundaryNode;
