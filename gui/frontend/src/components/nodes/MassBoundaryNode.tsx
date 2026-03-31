import { Handle, Position } from "reactflow";

const MassBoundaryNode = ({ data }: { data: any }) => {
	return (
		<div
			className={`px-4 py-2 shadow-md rounded-md bg-white border-2 border-orange-400`}
		>
			<div className="flex items-center">
				<div className="rounded-full w-10 h-10 flex justify-center items-center bg-orange-100">
					<svg
						width="24"
						height="24"
						viewBox="0 0 24 24"
						fill="none"
						stroke="currentColor"
						strokeWidth="2"
						strokeLinecap="round"
						strokeLinejoin="round"
						className="text-orange-600"
					>
						<title>Mass Boundary Icon</title>
						<path d="M2 12h20" />
						<path d="M19 9l3 3-3 3" />
						<circle cx="5" cy="12" r="3" />
					</svg>
				</div>
				<div className="ml-2">
					<div className="text-sm font-bold">
						{data.label ? `${data.label} (Mass)` : "Mass Boundary"}
					</div>
					<div className="text-xs text-gray-400">
						Fixed Mdot: {data.m_dot} kg/s
					</div>
				</div>
			</div>

			{data.result?.state && (
				<div className="mt-2 text-xs font-mono bg-green-50 p-1 rounded border border-green-200">
					<div className="flex justify-between">
						<span className="text-green-600">P:</span>
						<span>{(data.result.state.P / 1e5).toFixed(2)} bar</span>
					</div>
					<div className="flex justify-between">
						<span className="text-green-600">T:</span>
						<span>{data.result.state.T.toFixed(1)} K</span>
					</div>
				</div>
			)}

			{/* Cardinal Source Handles */}
			<Handle
				type="source"
				position={Position.Top}
				id="s-top"
				className="w-2 h-2 !bg-orange-500"
			/>
			<Handle
				type="source"
				position={Position.Bottom}
				id="s-bottom"
				className="w-2 h-2 !bg-orange-500"
			/>
			<Handle
				type="source"
				position={Position.Left}
				id="s-left"
				className="w-2 h-2 !bg-orange-500"
			/>
			<Handle
				type="source"
				position={Position.Right}
				id="s-right"
				className="w-2 h-2 !bg-orange-500"
			/>
		</div>
	);
};

export default MassBoundaryNode;
