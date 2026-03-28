import { Database } from "lucide-react";
import { Handle, Position } from "reactflow";

const PlenumNode = ({ data }: { data: any }) => {
	const isSolved = !!data.result;

	return (
		<div
			className={`px-4 py-2 shadow-md rounded-md bg-white border-2 ${isSolved ? "border-green-500" : "border-stone-400"}`}
		>
			<div className="flex items-center">
				<div className="rounded-full w-10 h-10 flex justify-center items-center bg-stone-100">
					<Database size={20} className="text-stone-600" />
				</div>
				<div className="ml-2">
					<div className="text-lg font-bold">Plenum</div>
					<div className="text-gray-500">ID: {data.id}</div>
				</div>
			</div>

			{isSolved && (
				<div className="mt-2 text-xs font-mono bg-green-50 p-1 rounded border border-green-200">
					<div>P: {(data.result.P / 1e5).toFixed(2)} bar</div>
					<div>T: {data.result.T.toFixed(1)} K</div>
				</div>
			)}

			<Handle
				type="target"
				position={Position.Top}
				className="w-2 h-2 !bg-stone-500"
			/>
			<Handle
				type="source"
				position={Position.Bottom}
				className="w-2 h-2 !bg-stone-500"
			/>
		</div>
	);
};

export default PlenumNode;
