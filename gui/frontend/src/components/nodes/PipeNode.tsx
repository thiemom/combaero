import { ArrowRight } from "lucide-react";
import { Handle, Position } from "reactflow";

const PipeNode = ({ data }: { data: any }) => {
	const isSolved = !!data.result;

	return (
		<div
			className={`px-3 py-1 shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 ${isSolved ? "border-green-400" : "border-stone-300"}`}
			style={{ minWidth: "120px" }}
		>
			<Handle
				type="target"
				position={Position.Left}
				className="w-2 h-2 !bg-stone-400"
			/>

			<div className="flex flex-col items-center justify-center p-1 bg-white rounded border border-stone-200">
				<ArrowRight size={14} className="text-stone-500" />
			</div>

			<div className="flex flex-col">
				<div className="text-[10px] font-bold text-gray-400 uppercase leading-none">
					Pipe
				</div>
				<div className="text-xs font-semibold truncate max-w-[80px]">
					{data.L}m / {data.D}m
				</div>
			</div>

			<Handle
				type="source"
				position={Position.Right}
				className="w-2 h-2 !bg-stone-400"
			/>
		</div>
	);
};

export default PipeNode;
