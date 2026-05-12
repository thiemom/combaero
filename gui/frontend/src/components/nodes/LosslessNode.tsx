import { Link } from "lucide-react";
import { Handle, Position } from "reactflow";

const LosslessNode = ({ data }: { data: any }) => {
	const isSolved = !!data.result;

	return (
		<div
			className={`px-3 py-1 shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 ${isSolved ? "border-green-400" : "border-stone-300"}`}
			style={{ minWidth: "120px" }}
		>
			<div className="flex flex-col items-center justify-center p-1 bg-white rounded border border-stone-200">
				<Link size={14} className="text-slate-500" />
			</div>

			<div className="flex flex-col">
				<div className="text-[10px] font-bold text-gray-400 uppercase leading-none">
					Lossless
				</div>
				<div className="text-xs font-semibold truncate max-w-[80px]">
					Connection
				</div>
			</div>

			<div className="absolute -left-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-blue-500">
				in
			</div>
			<div className="absolute -right-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-amber-500">
				out
			</div>

			<Handle
				type="target"
				position={Position.Left}
				className="w-2 h-2 !bg-stone-400"
			/>
			<Handle
				type="source"
				position={Position.Right}
				className="w-2 h-2 !bg-stone-400"
			/>
		</div>
	);
};

export default LosslessNode;
