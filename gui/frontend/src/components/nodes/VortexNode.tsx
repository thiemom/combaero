import { Tornado } from "lucide-react";
import { Handle, type NodeProps, Position } from "reactflow";
import { NodeDiagRows } from "../NodeDiagRows";

const VortexNode = ({ data, selected }: NodeProps) => {
	const isSolved = !!data.result;
	const hasLocalRpm = data.omega_rpm != null;

	return (
		<div
			className={`shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 px-3 py-1 ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-green-400"
						: "border-stone-300"
			}`}
			style={{ width: 140, height: 56 }}
		>
			<div className="flex items-center justify-center p-1 bg-violet-50 border border-violet-200 rounded shrink-0">
				<Tornado size={16} className="text-violet-500" />
			</div>

			<div className="flex flex-col items-start flex-1 min-w-0">
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap">
					{data.label ? data.label : "Vortex"}
				</div>
				<div className="text-[9px] font-mono text-violet-500 whitespace-nowrap">
					{hasLocalRpm ? `${data.omega_rpm} rpm` : "global rpm"}
				</div>
				{isSolved && <NodeDiagRows result={data.result} maxRows={1} />}
			</div>

			<div className="absolute -left-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-blue-500">
				i
			</div>
			<div className="absolute -right-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-green-500">
				o
			</div>

			<Handle type="target" position={Position.Left} id="flow-target" />
			<Handle type="source" position={Position.Right} id="flow-source" />
		</div>
	);
};

export default VortexNode;
