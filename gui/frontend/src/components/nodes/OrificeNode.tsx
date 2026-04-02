import { ChevronRight } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { getHandlePosition, getWrapperSize } from "../../utils/nodeUtils";

const NODE_W = 120;
const NODE_H = 48;
const WRAPPER = getWrapperSize(NODE_W, NODE_H);

const OrificeNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	return (
		<div
			className="relative"
			style={{
				width: WRAPPER,
				height: WRAPPER,
				display: "flex",
				alignItems: "center",
				justifyContent: "center",
			}}
		>
			<div
				className={`shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 px-3 py-1 ${
					selected
						? "border-blue-500 shadow-blue-100"
						: isSolved
							? "border-green-400"
							: "border-stone-300"
				}`}
				style={{
					width: NODE_W,
					height: NODE_H,
					transform: `rotate(${rotation}deg)`,
					transformOrigin: "center center",
				}}
			>
				<div className="flex items-center justify-center p-1 bg-white rounded border border-stone-200 shrink-0">
					<ChevronRight size={14} className="text-orange-400" />
				</div>

				<div
					className="flex flex-col items-start"
					style={{ transform: `rotate(${-rotation}deg)` }}
				>
					<div className="text-[10px] font-bold text-gray-400 uppercase leading-none whitespace-nowrap">
						{data.label ? data.label : "Orifice"}
					</div>
					<div className="text-[10px] font-semibold whitespace-nowrap">
						Cd: {data.Cd}
					</div>
				</div>
			</div>

			<Handle
				type="target"
				position={getHandlePosition(Position.Left, rotation)}
				id="flow-target"
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Right, rotation)}
				id="flow-source"
			/>
			<Handle
				type="target"
				position={getHandlePosition(Position.Top, rotation)}
				id="thermal-target"
			/>
			<Handle
				type="source"
				position={getHandlePosition(Position.Bottom, rotation)}
				id="thermal-source"
			/>
		</div>
	);
};

export default OrificeNode;
