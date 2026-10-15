import { Maximize2, Minimize2 } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import useStore from "../../store/useStore";
import { inheritedField } from "../../utils/inheritance";
import { iconClass } from "../../utils/nodeAppearance";
import { handleStyle, rotPos } from "../../utils/nodeUtils";
import { NodeDiagRows } from "../NodeDiagRows";

const AreaChangeNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	// Unset faces are inherited: show the solver's resolved areas (#559), not
	// placeholders (the old 0.01 / 0.02 always read "F1/F0: 2.00").
	const previewItems = useStore((s) => s.previewItems);
	const f0: number | null =
		data.F0 ?? inheritedField(previewItems, id, "F0")?.value ?? null;
	const f1: number | null =
		data.F1 ?? inheritedField(previewItems, id, "F1")?.value ?? null;
	const ratio: number | null =
		data.result?.ratio ?? (f0 && f1 ? f1 / f0 : null);
	const isExpansion = ratio == null || ratio >= 1;
	const Icon = isExpansion ? Maximize2 : Minimize2;

	return (
		<div
			className={`shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 px-3 py-1 ${
				selected
					? "border-blue-500 shadow-blue-100"
					: isSolved
						? "border-green-400"
						: "border-stone-300"
			}`}
			style={{
				width: 140,
				height: 48,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center p-1 bg-white rounded border border-stone-200 shrink-0">
				<Icon
					size={12}
					className={
						data.model_type === "conical"
							? "text-purple-500"
							: iconClass("area_change")
					}
				/>
			</div>

			<div
				className="flex flex-col items-start flex-1 min-w-0 overflow-hidden"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap truncate max-w-full">
					{data.label ? data.label : "Area Change"}
				</div>
				<div className="flex flex-row gap-1 text-[9px] font-mono text-gray-500 whitespace-nowrap truncate max-w-full">
					<div>F₁/F₀: {ratio != null ? ratio.toFixed(2) : "—"}</div>
				</div>
				{isSolved && <NodeDiagRows result={data.result} maxRows={1} />}
			</div>

			<div
				className="absolute -left-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-blue-500"
				style={{ transform: `rotate(${-rotation}deg)` }}
			>
				i
			</div>
			<div
				className="absolute -right-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-green-500"
				style={{ transform: `rotate(${-rotation}deg)` }}
			>
				o
			</div>
			<div
				className="absolute -top-3 left-1/2 -translate-x-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-orange-400"
				style={{ transform: `rotate(${-rotation}deg)` }}
			>
				q
			</div>
			<div
				className="absolute -bottom-3 left-1/2 -translate-x-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5 text-orange-400"
				style={{ transform: `rotate(${-rotation}deg)` }}
			>
				q
			</div>

			<Handle
				type="target"
				position={rotPos(Position.Left, rotation)}
				style={handleStyle(Position.Left, rotation)}
				id="flow-target"
			/>
			<Handle
				type="source"
				position={rotPos(Position.Right, rotation)}
				style={handleStyle(Position.Right, rotation)}
				id="flow-source"
			/>
			<Handle
				type="target"
				position={rotPos(Position.Top, rotation)}
				style={handleStyle(Position.Top, rotation)}
				id="thermal-target"
			/>
			<Handle
				type="source"
				position={rotPos(Position.Bottom, rotation)}
				style={handleStyle(Position.Bottom, rotation)}
				id="thermal-source"
			/>
		</div>
	);
};

export default AreaChangeNode;
