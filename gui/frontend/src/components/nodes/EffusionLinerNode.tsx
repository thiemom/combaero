import { Grip } from "lucide-react";
import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { handleStyle, rotPos } from "../../utils/nodeUtils";
import { NodeDiagRows } from "../NodeDiagRows";

// Port colours: blue = inlet, green = outlet (matches handle triangle colours).
const INLET = "#3b82f6";
const OUTLET = "#22c55e";

// Effusion liner (#471): a backside coolant duct whose wall bleeds through
// effusion holes into the hot gas. Three named ports:
//   coolant in  -> left  (target)
//   coolant out -> right (source)
//   discharge   -> top   (source), into the gas: a chamber's "s" handle
// The backend expands it into N stations (bleed + duct-fed panel).
const EffusionLinerNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	const inBase = Position.Left;
	const outBase = Position.Right;
	const dischargeBase = Position.Top;

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

	const bleed = data.result?.bleed_fraction;
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
				width: 150,
				height: 56,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center p-1 rounded border shrink-0 bg-orange-50 border-orange-200 text-orange-400">
				<Grip size={16} />
			</div>

			<div
				className="flex flex-col items-start flex-1 min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap">
					{data.label ? data.label : "Effusion Liner"}
				</div>
				<div className="text-[9px] font-mono whitespace-nowrap text-gray-500">
					{isSolved && bleed != null
						? `bleed ${(100 * Number(bleed)).toFixed(0)}%`
						: `${data.n_segments ?? 4} stations`}
				</div>
				{isSolved && <NodeDiagRows result={data.result} maxRows={1} />}
			</div>

			{/* Port labels: C in (left), C out (right), G discharge (top). */}
			<div
				className="absolute -left-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5"
				style={{ color: INLET, transform: `rotate(${-rotation}deg)` }}
			>
				C
			</div>
			<div
				className="absolute -right-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5"
				style={{ color: OUTLET, transform: `rotate(${-rotation}deg)` }}
			>
				C
			</div>
			<div
				className="absolute -top-3 left-1/2 -translate-x-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5"
				style={{ color: OUTLET, transform: `rotate(${-rotation}deg)` }}
			>
				G
			</div>

			<Handle
				type="target"
				position={rotPos(inBase, rotation)}
				style={handleStyle(inBase, rotation)}
				id="port-coolantin-target"
				title="Coolant in (backside duct)"
			/>
			<Handle
				type="source"
				position={rotPos(outBase, rotation)}
				style={handleStyle(outBase, rotation)}
				id="port-coolantout-source"
				title="Coolant out (backside duct)"
			/>
			<Handle
				type="source"
				position={rotPos(dischargeBase, rotation)}
				style={handleStyle(dischargeBase, rotation)}
				id="port-discharge-source"
				title="Effusion discharge into the gas (a chamber's s handle)"
			/>
		</div>
	);
};

export default EffusionLinerNode;
