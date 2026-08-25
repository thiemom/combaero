import { Wind } from "lucide-react";
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

// Critical-mode supersonic ejector: fixed 3-port topology (no merge/branch
// flip). Primary (motive) + secondary (suction) are inlets; one outlet.
//   Primary  -> left  (target)
//   Secondary-> bottom(target)
//   Outlet   -> right (source)
const EjectorNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const isSolved = !!data.result;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;

	const primaryBase = Position.Left;
	const secondaryBase = Position.Bottom;
	const outletBase = Position.Right;
	const primaryPos = rotPos(primaryBase, rotation);
	const secondaryPos = rotPos(secondaryBase, rotation);
	const outletPos = rotPos(outletBase, rotation);

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation triggers handle re-measurement
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, updateNodeInternals]);

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
				height: 56,
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center p-1 rounded border shrink-0 bg-sky-50 border-sky-200 text-sky-500">
				<Wind size={16} />
			</div>

			<div
				className="flex flex-col items-start flex-1 min-w-0"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap">
					{data.label ? data.label : "Ejector"}
				</div>
				<div className="text-[9px] font-mono whitespace-nowrap text-sky-500">
					{isSolved && data.result?.omega != null
						? `ω ${Number(data.result.omega).toFixed(3)}`
						: "critical"}
				</div>
				{isSolved && <NodeDiagRows result={data.result} maxRows={1} />}
			</div>

			{/* Port labels: P (primary) left, S (secondary) bottom, O (outlet) right. */}
			<div
				className="absolute -left-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5"
				style={{ color: INLET, transform: `rotate(${-rotation}deg)` }}
			>
				P
			</div>
			<div
				className="absolute -bottom-3 left-1/2 -translate-x-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5"
				style={{ color: INLET, transform: `rotate(${-rotation}deg)` }}
			>
				S
			</div>
			<div
				className="absolute -right-3 top-1/2 -translate-y-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5"
				style={{ color: OUTLET, transform: `rotate(${-rotation}deg)` }}
			>
				O
			</div>

			{/* Primary inlet (motive) -- left */}
			<Handle
				type="target"
				position={primaryPos}
				style={handleStyle(primaryBase, rotation)}
				id="port-primary-target"
			/>
			{/* Secondary inlet (suction) -- bottom */}
			<Handle
				type="target"
				position={secondaryPos}
				style={handleStyle(secondaryBase, rotation)}
				id="port-secondary-target"
			/>
			{/* Outlet -- right */}
			<Handle
				type="source"
				position={outletPos}
				style={handleStyle(outletBase, rotation)}
				id="port-outlet-source"
			/>
		</div>
	);
};

export default EjectorNode;
