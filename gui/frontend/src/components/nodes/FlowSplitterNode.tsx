import { useEffect } from "react";
import {
	Handle,
	type NodeProps,
	Position,
	useUpdateNodeInternals,
} from "reactflow";
import { handleStyle, rotPos } from "../../utils/nodeUtils";
import {
	resolvedFractions,
	SPLITTER_WIDTH,
	splitterFractions,
	splitterHeight,
	splitterReference,
} from "../../utils/splitter";
import NodeIcon from "../NodeIcon";

const OUTLET = "#22c55e";

/**
 * A flow splitter (#903): one inlet (left) to N outlets (right), each at a
 * given fraction of the throughput, a mass-flow controller per outlet. The
 * card grows with its outlets; each outlet's fraction is drawn at its
 * handle, the remainder as "rest", the reference (valve open) with "open".
 * After a solve each metered outlet shows the drop its valve takes. No
 * surface: no thermal handles.
 */
const FlowSplitterNode = ({ id, data, selected }: NodeProps) => {
	const rotation = data.rotation || 0;
	const updateNodeInternals = useUpdateNodeInternals();
	const textRotation = rotation === 90 || rotation === 180 ? 180 : 0;
	const fractions = splitterFractions(data);
	const n = fractions.length;
	const shown = resolvedFractions(fractions);
	const ref = splitterReference(data);
	const dPt: (number | null)[] | undefined = data.result?.rows_valve_dPt;

	// biome-ignore lint/correctness/useExhaustiveDependencies: rotation and the outlet count move the handles
	useEffect(() => {
		updateNodeInternals(id);
	}, [id, rotation, n, updateNodeInternals]);

	const top = (i: number) => `calc(${(100 * (i + 0.5)) / n}% - 6px)`;
	return (
		<div
			className={`shadow-sm rounded bg-stone-50 border-2 flex items-center gap-2 pl-3 pr-12 py-1 ${
				selected
					? "border-blue-500 shadow-blue-100"
					: data.result
						? "border-green-400"
						: "border-stone-300"
			}`}
			style={{
				width: SPLITTER_WIDTH,
				height: splitterHeight(n),
				transform: `rotate(${rotation}deg)`,
				transformOrigin: "center center",
			}}
		>
			<div className="flex items-center justify-center p-1 bg-white rounded border border-stone-200 shrink-0">
				<NodeIcon type="flow_splitter" size={14} />
			</div>
			<div
				className="flex flex-col items-start flex-1 min-w-0 overflow-hidden"
				style={{ transform: `rotate(${textRotation}deg)` }}
			>
				<div className="text-[10px] font-bold uppercase leading-none whitespace-nowrap truncate max-w-full">
					{data.label ? data.label : "Splitter"}
				</div>
				<div className="text-[9px] font-mono text-gray-500 whitespace-nowrap truncate max-w-full">
					{typeof data.result?.m_dot === "number"
						? `${Number(data.result.m_dot).toPrecision(4)} kg/s`
						: `${n} outlets`}
				</div>
			</div>

			{fractions.map((f, i) => {
				const drop = dPt?.[i];
				const text =
					i === ref
						? "open"
						: typeof drop === "number"
							? `${(drop / 1000).toPrecision(3)} kPa`
							: "";
				return (
					<div
						key={`label_out${i}`}
						className="absolute right-1.5 text-[8px] font-mono leading-none text-right select-none pointer-events-none"
						style={{
							top: `calc(${(100 * (i + 0.5)) / n}% - 8px)`,
							transform: `rotate(${textRotation}deg)`,
						}}
						title={
							i === ref
								? "The reference outlet: its valve is open (Pt out = Pt in)"
								: "The drop this outlet's valve takes"
						}
					>
						<div style={{ color: OUTLET }} className="font-bold">
							{f == null ? "rest" : `${(100 * shown[i]).toPrecision(3)}%`}
						</div>
						<div className="text-gray-500">{text}</div>
					</div>
				);
			})}

			<Handle
				type="target"
				position={rotPos(Position.Left, rotation)}
				style={handleStyle(Position.Left, rotation)}
				id="port-in-target"
				title="Inlet"
			/>
			{fractions.map((_, i) => (
				<Handle
					key={`port-out${i}`}
					type="source"
					position={rotPos(Position.Right, rotation)}
					style={{ ...handleStyle(Position.Right, rotation), top: top(i) }}
					id={`port-out${i}-source`}
					title={`Outlet ${i + 1}`}
				/>
			))}
		</div>
	);
};

export default FlowSplitterNode;
