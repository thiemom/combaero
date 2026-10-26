import { memo, useMemo } from "react";
import { Handle, Position } from "reactflow";
import useStore from "../../store/useStore";
import { isConvergedResult, resolveKey } from "../../utils/diagnostics";
import { nodeTypeLabel } from "../../utils/nodeUtils";
import { withPortArm } from "../../utils/portArms";
import { quantityMeta } from "../../utils/quantities";
import NodeIcon from "../NodeIcon";

function ProbeSlot({
	result,
	quantityKey,
	targetType,
	targetData,
}: {
	result: any;
	quantityKey: string | undefined;
	targetType: string | undefined;
	targetData: any;
}) {
	if (!quantityKey) {
		return (
			<div className="flex justify-between text-[10px] font-mono text-stone-400">
				<span className="italic">— select quantity —</span>
			</div>
		);
	}

	const meta = quantityMeta(quantityKey);
	const value = resolveKey(result, quantityKey);

	return (
		<div className="flex justify-between items-baseline gap-2 text-[11px] font-mono">
			<span className="text-stone-500 min-w-[24px]">
				{withPortArm(meta.label, quantityKey, targetType, targetData)}
			</span>
			<span className="font-bold text-right">
				{value !== undefined ? meta.format(value) : "—"}
			</span>
			<span className="text-stone-400 text-[9px] font-light">
				{meta.unit === "—" ? "" : meta.unit}
			</span>
		</div>
	);
}

function ProbeNode({ data, selected }: { data: any; selected?: boolean }) {
	const { nodes, edges } = useStore();

	const targetElement = useMemo(
		() =>
			nodes.find((n) => n.id === data.target_id) ||
			edges.find((e) => e.id === data.target_id),
		[nodes, edges, data.target_id],
	);

	const result = targetElement?.data?.result;
	const label = data.label || "Probe";
	const isWall = (targetElement as any)?.data?.type === "thermal";
	const targetLabel = isWall
		? (targetElement?.data?.label ?? "Thermal Wall")
		: targetElement?.data?.label ||
			nodeTypeLabel(targetElement?.type) ||
			(data.target_id ? "Unknown target" : "No target");
	const hasTarget = !!data.target_id;
	const missing = hasTarget && !targetElement;
	// The backend returns the last iterate of a solve that did not converge
	// (success: false): shown, but flagged so it does not read as an answer.
	const notConverged = !!result && !isConvergedResult(result);

	return (
		<div
			className={`relative bg-white border rounded-lg shadow-sm min-w-[140px] text-left select-none transition-shadow ${
				selected
					? "border-blue-400 shadow-blue-100 shadow-md"
					: "border-stone-200"
			}`}
		>
			{/* Header */}
			<div className="flex items-center gap-1.5 px-2.5 pt-2 pb-1">
				<NodeIcon type="probe" size={10} />
				<span className="text-[10px] font-bold text-stone-700 truncate">
					{label}
				</span>
			</div>

			{/* Target badge */}
			<div className="px-2.5 pb-1 flex flex-wrap gap-1">
				<span
					className={`text-[8px] px-1.5 py-0.5 rounded-full ${
						hasTarget && !missing
							? "bg-stone-100 text-stone-500"
							: "bg-amber-50 text-amber-500"
					}`}
				>
					{!hasTarget
						? "no target"
						: missing
							? "target deleted"
							: `→ ${targetLabel}`}
				</span>
				{notConverged && (
					<span
						className="text-[8px] px-1.5 py-0.5 rounded-full bg-red-50 text-red-500"
						title="The last solve did not converge: these are its last iterate, not a solution."
					>
						not converged
					</span>
				)}
			</div>

			{/* Divider */}
			<div className="border-t border-stone-100 mx-2 mb-1" />

			{/* Quantity slots */}
			<div
				className={`px-2.5 pb-2 flex flex-col gap-1 ${notConverged ? "opacity-60" : ""}`}
			>
				<ProbeSlot
					result={result}
					quantityKey={data.slot1_key}
					targetType={targetElement?.type}
					targetData={targetElement?.data}
				/>
				<ProbeSlot
					result={result}
					quantityKey={data.slot2_key}
					targetType={targetElement?.type}
					targetData={targetElement?.data}
				/>
			</div>

			{/* Drawing a connection to this handle sets the probe's target
			    (useStore.onConnect); it never becomes a flow edge. */}
			<Handle
				type="target"
				position={Position.Left}
				id="probe-in"
				style={{ opacity: 0, width: 6, height: 6 }}
			/>
		</div>
	);
}

export default memo(ProbeNode);
