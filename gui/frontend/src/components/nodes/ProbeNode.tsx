import { memo, useMemo } from "react";
import { Handle, Position } from "reactflow";
import useStore from "../../store/useStore";

import { QUANTITY_CATALOGUE } from "../../utils/quantities";

function resolveQuantity(result: any, key: string): number | undefined {
	if (!result) return undefined;
	if (key.startsWith("state.")) {
		const subKey = key.slice(6);
		return result.state?.[subKey] ?? result[subKey];
	}
	return result[key] ?? result.state?.[key];
}

function ProbeSlot({
	result,
	quantityKey,
}: {
	result: any;
	quantityKey: string | undefined;
}) {
	if (!quantityKey) {
		return (
			<div className="flex justify-between text-[10px] font-mono text-stone-400">
				<span className="italic">— select quantity —</span>
			</div>
		);
	}

	const meta = QUANTITY_CATALOGUE[quantityKey];
	const value = resolveQuantity(result, quantityKey);

	return (
		<div className="flex justify-between items-baseline gap-2 text-[11px] font-mono">
			<span className="text-stone-500 min-w-[24px]">
				{meta?.label ?? quantityKey}
			</span>
			<span className="font-bold text-right">
				{value !== undefined && value !== null
					? (meta?.format(value) ?? String(value))
					: "—"}
			</span>
			<span className="text-stone-400 text-[9px]">{meta?.unit ?? ""}</span>
		</div>
	);
}

function ProbeNode({ data, selected }: { data: any; selected?: boolean }) {
	const { nodes } = useStore();

	const targetNode = useMemo(
		() => nodes.find((n) => n.id === data.target_id),
		[nodes, data.target_id],
	);

	const result = targetNode?.data?.result;
	const label = data.label || "Probe";
	const targetLabel =
		targetNode?.data?.label ||
		targetNode?.type?.replace(/_/g, " ") ||
		(data.target_id ? "Unknown target" : "No target");
	const hasTarget = !!data.target_id;

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
				<span className="text-[9px] text-blue-400 font-bold">⊙</span>
				<span className="text-[10px] font-bold text-stone-700 truncate">
					{label}
				</span>
			</div>

			{/* Target badge */}
			<div className="px-2.5 pb-1">
				<span
					className={`text-[8px] px-1.5 py-0.5 rounded-full ${
						hasTarget
							? "bg-stone-100 text-stone-500"
							: "bg-amber-50 text-amber-500"
					}`}
				>
					{hasTarget ? `→ ${targetLabel}` : "no target"}
				</span>
			</div>

			{/* Divider */}
			<div className="border-t border-stone-100 mx-2 mb-1" />

			{/* Quantity slots */}
			<div className="px-2.5 pb-2 flex flex-col gap-1">
				<ProbeSlot result={result} quantityKey={data.slot1_key} />
				<ProbeSlot result={result} quantityKey={data.slot2_key} />
			</div>

			{/* React Flow handle (invisible, for probe-link edges) */}
			<Handle
				type="target"
				position={Position.Left}
				id="probe-in"
				style={{ opacity: 0, width: 6, height: 6 }}
			/>
		</div>
	);
}

export const QUANTITY_CATALOGUE_KEYS = Object.keys(QUANTITY_CATALOGUE);
export { QUANTITY_CATALOGUE };
export default memo(ProbeNode);
