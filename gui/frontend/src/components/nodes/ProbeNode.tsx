import { memo, useMemo } from "react";
import { Handle, Position } from "reactflow";
import useStore from "../../store/useStore";

// Quantity catalogue: key → {label, unit, format}
const QUANTITY_CATALOGUE: Record<
	string,
	{ label: string; unit: string; format: (v: number) => string }
> = {
	// Node state keys (accessed via result.state.*)
	"state.m_dot": { label: "ṁ", unit: "kg/s", format: (v) => v.toFixed(4) },
	"state.P": { label: "P", unit: "Pa", format: (v) => v.toFixed(0) },
	"state.P_total": {
		label: "Pt",
		unit: "Pa",
		format: (v) => v.toFixed(0),
	},
	"state.T": { label: "T", unit: "K", format: (v) => v.toFixed(2) },
	"state.T_total": { label: "Tt", unit: "K", format: (v) => v.toFixed(2) },
	"state.mach": { label: "Ma", unit: "—", format: (v) => v.toFixed(4) },
	"state.a": { label: "a", unit: "m/s", format: (v) => v.toFixed(2) },
	"state.rho": { label: "ρ", unit: "kg/m³", format: (v) => v.toFixed(4) },
	"state.h": { label: "h", unit: "J/kg", format: (v) => v.toFixed(0) },
	"state.s": { label: "s", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	"state.cp": { label: "cp", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	"state.gamma": { label: "γ", unit: "—", format: (v) => v.toFixed(4) },
	"state.mu": { label: "μ", unit: "Pa·s", format: (v) => v.toExponential(2) },
	"state.Pr": { label: "Pr", unit: "—", format: (v) => v.toFixed(3) },
	"state.Re": {
		label: "Re",
		unit: "—",
		format: (v) => v.toLocaleString("en", { maximumFractionDigits: 0 }),
	},
	"state.velocity": { label: "v", unit: "m/s", format: (v) => v.toFixed(3) },
	"state.phi": { label: "Φ", unit: "—", format: (v) => v.toFixed(3) },
	"state.mw": { label: "MW", unit: "g/mol", format: (v) => v.toFixed(2) },
	"state.nu": { label: "ν", unit: "m²/s", format: (v) => v.toExponential(2) },
	"state.u": { label: "u", unit: "J/kg", format: (v) => v.toFixed(0) },
	"state.cv": { label: "cv", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	// Element result keys (accessed via result.*)
	m_dot: { label: "ṁ", unit: "kg/s", format: (v) => v.toFixed(4) },
	mach: { label: "Ma", unit: "—", format: (v) => v.toFixed(4) },
	mach_in: { label: "Ma in", unit: "—", format: (v) => v.toFixed(4) },
	mach_out: { label: "Ma out", unit: "—", format: (v) => v.toFixed(4) },
	phi: { label: "Φ", unit: "—", format: (v) => v.toFixed(3) },
	p_ratio_total: { label: "P ratio", unit: "—", format: (v) => v.toFixed(4) },
	p_ratio: { label: "P ratio", unit: "—", format: (v) => v.toFixed(4) },
	Re: {
		label: "Re",
		unit: "—",
		format: (v) => v.toLocaleString("en", { maximumFractionDigits: 0 }),
	},
	Cd: { label: "Cd", unit: "—", format: (v) => v.toFixed(4) },
	h: { label: "h", unit: "J/kg", format: (v) => v.toFixed(0) },
	s: { label: "s", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	rho: { label: "ρ", unit: "kg/m³", format: (v) => v.toFixed(4) },
	a: { label: "a", unit: "m/s", format: (v) => v.toFixed(2) },
	mu: { label: "μ", unit: "Pa·s", format: (v) => v.toExponential(2) },
	Pr: { label: "Pr", unit: "—", format: (v) => v.toFixed(3) },
};

function resolveQuantity(result: any, key: string): number | undefined {
	if (!result) return undefined;
	if (key.startsWith("state.")) {
		const subKey = key.slice(6);
		return result.state?.[subKey];
	}
	return result[key];
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
