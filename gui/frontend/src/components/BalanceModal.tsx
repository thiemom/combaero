import { X } from "lucide-react";
import type React from "react";
import useStore from "../store/useStore";
import {
	balanceIsOpen,
	describeElement,
	displayName,
	type Ledger,
	type SecondLawFlag,
} from "../utils/secondLaw";

/**
 * The network's control volume after a solve (#545): what enters and leaves
 * at each boundary, the heat terms, and how well mass and energy close; and
 * the elements whose control volume destroys entropy (creates total
 * pressure). Both are reports: the solve's success does not depend on them.
 */
const BalanceModal: React.FC<{
	ledger: Ledger;
	failures: SecondLawFlag[];
	onClose: () => void;
}> = ({ ledger, failures, onClose }) => {
	const { nodes } = useStore();
	const name = (id: string) => displayName(id, nodes);
	const kW = (w: number) => (w / 1000).toFixed(3);
	const open = balanceIsOpen(ledger);

	return (
		<div className="fixed inset-0 z-[200] flex items-center justify-center bg-black/40">
			<div
				role="dialog"
				aria-modal="true"
				aria-label="Network balance"
				className="bg-white rounded-lg shadow-xl border border-gray-200 p-4 w-[680px] max-w-[95vw] max-h-[85vh] overflow-y-auto"
				onKeyDown={(e) => e.key === "Escape" && onClose()}
			>
				<div className="flex items-center justify-between mb-3">
					<span className="font-semibold text-sm text-gray-700">
						Network balance
					</span>
					<button
						type="button"
						onClick={onClose}
						className="hover:bg-gray-100 p-1 rounded-full"
						aria-label="Close"
					>
						<X size={16} />
					</button>
				</div>

				<div className="overflow-x-auto mb-3">
					<table className="w-full table-fixed text-xs font-mono">
						<colgroup>
							<col className="w-[34%]" />
							<col className="w-[34%]" />
							<col className="w-[16%]" />
							<col className="w-[16%]" />
						</colgroup>
						<thead>
							<tr className="text-left text-gray-500 border-b">
								<th className="py-1 pr-2 font-semibold">Boundary</th>
								<th className="py-1 pr-2 font-semibold">via</th>
								<th className="py-1 font-semibold text-right whitespace-nowrap">
									ṁ in [kg/s]
								</th>
								<th className="py-1 font-semibold text-right whitespace-nowrap">
									H in [kW]
								</th>
							</tr>
						</thead>
						<tbody>
							{ledger.streams.map((s) => (
								<tr
									key={`${s.node}-${s.via}`}
									className="border-b border-gray-100 align-top"
								>
									<td className="py-0.5 pr-2 break-words">{name(s.node)}</td>
									<td className="py-0.5 pr-2 break-words">
										{s.direct
											? `direct to ${name(s.via_owner)}`
											: name(s.via_owner)}
									</td>
									<td className="py-0.5 text-right whitespace-nowrap">
										{s.m_dot_in.toFixed(5)}
									</td>
									<td className="py-0.5 text-right whitespace-nowrap">
										{kW(s.H_in)}
									</td>
								</tr>
							))}
						</tbody>
					</table>
				</div>

				<div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs font-mono mb-3">
					<span className="text-gray-500">Heat given</span>
					<span className="text-right">{kW(ledger.Q_given)} kW</span>
					<span className="text-gray-500">Wall heat leaving at boundaries</span>
					<span className="text-right">{kW(ledger.Q_wall_out)} kW</span>
					<span className="text-gray-500">Heat withheld (stagnant nodes)</span>
					<span className="text-right">{kW(ledger.Q_withheld)} kW</span>
					<span className="text-gray-500">Mass closure (rel.)</span>
					<span className="text-right">
						{ledger.mass_closure_rel.toExponential(2)}
					</span>
					<span className="text-gray-500">Energy closure (rel.)</span>
					<span className="text-right">
						{ledger.energy_closure_rel.toExponential(2)}
					</span>
				</div>
				<p
					className={`text-[10px] leading-snug mb-3 ${open ? "text-amber-700" : "text-gray-500"}`}
				>
					{open
						? "The balance does not close: energy or mass is lost or created somewhere in the network. A node held at the temperature clamp discards energy this way."
						: "Closure: what enters plus the heat given, minus what leaves and the heat leaving or withheld. Enthalpy is absolute (formation included), so a combustor closes it as it is."}
				</p>

				<div className="text-xs font-semibold text-gray-700 mb-1">
					Second law per element
				</div>
				{failures.length === 0 ? (
					<p className="text-[10px] text-gray-500 leading-snug">
						No element creates total pressure: every element's control volume
						generates entropy, or is not checked (listed in its result).
					</p>
				) : (
					<ul className="text-[10px] text-amber-800 leading-snug list-disc pl-4">
						{failures.map((f) => (
							<li key={f.element_id}>
								{name(f.owner_id)}: {describeElement(f)}, σ/(R ṁ) ={" "}
								{f.sigma_norm.toExponential(2)}
							</li>
						))}
					</ul>
				)}
			</div>
		</div>
	);
};

export default BalanceModal;
