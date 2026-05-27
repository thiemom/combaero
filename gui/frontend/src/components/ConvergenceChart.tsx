import { X } from "lucide-react";
import type React from "react";

export interface ConvergencePoint {
	eval: number;
	t: number;
	norm: number;
}

export interface WorstResidual {
	name: string;
	residual: number;
}

interface Props {
	history: ConvergencePoint[];
	worstResiduals?: WorstResidual[];
	lmStartedAt?: number | null;
	onClose: () => void;
}

const W = 480;
const H = 200;
const PAD = { top: 12, right: 16, bottom: 36, left: 56 };
const CHART_W = W - PAD.left - PAD.right;
const CHART_H = H - PAD.top - PAD.bottom;

function toLog(v: number): number {
	return v > 0 ? Math.log10(v) : -20;
}

const ConvergenceChart: React.FC<Props> = ({
	history,
	worstResiduals,
	lmStartedAt,
	onClose,
}) => {
	if (!history || history.length === 0) return null;

	const norms = history.map((p) => p.norm).filter((n) => n > 0);
	const logMin = Math.floor(Math.min(...norms.map(toLog)));
	const logMax = Math.ceil(Math.max(...norms.map(toLog)));
	const logRange = logMax - logMin || 1;

	const tMax = history[history.length - 1].t || 1;
	const lmT =
		lmStartedAt != null
			? (history.find((p) => p.eval >= lmStartedAt)?.t ?? null)
			: null;

	const px = (t: number) => PAD.left + (t / tMax) * CHART_W;
	const py = (norm: number) =>
		PAD.top + CHART_H - ((toLog(norm) - logMin) / logRange) * CHART_H;

	const makePolyline = (pts: ConvergencePoint[]) =>
		pts
			.filter((p) => p.norm > 0)
			.map((p) => `${px(p.t).toFixed(1)},${py(p.norm).toFixed(1)}`)
			.join(" ");

	const hybrPts = lmT != null ? history.filter((p) => p.t <= lmT) : history;
	const lmPts = lmT != null ? history.filter((p) => p.t >= lmT) : [];

	const decades: number[] = [];
	for (let d = logMin; d <= logMax; d++) decades.push(d);

	const top5 = (worstResiduals ?? []).slice(0, 5);

	return (
		<div className="fixed inset-0 z-[200] flex items-center justify-center bg-black/40">
			<div
				role="dialog"
				aria-modal="true"
				aria-label="Convergence history"
				className="bg-white rounded-lg shadow-xl border border-gray-200 p-4 w-[540px] max-w-[95vw]"
				onKeyDown={(e) => e.key === "Escape" && onClose()}
			>
				<div className="flex items-center justify-between mb-3">
					<span className="font-semibold text-sm text-gray-700">
						Convergence history
					</span>
					<button
						type="button"
						onClick={onClose}
						className="hover:bg-gray-100 p-1 rounded-full transition-colors"
						aria-label="Close"
					>
						<X size={15} />
					</button>
				</div>

				<svg
					width={W}
					height={H}
					style={{ display: "block", overflow: "visible" }}
				>
					<title>Convergence history</title>
					{decades.map((d) => {
						const y = py(10 ** d);
						return (
							<g key={d}>
								<line
									x1={PAD.left}
									x2={PAD.left + CHART_W}
									y1={y}
									y2={y}
									stroke="#e5e7eb"
									strokeWidth={1}
								/>
								<text
									x={PAD.left - 4}
									y={y}
									textAnchor="end"
									dominantBaseline="middle"
									fontSize={9}
									fill="#9ca3af"
								>
									{d === 0 ? "1" : `1e${d}`}
								</text>
							</g>
						);
					})}

					{lmT != null && (
						<line
							x1={px(lmT)}
							x2={px(lmT)}
							y1={PAD.top}
							y2={PAD.top + CHART_H}
							stroke="#6366f1"
							strokeWidth={1}
							strokeDasharray="3 2"
						/>
					)}

					{hybrPts.length > 1 && (
						<polyline
							points={makePolyline(hybrPts)}
							fill="none"
							stroke="#f97316"
							strokeWidth={1.5}
							strokeLinejoin="round"
						/>
					)}
					{lmPts.length > 1 && (
						<polyline
							points={makePolyline(lmPts)}
							fill="none"
							stroke="#3b82f6"
							strokeWidth={1.5}
							strokeLinejoin="round"
						/>
					)}

					<text
						x={PAD.left + CHART_W / 2}
						y={H - 4}
						textAnchor="middle"
						fontSize={9}
						fill="#9ca3af"
					>
						t (s)
					</text>
					{[0, 0.25, 0.5, 0.75, 1.0].map((f) => (
						<text
							key={f}
							x={px(f * tMax)}
							y={PAD.top + CHART_H + 10}
							textAnchor="middle"
							fontSize={9}
							fill="#9ca3af"
						>
							{(f * tMax).toFixed(f === 0 ? 0 : 1)}
						</text>
					))}

					<text
						x={PAD.left - 44}
						y={PAD.top + CHART_H / 2}
						textAnchor="middle"
						fontSize={9}
						fill="#9ca3af"
						transform={`rotate(-90, ${PAD.left - 44}, ${PAD.top + CHART_H / 2})`}
					>
						||F||
					</text>
				</svg>

				<div className="flex items-center gap-4 mt-1 mb-3 text-[10px] text-gray-500">
					<span className="flex items-center gap-1">
						<span
							style={{
								display: "inline-block",
								width: 20,
								height: 2,
								background: "#f97316",
							}}
						/>
						hybr
					</span>
					{lmT != null && (
						<span className="flex items-center gap-1">
							<span
								style={{
									display: "inline-block",
									width: 20,
									height: 2,
									background: "#3b82f6",
								}}
							/>
							LM fallback
						</span>
					)}
					<span className="ml-auto text-gray-400">
						{history.length} samples
					</span>
				</div>

				{top5.length > 0 && (
					<div>
						<div className="text-[10px] font-semibold text-gray-500 uppercase tracking-wider mb-1">
							Worst residuals
						</div>
						<table className="w-full text-[10px] font-mono">
							<tbody>
								{top5.map((r) => (
									<tr key={r.name} className="border-t border-gray-100">
										<td className="py-0.5 pr-3 text-gray-600 truncate max-w-[320px]">
											{r.name}
										</td>
										<td className="py-0.5 text-right text-gray-800">
											{r.residual.toExponential(2)}
										</td>
									</tr>
								))}
							</tbody>
						</table>
					</div>
				)}
			</div>
		</div>
	);
};

export default ConvergenceChart;
