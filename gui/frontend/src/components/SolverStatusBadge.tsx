import { CheckCircle, X, XCircle } from "lucide-react";
import type React from "react";
import { useState } from "react";
import useStore from "../store/useStore";
import ConvergenceChart from "./ConvergenceChart";

const SolverStatusBadge: React.FC = () => {
	const { solveResults, dismissSolveStatus } = useStore();
	const [showChart, setShowChart] = useState(false);

	if (!solveResults) return null;

	const isSuccess = solveResults.success;
	const message =
		solveResults.message ||
		(isSuccess ? "Solve completed" : "Failed to converge");
	const finalNorm = solveResults.final_norm;
	const hasHistory = (solveResults.convergence_history?.length ?? 0) > 0;

	if (isSuccess) {
		return (
			<>
				<div className="absolute top-4 left-1/2 -translate-x-1/2 z-50 flex items-center gap-3 px-4 py-2 rounded shadow-md border bg-green-50 border-green-200 text-green-800">
					<div className="flex items-center gap-2">
						<CheckCircle size={20} />
						<span className="font-semibold text-sm">Solver: Success</span>
					</div>
					<div className="flex flex-col border-l pl-3 ml-1 border-green-300">
						{finalNorm !== undefined && finalNorm !== null && (
							<span className="text-xs font-mono">
								||F|| = {finalNorm.toExponential(3)}
							</span>
						)}
						<span className="text-xs opacity-80">{message}</span>
						{hasHistory && (
							<button
								type="button"
								onClick={() => setShowChart(true)}
								className="text-[10px] text-green-600 underline hover:text-green-800 text-left mt-0.5"
							>
								View convergence
							</button>
						)}
					</div>
					<button
						type="button"
						onClick={dismissSolveStatus}
						className="ml-2 hover:bg-black/5 p-1 rounded-full transition-colors"
						aria-label="Dismiss"
					>
						<X size={16} />
					</button>
				</div>
				{showChart && (
					<ConvergenceChart
						history={solveResults.convergence_history ?? []}
						worstResiduals={solveResults.worst_residuals}
						lmStartedAt={solveResults.lm_started_at_eval}
						onClose={() => setShowChart(false)}
					/>
				)}
			</>
		);
	}

	return (
		<>
			<div className="absolute top-4 left-1/2 -translate-x-1/2 z-50 w-[480px] max-w-[90vw] rounded shadow-lg border bg-red-50 border-red-200 text-red-800">
				<div className="flex items-center gap-2 px-4 py-2 border-b border-red-200">
					<XCircle size={20} className="shrink-0" />
					<span className="font-semibold text-sm flex-1">Solver: Failed</span>
					{finalNorm !== undefined && finalNorm !== null && (
						<span className="text-xs font-mono opacity-70 mr-2">
							||F|| = {finalNorm.toExponential(3)}
						</span>
					)}
					<button
						type="button"
						onClick={dismissSolveStatus}
						className="hover:bg-red-100 p-1 rounded-full transition-colors shrink-0"
						aria-label="Dismiss"
					>
						<X size={16} />
					</button>
				</div>
				<div className="px-4 py-2 max-h-40 overflow-y-auto">
					<pre className="text-xs font-mono whitespace-pre-wrap break-words">
						{message}
					</pre>
				</div>
				{hasHistory && (
					<div className="px-4 pb-2">
						<button
							type="button"
							onClick={() => setShowChart(true)}
							className="text-[10px] text-red-500 underline hover:text-red-700"
						>
							View convergence
						</button>
					</div>
				)}
			</div>
			{showChart && (
				<ConvergenceChart
					history={solveResults.convergence_history ?? []}
					worstResiduals={solveResults.worst_residuals}
					lmStartedAt={solveResults.lm_started_at_eval}
					onClose={() => setShowChart(false)}
				/>
			)}
		</>
	);
};

export default SolverStatusBadge;
