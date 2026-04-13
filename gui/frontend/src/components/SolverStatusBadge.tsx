import { CheckCircle, XCircle, X } from "lucide-react";
import type React from "react";
import useStore from "../store/useStore";

const SolverStatusBadge: React.FC = () => {
	const { solveResults, dismissSolveStatus } = useStore();

	if (!solveResults) return null;

	const isSuccess = solveResults.success;
	const message = solveResults.message || (isSuccess ? "Solve completed" : "Failed to converge");
	const finalNorm = solveResults.final_norm;

	return (
		<div
			className={`absolute top-4 left-1/2 -translate-x-1/2 z-50 flex items-center gap-3 px-4 py-2 rounded shadow-md border ${
				isSuccess ? "bg-green-50 border-green-200 text-green-800" : "bg-red-50 border-red-200 text-red-800"
			}`}
		>
			<div className="flex items-center gap-2">
				{isSuccess ? <CheckCircle size={20} /> : <XCircle size={20} />}
				<span className="font-semibold text-sm">
					{isSuccess ? "Solver: Success" : "Solver: Failed"}
				</span>
			</div>
			
			<div className="flex flex-col border-l pl-3 ml-1 border-current/20">
				{finalNorm !== undefined && finalNorm !== null && (
					<span className="text-xs font-mono">
						||F|| = {finalNorm.toExponential(3)}
					</span>
				)}
				<span className="text-xs opacity-80 max-w-xs truncate" title={message}>
					{message}
				</span>
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
	);
};

export default SolverStatusBadge;
