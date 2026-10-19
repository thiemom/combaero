import { CheckCircle, X, XCircle } from "lucide-react";
import type React from "react";
import { useState } from "react";
import useStore from "../store/useStore";
import { summary } from "../utils/incompressible";
import { balanceIsOpen, type Ledger } from "../utils/secondLaw";
import BalanceModal from "./BalanceModal";
import ConvergenceChart from "./ConvergenceChart";

/** One line on the second law and the balance (#545). A failing element or an open balance also fails the solve (#598, #552). */
const BalanceLine: React.FC<{ onOpen: () => void }> = ({ onOpen }) => {
	const { solveResults } = useStore();
	const ledger = solveResults?.ledger as Ledger | undefined;
	if (!ledger || ledger.m_dot_in === undefined) return null;
	const n = solveResults?.second_law_failures?.length ?? 0;
	const open = balanceIsOpen(ledger);
	const warn = n > 0 || open;
	const text = [
		n > 0 ? `${n} element${n > 1 ? "s" : ""} create total pressure` : null,
		open ? "balance does not close" : null,
	]
		.filter(Boolean)
		.join("; ");
	return (
		<button
			type="button"
			onClick={onOpen}
			className={`text-[10px] underline text-left mt-0.5 ${warn ? "text-amber-700 hover:text-amber-900" : "opacity-70 hover:opacity-100"}`}
		>
			{warn ? `Second law: ${text}` : "View balance"}
		</button>
	);
};

const SolverStatusBadge: React.FC = () => {
	const { solveResults, dismissSolveStatus, solverSettings } = useStore();
	const [showChart, setShowChart] = useState(false);
	const [showBalance, setShowBalance] = useState(false);

	if (!solveResults) return null;

	const isSuccess = solveResults.success;
	const incompressibleNote = summary(solveResults.incompressible_warnings);
	const message =
		solveResults.message ||
		(isSuccess ? "Solve completed" : "Failed to converge");
	const finalNorm = solveResults.final_norm;
	const hasHistory = (solveResults.convergence_history?.length ?? 0) > 0;
	// The solver may auto-upgrade "default" (e.g. to analytical_pt_prop on
	// MPCE networks); show what actually ran, flagged when it differs from
	// the requested strategy.
	const usedInit = solveResults.solver_settings_used?.init_strategy as
		| string
		| undefined;
	const initLabel =
		usedInit &&
		`init: ${usedInit}${usedInit !== solverSettings.init_strategy ? " (auto)" : ""}`;

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
						{initLabel && (
							<span className="text-[10px] font-mono opacity-70">
								{initLabel}
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
						<BalanceLine onOpen={() => setShowBalance(true)} />
						{incompressibleNote && (
							<span
								className="text-[10px] text-amber-700 mt-0.5"
								title="Report only: the solution is unchanged. Above this Mach the incompressible closures drift from the compressible result (about 4% on a channel at Ma 0.2). A rule of thumb, not a guarantee -- the application stays your call."
							>
								{incompressibleNote}
							</span>
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
				{showBalance && solveResults.ledger?.m_dot_in !== undefined && (
					<BalanceModal
						ledger={solveResults.ledger as Ledger}
						failures={solveResults.second_law_failures ?? []}
						onClose={() => setShowBalance(false)}
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
					{initLabel && (
						<div className="text-[10px] font-mono opacity-70 pb-1">
							{initLabel}
						</div>
					)}
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
				<div className="px-4 pb-2">
					<BalanceLine onOpen={() => setShowBalance(true)} />
				</div>
			</div>
			{showChart && (
				<ConvergenceChart
					history={solveResults.convergence_history ?? []}
					worstResiduals={solveResults.worst_residuals}
					lmStartedAt={solveResults.lm_started_at_eval}
					onClose={() => setShowChart(false)}
				/>
			)}
			{showBalance && solveResults.ledger?.m_dot_in !== undefined && (
				<BalanceModal
					ledger={solveResults.ledger as Ledger}
					failures={solveResults.second_law_failures ?? []}
					onClose={() => setShowBalance(false)}
				/>
			)}
		</>
	);
};

export default SolverStatusBadge;
