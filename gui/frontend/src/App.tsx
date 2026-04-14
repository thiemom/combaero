import axios from "axios";
import { Download, Play, Zap } from "lucide-react";
import { useRef, useState } from "react";
import { ReactFlowProvider } from "reactflow";
import { exportResults, solveNetwork } from "./api";
import Inspector from "./components/Inspector";
import NetworkCanvas from "./components/NetworkCanvas";
import Sidebar from "./components/Sidebar";
import useStore from "./store/useStore";

const App = () => {
	const reactFlowWrapper = useRef<HTMLDivElement>(null);
	const { nodes, edges, solverSettings, setSolveResults } = useStore();
	const [isSolving, setIsSolving] = useState(false);

	const validateNetwork = () => {
		const errors: string[] = [];
		for (const node of nodes) {
			if (node.type === "pressure_boundary" || node.type === "mass_boundary") {
				if (
					node.type === "pressure_boundary" &&
					(node.data.P_total || 0) <= 0
				) {
					errors.push(`${node.id}: Total Pressure must be > 0`);
				}
				if ((node.data.T_total || 0) <= 0) {
					errors.push(`${node.id}: Temperature must be > 0`);
				}
				if (node.data.composition?.source === "custom") {
					const values = Object.values(
						node.data.composition.custom_fractions || {},
					) as number[];
					const sum = values.reduce((a: number, b: number) => a + (b || 0), 0);
					if (Math.abs(sum - 1.0) > 1e-3) {
						errors.push(
							`${node.id}: Custom composition sum is ${sum.toFixed(3)} (expected 1.0)`,
						);
					}
				}
			}
			if (node.type === "channel") {
				if ((node.data.D || 0) <= 0)
					errors.push(`${node.id}: Channel diameter must be > 0`);
				if ((node.data.L || 0) <= 0)
					errors.push(`${node.id}: Channel length must be > 0`);
			}
			if (node.type === "orifice") {
				if ((node.data.area || 0) <= 0)
					errors.push(`${node.id}: Orifice area must be > 0`);
			}
		}
		return errors;
	};

	const handleSolve = async () => {
		const errors = validateNetwork();
		if (errors.length > 0) {
			alert(
				`Please fix the following errors before solving:\n\n${errors.join("\n")}`,
			);
			return;
		}

		setSolveResults(null);
		setIsSolving(true);
		try {
			const results = await solveNetwork({
				nodes,
				edges,
				solver_settings: solverSettings,
			});
			setSolveResults(results);
		} catch (err) {
			console.error("Solve failed:", err);
			if (axios.isAxiosError(err)) {
				const detail = err.response?.data?.detail;
				alert(
					detail
						? `Solve failed: ${detail}`
						: "Solve failed. Check console for details.",
				);
				return;
			}
			alert("Solve failed. Check console for details.");
		} finally {
			setIsSolving(false);
		}
	};

	const handleExport = async () => {
		try {
			await exportResults({ nodes, edges, solver_settings: solverSettings });
		} catch (err) {
			console.error("Export failed:", err);
		}
	};

	return (
		<div className="flex flex-col h-screen w-screen overflow-hidden text-slate-900">
			{/* Header */}
			<header className="h-14 border-b bg-white flex items-center justify-between px-6 shadow-sm z-10">
				<div className="flex items-center gap-2">
					<Zap className="text-orange-500 fill-orange-500" size={24} />
					<h1 className="text-xl font-bold tracking-tight">
						CombAero{" "}
						<span className="font-light text-slate-400">Network Designer</span>
					</h1>
				</div>

				<div className="flex items-center gap-3">
					<button
						type="button"
						onClick={handleExport}
						className="flex items-center gap-2 px-4 py-1.5 border rounded-md hover:bg-slate-50 transition-colors text-sm font-medium"
					>
						<Download size={16} /> Export CSV
					</button>
					<button
						type="button"
						onClick={handleSolve}
						disabled={isSolving}
						className="flex items-center gap-2 px-6 py-1.5 bg-orange-500 text-white rounded-md hover:bg-orange-600 transition-colors shadow-sm text-sm font-bold disabled:opacity-60 disabled:cursor-not-allowed"
					>
						{isSolving ? (
							<>
								<svg
									aria-label="Solving"
									role="img"
									className="animate-spin"
									width={16}
									height={16}
									viewBox="0 0 24 24"
									fill="none"
									stroke="currentColor"
									strokeWidth={2.5}
								>
									<title>Solving</title>
									<circle cx="12" cy="12" r="10" strokeOpacity={0.25} />
									<path d="M12 2a10 10 0 0 1 10 10" />
								</svg>
								Solving…
							</>
						) : (
							<>
								<Play size={16} fill="white" /> Solve Network
							</>
						)}
					</button>
				</div>
			</header>

			{/* Main Content */}
			<div className="flex flex-grow overflow-hidden" ref={reactFlowWrapper}>
				<Sidebar />
				<ReactFlowProvider>
					<NetworkCanvas />
				</ReactFlowProvider>
				<Inspector />
			</div>

			{/* Footer / Status Bar */}
			<footer className="h-8 border-t bg-stone-100 flex items-center px-4 text-[10px] text-stone-500 uppercase tracking-widest font-bold">
				{isSolving
					? "Engine Solver Status: Running…"
					: "Engine Solver Status: Idle"}{" "}
				| Nodes: {nodes.length} | Elements: {edges.length}
			</footer>
		</div>
	);
};

export default App;
