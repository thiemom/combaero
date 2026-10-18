import axios from "axios";
import { Download, Play, Zap } from "lucide-react";
import { useEffect, useRef, useState } from "react";
import { ReactFlowProvider } from "reactflow";
import { isBackendReachable, solveNetwork } from "./api";
import Inspector from "./components/Inspector";
import NetworkCanvas from "./components/NetworkCanvas";
import Sidebar from "./components/Sidebar";
import useStore from "./store/useStore";
import { readAutosave, startAutosave } from "./utils/autosave";
import { solveErrorMessage } from "./utils/solveError";
import { validateNetwork } from "./utils/validation";

const App = () => {
	const reactFlowWrapper = useRef<HTMLDivElement>(null);
	const {
		nodes,
		edges,
		solverSettings,
		setSolveResults,
		isExporting,
		exportNetworkResults,
		refreshPreview,
	} = useStore();
	const [isSolving, setIsSolving] = useState(false);
	// When the canvas came back from the autosave (#596), so the restored
	// network is never mistaken for a fresh start.
	const [restoredAt, setRestoredAt] = useState<string | null>(null);

	// Restore the last session once, then keep it saved. Only into an empty
	// canvas: StrictMode runs this twice in development, and a restore must
	// never overwrite work.
	useEffect(() => {
		const store = useStore.getState();
		if (store.nodes.length === 0 && store.edges.length === 0) {
			const rec = readAutosave();
			if (rec) {
				store.loadNetwork(rec);
				const t = new Date(rec.savedAt);
				setRestoredAt(
					Number.isNaN(t.getTime())
						? "earlier"
						: t.toLocaleString([], {
								dateStyle: "short",
								timeStyle: "short",
							}),
				);
			}
		}
		return startAutosave(useStore);
	}, []);

	// Keep the geometry-inheritance preview (#559) fresh as the graph changes,
	// debounced so a drag or a burst of edits triggers one resolve, not many.
	// The preview is advisory: the Inspector reads it for inherited placeholders
	// and to flag unresolved geometry before Solve. nodes/edges are intentional
	// re-run triggers -- refreshPreview reads the current graph from the store via
	// get(), so biome cannot see that the effect depends on them.
	// biome-ignore lint/correctness/useExhaustiveDependencies: nodes/edges are intentional re-run triggers; refreshPreview reads the graph via store get()
	useEffect(() => {
		const t = setTimeout(() => {
			void refreshPreview();
		}, 300);
		return () => clearTimeout(t);
	}, [nodes, edges, refreshPreview]);

	const handleSolve = async () => {
		const errors = validateNetwork(nodes);
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
			// A request error may mean the backend is down rather than the solve
			// failing; a /health probe tells the two apart.
			const backendReachable = axios.isAxiosError(err)
				? await isBackendReachable()
				: true;
			setSolveResults({
				success: false,
				...solveErrorMessage(err, backendReachable),
			});
		} finally {
			setIsSolving(false);
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
						onClick={exportNetworkResults}
						disabled={isExporting || nodes.length === 0}
						className="flex items-center gap-2 px-4 py-1.5 border rounded-md hover:bg-slate-50 transition-colors text-sm font-medium disabled:opacity-60 disabled:cursor-not-allowed min-w-[124px] justify-center"
					>
						{isExporting ? (
							<>
								<svg
									aria-label="Exporting"
									role="img"
									className="animate-spin text-slate-400"
									width={16}
									height={16}
									viewBox="0 0 24 24"
									fill="none"
									stroke="currentColor"
									strokeWidth={2.5}
								>
									<title>Exporting</title>
									<circle cx="12" cy="12" r="10" strokeOpacity={0.25} />
									<path d="M12 2a10 10 0 0 1 10 10" />
								</svg>
								Exporting…
							</>
						) : (
							<>
								<Download size={16} /> Export CSV
							</>
						)}
					</button>
					<button
						type="button"
						onClick={handleSolve}
						disabled={isSolving || nodes.length === 0}
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

			{restoredAt && (
				<div className="flex items-center justify-between gap-4 px-6 py-1 bg-amber-50 border-b border-amber-200 text-xs text-amber-800">
					<span>
						Restored your last session (autosaved {restoredAt}). Use Clear in
						the sidebar to start empty.
					</span>
					<button
						type="button"
						onClick={() => setRestoredAt(null)}
						className="px-2 rounded hover:bg-amber-100"
						aria-label="Dismiss restore note"
					>
						×
					</button>
				</div>
			)}

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
					: isExporting
						? "Engine Solver Status: Exporting Data…"
						: "Engine Solver Status: Idle"}{" "}
				| Nodes: {nodes.length} | Elements: {edges.length}
			</footer>
		</div>
	);
};

export default App;
