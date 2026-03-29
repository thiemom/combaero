import axios from "axios";
import { Download, Play, Zap } from "lucide-react";
import type React from "react";
import { useCallback, useRef } from "react";
import { ReactFlowProvider } from "reactflow";
import { exportResults, solveNetwork } from "./api";
import Inspector from "./components/Inspector";
import NetworkCanvas from "./components/NetworkCanvas";
import Sidebar from "./components/Sidebar";
import useStore from "./store/useStore";

const App = () => {
	const reactFlowWrapper = useRef<HTMLDivElement>(null);
	const { nodes, edges, setNodes, setSolveResults } = useStore();

	const onDragOver = useCallback((event: React.DragEvent) => {
		event.preventDefault();
		event.dataTransfer.dropEffect = "move";
	}, []);

	const onDrop = useCallback(
		(event: React.DragEvent) => {
			event.preventDefault();

			const type = event.dataTransfer.getData("application/reactflow");

			if (typeof type === "undefined" || !type) {
				return;
			}

			const position = { x: event.clientX - 300, y: event.clientY - 100 };
			const id = `node_${Date.now()}`;
			let data = { id };

			if (type === "mass_boundary") {
				data = { ...data, m_dot: 1.0, T_total: 300 } as any;
			} else if (type === "pressure_boundary") {
				data = { ...data, P_total: 101325 } as any;
			} else if (type === "pipe") {
				data = { ...data, L: 1.0, D: 0.1, roughness: 1e-5 } as any;
			} else if (type === "orifice") {
				data = { ...data, area: 0.01, Cd: 0.6 } as any;
			}

			const newNode = {
				id,
				type,
				position,
				data,
			};

			setNodes(nodes.concat(newNode));
		},
		[nodes, setNodes],
	);

	const handleSolve = async () => {
		try {
			const results = await solveNetwork({ nodes, edges });
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
		}
	};

	const handleExport = async () => {
		try {
			await exportResults({ nodes, edges });
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
						className="flex items-center gap-2 px-6 py-1.5 bg-orange-500 text-white rounded-md hover:bg-orange-600 transition-colors shadow-sm text-sm font-bold"
					>
						<Play size={16} fill="white" /> Solve Network
					</button>
				</div>
			</header>

			{/* Main Content */}
			<div className="flex flex-grow overflow-hidden" ref={reactFlowWrapper}>
				<Sidebar />
				<ReactFlowProvider>
					<section
						className="flex-grow h-full"
						onDragOver={onDragOver}
						onDrop={onDrop}
						aria-label="Flow Network Canvas"
					>
						<NetworkCanvas />
					</section>
				</ReactFlowProvider>
				<Inspector />
			</div>

			{/* Footer / Status Bar */}
			<footer className="h-8 border-t bg-stone-100 flex items-center px-4 text-[10px] text-stone-500 uppercase tracking-widest font-bold">
				Engine Solver Status: Idle | Nodes: {nodes.length} | Elements:{" "}
				{edges.length}
			</footer>
		</div>
	);
};

export default App;
