import {
	ArrowRight,
	ChevronRight,
	Crosshair,
	Database,
	Flame,
	Link,
	LogIn,
	LogOut,
	Save,
	Upload,
	Wind,
} from "lucide-react";
import type React from "react";
import { useRef, useState } from "react";
import useStore from "../store/useStore";

const Sidebar = () => {
	const fileInputRef = useRef<HTMLInputElement>(null);
	const [filename, setFilename] = useState("network");

	const { saveNetwork, loadNetwork, solverSettings, updateSolverSettings } =
		useStore();

	const onSave = () => saveNetwork(filename.trim() || "network");

	const onDragStart = (event: React.DragEvent, nodeType: string) => {
		event.dataTransfer.setData("application/reactflow", nodeType);
		event.dataTransfer.effectAllowed = "move";
	};

	const handleBlur = () => {
		if (!filename.trim()) {
			setFilename("network");
		}
	};

	const handleFileChange = (event: React.ChangeEvent<HTMLInputElement>) => {
		const file = event.target.files?.[0];
		if (file) {
			const reader = new FileReader();
			reader.onload = (e) => {
				try {
					const json = JSON.parse(e.target?.result as string);
					loadNetwork(json);
				} catch (error) {
					console.error("Failed to parse network JSON:", error);
					alert("Invalid network file");
				}
			};
			reader.readAsText(file);
		}
	};

	return (
		<aside className="w-64 border-r bg-white p-4 flex flex-col gap-4 overflow-y-auto">
			<div className="text-sm font-bold text-gray-500 uppercase tracking-wider">
				Project
			</div>

			<div className="flex gap-2">
				<button
					type="button"
					onClick={onSave}
					className="flex-1 flex items-center justify-center gap-2 p-2 border rounded hover:bg-stone-50 transition-colors text-sm bg-white"
				>
					<Save size={16} className="text-blue-600" />
					<span>Save</span>
				</button>
				<button
					type="button"
					onClick={() => fileInputRef.current?.click()}
					className="flex-1 flex items-center justify-center gap-2 p-2 border rounded hover:bg-stone-50 transition-colors text-sm bg-white"
				>
					<Upload size={16} className="text-emerald-600" />
					<span>Load</span>
				</button>
				<input
					type="file"
					ref={fileInputRef}
					onChange={handleFileChange}
					className="hidden"
					accept=".json"
				/>
			</div>

			<div className="mt-4 text-sm font-bold text-gray-500 uppercase tracking-wider">
				Nodes
			</div>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "mass_boundary")}
				draggable
			>
				<LogIn size={18} className="text-orange-500" />
				<span>Mass Boundary</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "plenum")}
				draggable
			>
				<Database size={18} className="text-stone-500" />
				<span>Plenum</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "pressure_boundary")}
				draggable
			>
				<LogOut size={18} className="text-blue-500" />
				<span>Pressure Boundary</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "combustor")}
				draggable
			>
				<Flame size={18} className="text-red-500" />
				<span>Combustor</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "momentum_chamber")}
				draggable
			>
				<Wind size={18} className="text-indigo-500" />
				<span>Momentum Chamber</span>
			</button>

			<div className="mt-4 text-sm font-bold text-gray-500 uppercase tracking-wider">
				Elements
			</div>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "channel")}
				draggable
			>
				<ArrowRight size={18} className="text-stone-500" />
				<span>Channel</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "orifice")}
				draggable
			>
				<ChevronRight size={18} className="text-orange-400" />
				<span>Orifice</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "lossless_connection")}
				draggable
			>
				<Link size={18} className="text-slate-400" />
				<span>Lossless Connection</span>
			</button>

			{/* Diagnostics divider */}
			<div className="text-[10px] font-bold text-stone-400 uppercase tracking-wider mt-2 mb-0.5">
				Diagnostics
			</div>
			<button
				type="button"
				className="flex items-center gap-2 p-2 border border-blue-200 rounded cursor-grab hover:bg-blue-50 transition-colors w-full text-left bg-white text-blue-700"
				onDragStart={(event) => onDragStart(event, "probe")}
				draggable
			>
				<Crosshair size={18} className="text-blue-400" />
				<span>Probe</span>
			</button>

			<div className="mt-auto pt-6 border-t border-stone-200">
				<div className="flex flex-col gap-2 mb-4">
					<label className="text-[10px] font-bold text-gray-400 uppercase tracking-tighter">
						Filename
					</label>
					<div className="flex gap-1">
						<input
							type="text"
							value={filename}
							onChange={(e) => setFilename(e.target.value)}
							onBlur={handleBlur}
							className="p-1 border rounded text-xs w-full bg-white outline-none focus:ring-1 focus:ring-stone-200"
							placeholder="MyNetwork"
						/>
						<span className="text-gray-400 text-xs self-center">.json</span>
					</div>
				</div>
				<button
					type="button"
					onClick={onSave}
					className="flex items-center justify-center gap-2 p-2 bg-stone-800 text-white rounded hover:bg-stone-700 transition-colors w-full"
				>
					<Save size={18} />
					<span>Save Network</span>
				</button>
			</div>

			<div className="mt-6 text-sm font-bold text-gray-500 uppercase tracking-wider">
				Solver Configuration
			</div>
			<div className="flex flex-col gap-4 mt-2 p-3 bg-stone-50 rounded border border-stone-200">
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						Global Regime
					</label>
					<select
						className="p-1 border rounded text-xs bg-white"
						value={solverSettings.global_regime}
						onChange={(e) =>
							updateSolverSettings({ global_regime: e.target.value as any })
						}
					>
						<option value="compressible">Compressible (Full)</option>
						<option value="incompressible">Incompressible (Fast)</option>
					</select>
				</div>
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						Init Strategy
					</label>
					<select
						className="p-1 border rounded text-xs bg-white"
						value={solverSettings.init_strategy}
						onChange={(e) =>
							updateSolverSettings({ init_strategy: e.target.value as any })
						}
					>
						<option value="default">Default Cold Start</option>
						<option value="incompressible_warmstart">Incomp. Warmstart</option>
						<option value="homotopy">Load-Stepping (Homotopy)</option>
					</select>
				</div>
			</div>
		</aside>
	);
};

export default Sidebar;
