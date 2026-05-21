import {
	Activity,
	ArrowRight,
	ChevronRight,
	Crosshair,
	Database,
	Flame,
	Link,
	LogIn,
	LogOut,
	Maximize2,
	Save,
	Split,
	Square,
	Tornado,
	Upload,
	Wind,
} from "lucide-react";
import type React from "react";
import { useEffect, useRef, useState } from "react";
import { checkContinuationAvailable } from "../api";
import useStore from "../store/useStore";
import { DiagnosticsPanel } from "./DiagnosticsPanel";
import NumericInput from "./NumericInput";

const Sidebar = () => {
	const fileInputRef = useRef<HTMLInputElement>(null);
	const [filename, setFilename] = useState("network");

	const {
		saveNetwork,
		loadNetwork,
		solverSettings,
		updateSolverSettings,
		unitPreferences,
		setRotSpeedUnit,
		nodes,
		setNodes,
		solveResults,
	} = useStore();

	const [isContinuationAvailable, setIsContinuationAvailable] = useState(false);

	useEffect(() => {
		// Use dependencies to satisfy linter while ensuring re-check on change
		const _trigger = [solverSettings, solveResults];
		void _trigger;
		const check = async () => {
			try {
				const available = await checkContinuationAvailable();
				setIsContinuationAvailable(available);
			} catch (err) {
				console.error("Failed to check continuation availability:", err);
			}
		};
		check();
	}, [solverSettings, solveResults]);

	const resetAllInitialGuesses = () => {
		let touched = 0;
		const next = nodes.map((n) => {
			if (
				n.data?.initial_guess &&
				Object.keys(n.data.initial_guess).length > 0
			) {
				touched += 1;
				const { initial_guess: _ig, ...rest } = n.data;
				return { ...n, data: { ...rest, initial_guess: {} } };
			}
			return n;
		});
		if (touched === 0) return;
		setNodes(next);
	};

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

			<div className="flex flex-col gap-1">
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
				onDragStart={(event) => onDragStart(event, "wall")}
				draggable
			>
				<Square size={18} className="text-stone-600 fill-stone-600" />
				<span>Wall</span>
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
			<div className="text-[9px] text-stone-400 flex gap-2">
				<span>
					<span className="font-bold text-blue-500">i</span> flow in
				</span>
				<span>
					<span className="font-bold text-green-500">o</span> flow out
				</span>
				<span>
					<span className="font-bold text-orange-400">q</span> thermal
				</span>
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
				onDragStart={(event) => onDragStart(event, "area_change")}
				draggable
			>
				<Maximize2 size={18} className="text-blue-500" />
				<span>Area Change</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "discrete_loss")}
				draggable
			>
				<Activity size={18} className="text-purple-500" />
				<span>Discrete Loss</span>
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

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "tee_junction")}
				draggable
			>
				<Split size={18} className="text-teal-500" />
				<span>Tee Junction</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "vortex")}
				draggable
			>
				<Tornado size={18} className="text-violet-500" />
				<span>Vortex</span>
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

			{/* Display field picker */}
			<div className="border-t border-stone-100 pt-3 mt-1">
				<DiagnosticsPanel />
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
						<option
							value="continuation"
							disabled={!isContinuationAvailable}
							title={
								!isContinuationAvailable
									? "Run a standard solve first to enable continuation"
									: ""
							}
						>
							Continuation {!isContinuationAvailable ? "(Locked)" : ""}
						</option>
					</select>
				</div>
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						Solver Method
					</label>
					<select
						className="p-1 border rounded text-xs bg-white"
						value={solverSettings.method}
						onChange={(e) => updateSolverSettings({ method: e.target.value })}
					>
						<option value="hybr">hybr (Powell)</option>
						<option value="lm">lm (Levenberg-M.)</option>
						<option value="krylov">krylov (Iterative)</option>
						<option value="broyden1">broyden1 (Quasi-N.)</option>
						<option value="broyden2">broyden2 (Quasi-N.)</option>
						<option value="anderson">anderson</option>
						<option value="linearmixing">linearmixing</option>
						<option value="diagbroyden">diagbroyden</option>
						<option value="excitingmixing">excitingmixing</option>
						<option value="df-sane">df-sane</option>
					</select>
				</div>
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						Max Solver Time (s)
					</label>
					<NumericInput
						className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200"
						value={solverSettings.timeout}
						onChange={(val) => updateSolverSettings({ timeout: val })}
						onClear={() => updateSolverSettings({ timeout: null })}
						min={0}
						placeholder="None"
					/>
				</div>

				<div className="border-t border-stone-200 pt-2 flex flex-col gap-2">
					<button
						type="button"
						onClick={resetAllInitialGuesses}
						className="w-full text-[10px] font-bold uppercase tracking-wider bg-white hover:bg-amber-50 text-amber-700 border border-amber-200 rounded px-2 py-1.5 transition-colors"
						title="Clear initial-guess overrides on every node and element; solver reverts to auto-seeding"
					>
						Reset All Initial Guesses
					</button>

					{solveResults && (solveResults as any).isTopologyMismatch && (
						<button
							type="button"
							onClick={() => updateSolverSettings({ init_strategy: "default" })}
							className="w-full text-[10px] font-bold uppercase tracking-wider bg-red-50 hover:bg-red-100 text-red-700 border border-red-200 rounded px-2 py-1.5 transition-colors"
						>
							Reset Solver (Cold Start)
						</button>
					)}
				</div>
			</div>

			<div className="mt-6 text-sm font-bold text-gray-500 uppercase tracking-wider">
				Global Settings
			</div>
			<div className="flex flex-col gap-4 mt-2 p-3 bg-stone-50 rounded border border-stone-200">
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						Shaft Speed
					</label>
					<div className="flex items-center gap-1 w-full">
						<NumericInput
							className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200 flex-1 min-w-0"
							value={
								solverSettings.omega_rpm != null
									? solverSettings.omega_rpm *
										(unitPreferences.rotSpeed === "Hz" ? 1 / 60 : 1)
									: null
							}
							onChange={(val) =>
								updateSolverSettings({
									omega_rpm: val * (unitPreferences.rotSpeed === "Hz" ? 60 : 1),
								})
							}
							onClear={() => updateSolverSettings({ omega_rpm: null })}
							min={0}
							placeholder="None"
						/>
						<select
							value={unitPreferences.rotSpeed}
							onChange={(e) => setRotSpeedUnit(e.target.value as "rpm" | "Hz")}
							className="w-12 h-7 border rounded text-[9px] bg-white outline-none shrink-0"
						>
							<option value="rpm">rpm</option>
							<option value="Hz">Hz</option>
						</select>
					</div>
				</div>
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						Nu Multiplier
					</label>
					<NumericInput
						className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200"
						value={solverSettings.Nu_multiplier ?? null}
						onChange={(val) => updateSolverSettings({ Nu_multiplier: val })}
						onClear={() => updateSolverSettings({ Nu_multiplier: null })}
						min={0}
						placeholder="1.0"
					/>
					<span className="text-[9px] text-stone-400">
						Scales heat transfer on all elements multiplicatively
					</span>
				</div>
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						f Multiplier
					</label>
					<NumericInput
						className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200"
						value={solverSettings.f_multiplier ?? null}
						onChange={(val) => updateSolverSettings({ f_multiplier: val })}
						onClear={() => updateSolverSettings({ f_multiplier: null })}
						min={0}
						placeholder="1.0"
					/>
					<span className="text-[9px] text-stone-400">
						Scales friction factor on all elements multiplicatively
					</span>
				</div>
				<div className="border-t border-stone-200 my-1" />
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						Ambient Temperature
					</label>
					<div className="flex items-center gap-1 w-full">
						<NumericInput
							className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200 flex-1 min-w-0"
							value={solverSettings.ambient_T ?? null}
							onChange={(val) => updateSolverSettings({ ambient_T: val })}
							onClear={() => updateSolverSettings({ ambient_T: null })}
							placeholder="288.15"
						/>
						<span className="text-[10px] text-stone-400 w-10 text-right shrink-0">
							K
						</span>
					</div>
				</div>
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						Ambient Pressure
					</label>
					<div className="flex items-center gap-1 w-full">
						<NumericInput
							className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200 flex-1 min-w-0"
							value={
								solverSettings.ambient_P != null
									? solverSettings.ambient_P *
										(unitPreferences.pressure === "kPa"
											? 1 / 1000
											: unitPreferences.pressure === "MPa"
												? 1 / 1e6
												: 1)
									: null
							}
							onChange={(val) =>
								updateSolverSettings({
									ambient_P:
										val *
										(unitPreferences.pressure === "kPa"
											? 1000
											: unitPreferences.pressure === "MPa"
												? 1e6
												: 1),
								})
							}
							onClear={() => updateSolverSettings({ ambient_P: null })}
							placeholder={
								unitPreferences.pressure === "kPa"
									? "101.325"
									: unitPreferences.pressure === "MPa"
										? "0.101325"
										: "101325"
							}
						/>
						<span className="text-[10px] text-stone-400 w-10 text-right shrink-0">
							{unitPreferences.pressure}
						</span>
					</div>
				</div>
				<div className="flex flex-col gap-1">
					<label className="text-[10px] font-bold text-gray-400 uppercase">
						Ambient Rel. Humidity
					</label>
					<div className="flex items-center gap-1 w-full">
						<NumericInput
							className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200 flex-1 min-w-0"
							value={solverSettings.ambient_RH ?? null}
							onChange={(val) => updateSolverSettings({ ambient_RH: val })}
							onClear={() => updateSolverSettings({ ambient_RH: null })}
							min={0}
							max={1}
							placeholder="0.6"
						/>
						<span className="text-[10px] text-stone-400 w-10 text-right shrink-0">
							0–1
						</span>
					</div>
				</div>
			</div>
		</aside>
	);
};

export default Sidebar;
