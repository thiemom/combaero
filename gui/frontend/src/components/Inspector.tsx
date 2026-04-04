import { Activity, Crosshair, RotateCw } from "lucide-react";
import useStore from "../store/useStore";
import CompositionEditor from "./CompositionEditor";
import NumericInput from "./NumericInput";
import UnitInput from "./UnitInput";

const Inspector = () => {
	const {
		nodes,
		edges,
		updateNodeData,
		updateEdgeData,
		unitPreferences,
		speciesMetadata,
	} = useStore();

	const handleExport = async () => {
		try {
			const res = await fetch("http://localhost:8000/export", {
				method: "POST",
				headers: { "Content-Type": "application/json" },
				body: JSON.stringify({ nodes, edges }),
			});
			const blob = await res.blob();
			const url = window.URL.createObjectURL(blob);
			const a = document.createElement("a");
			a.href = url;
			a.download = "results.csv";
			a.click();
		} catch (error) {
			console.error("Export failed:", error);
		}
	};

	const validateNetwork = () => {
		const errors: string[] = [];
		for (const node of nodes) {
			if (node.type === "pressure_boundary" || node.type === "mass_boundary") {
				if (
					(node.data.P_total || 0) <= 0 &&
					node.type === "pressure_boundary"
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
			if (node.type === "pipe") {
				if ((node.data.D || 0) <= 0)
					errors.push(`${node.id}: Diameter must be > 0`);
				if ((node.data.L || 0) <= 0)
					errors.push(`${node.id}: Length must be > 0`);
			}
			if (node.type === "orifice") {
				if ((node.data.area || 0) <= 0)
					errors.push(`${node.id}: Area must be > 0`);
			}
		}
		return errors;
	};

	const validationErrors = validateNetwork();

	// Get selected node
	const selectedNode = nodes.find((n) => n.selected);
	const selectedEdge = edges.find((e) => e.selected);

	if (!selectedNode && !selectedEdge) {
		return (
			<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 italic text-gray-400">
				Select a node or edge to edit properties
			</aside>
		);
	}

	// ── Probe Inspector ──────────────────────────────────────────
	if (selectedNode?.type === "probe") {
		const probeData = selectedNode.data;
		const targetNode = nodes.find((n) => n.id === probeData.target_id);
		const result = targetNode?.data?.result;

		// Build available quantity list from target's result keys
		const availableKeys: string[] = [];
		if (result) {
			if (result.state) {
				for (const k of Object.keys(result.state)) {
					if (typeof result.state[k] === "number")
						availableKeys.push(`state.${k}`);
				}
			}
			for (const k of Object.keys(result)) {
				if (k !== "state" && typeof result[k] === "number")
					availableKeys.push(k);
			}
		}

		const slotSelect = (
			slot: "slot1_key" | "slot2_key",
			current: string | undefined,
		) => (
			<select
				className="p-1.5 border rounded text-xs bg-white w-full"
				value={current ?? ""}
				onChange={(e) =>
					updateNodeData(selectedNode.id, {
						[slot]: e.target.value || undefined,
					})
				}
			>
				<option value="">— select quantity —</option>
				{availableKeys.map((k) => (
					<option key={k} value={k}>
						{k}
					</option>
				))}
			</select>
		);

		return (
			<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto">
				<h2 className="text-lg font-bold border-b pb-2 flex items-center gap-2 text-blue-600">
					<Crosshair size={18} />
					<span>Probe</span>
				</h2>

				{/* Label */}
				<div className="flex flex-col gap-1">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Label
					</label>
					<input
						type="text"
						value={probeData.label ?? "Probe"}
						onChange={(e) =>
							updateNodeData(selectedNode.id, { label: e.target.value })
						}
						className="p-1.5 border rounded text-xs w-full"
					/>
				</div>

				{/* Target */}
				<div className="flex flex-col gap-1">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Target Node / Element
					</label>
					<select
						className="p-1.5 border rounded text-xs bg-white w-full"
						value={probeData.target_id ?? ""}
						onChange={(e) =>
							updateNodeData(selectedNode.id, {
								target_id: e.target.value || undefined,
							})
						}
					>
						<option value="">— select target —</option>
						{nodes
							.filter((n) => n.id !== selectedNode.id && n.type !== "probe")
							.map((n) => (
								<option key={n.id} value={n.id}>
									{n.data?.label || n.type?.replace(/_/g, " ") || n.id}
								</option>
							))}
					</select>
					{probeData.target_id && !targetNode && (
						<span className="text-amber-500 text-[10px]">
							⚠ Target not found
						</span>
					)}
				</div>

				{/* Slot selectors */}
				<div className="flex flex-col gap-3 border-t pt-3">
					<div className="text-[10px] font-bold text-stone-400 uppercase">
						Display Slots
					</div>
					{!probeData.target_id && (
						<p className="text-[10px] text-stone-400 italic">
							Select a target first to see available quantities.
						</p>
					)}
					{probeData.target_id && (
						<>
							<div className="flex flex-col gap-1">
								<label className="text-[10px] text-stone-500">Slot 1</label>
								{slotSelect("slot1_key", probeData.slot1_key)}
							</div>
							<div className="flex flex-col gap-1">
								<label className="text-[10px] text-stone-500">Slot 2</label>
								{slotSelect("slot2_key", probeData.slot2_key)}
							</div>
						</>
					)}
					{probeData.target_id && availableKeys.length === 0 && (
						<p className="text-[10px] text-stone-400 italic">
							Run Solve to populate available quantities.
						</p>
					)}
				</div>
			</aside>
		);
	}

	if (selectedNode) {
		return (
			<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto">
				{validationErrors.length > 0 && (
					<div className="bg-amber-50 border border-amber-200 p-2 rounded text-[10px] text-amber-800 flex flex-col gap-1">
						<div className="font-bold uppercase">Network Warnings</div>
						<ul className="list-disc pl-3">
							{validationErrors.map((err, i) => (
								<li key={i}>{err}</li>
							))}
						</ul>
					</div>
				)}
				<h2 className="text-lg font-bold border-b pb-2 capitalize flex justify-between items-center">
					<span>{(selectedNode.type || "unknown").replace("_", " ")} Node</span>
					<div className="flex gap-2">
						<button
							type="button"
							onClick={() => {
								const currentRotation = selectedNode.data.rotation || 0;
								updateNodeData(selectedNode.id, {
									rotation: (currentRotation + 90) % 360,
								});
							}}
							className="flex items-center gap-1 text-[10px] bg-blue-50 hover:bg-blue-100 text-blue-600 px-2 py-1 rounded border border-blue-200 transition-colors uppercase font-bold"
							title="Rotate node 90° clockwise"
						>
							<RotateCw size={12} />
							Rotate
						</button>
						<button
							type="button"
							onClick={handleExport}
							className="text-[10px] bg-stone-100 hover:bg-stone-200 px-2 py-1 rounded border border-stone-300 transition-colors uppercase font-bold"
						>
							Export CSV
						</button>
					</div>
				</h2>

				<div className="flex flex-col gap-2 mb-4">
					<label
						htmlFor={`node_label_${selectedNode.id}`}
						className="text-xs font-bold text-stone-500 uppercase"
					>
						Label / Name
					</label>
					<input
						id={`node_label_${selectedNode.id}`}
						type="text"
						value={selectedNode.data.label || ""}
						onChange={(e) =>
							updateNodeData(selectedNode.id, {
								label: e.target.value,
							})
						}
						placeholder="e.g. Compressor Outlet"
						className="p-2 border rounded border-stone-200 outline-none focus:ring-2 focus:ring-blue-100"
					/>
				</div>

				<div className="flex flex-col gap-2 opacity-60">
					<label
						htmlFor={`node_id_${selectedNode.id}`}
						className="text-xs font-bold text-gray-500 uppercase"
					>
						Node ID
					</label>
					<input
						id={`node_id_${selectedNode.id}`}
						type="text"
						value={selectedNode.id}
						className="p-1 border rounded bg-gray-50 cursor-not-allowed text-xs"
						disabled
					/>
				</div>
				{selectedNode.type === "plenum" && (
					<InitialGuessEditor node={selectedNode} />
				)}

				{selectedNode.type === "mass_boundary" && (
					<>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`m_dot_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Mass Flow (kg/s)
							</label>
							<NumericInput
								id={`m_dot_${selectedNode.id}`}
								value={selectedNode.data.m_dot || 1.0}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										m_dot: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`T_total_mb_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								T_total (K)
							</label>
							<NumericInput
								id={`T_total_mb_${selectedNode.id}`}
								value={selectedNode.data.T_total || 300}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										T_total: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<CompositionEditor
							key={selectedNode.id}
							nodeId={selectedNode.id}
							data={selectedNode.data}
						/>
						<InitialGuessEditor node={selectedNode} />
					</>
				)}

				{selectedNode.type === "pressure_boundary" && (
					<>
						<UnitInput
							label="P_total"
							value={selectedNode.data.P_total || 101325}
							onChange={(val) =>
								updateNodeData(selectedNode.id, { P_total: val })
							}
						/>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`T_total_pb_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								T_total (K)
							</label>
							<NumericInput
								id={`T_total_pb_${selectedNode.id}`}
								value={selectedNode.data.T_total || 300}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										T_total: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<CompositionEditor
							key={selectedNode.id}
							nodeId={selectedNode.id}
							data={selectedNode.data}
						/>
						<InitialGuessEditor node={selectedNode} />
					</>
				)}

				{selectedNode.type === "pipe" && (
					<>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`L_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Length (m)
							</label>
							<NumericInput
								id={`L_${selectedNode.id}`}
								value={selectedNode.data.L || 1.0}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										L: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`D_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Diameter (m)
							</label>
							<NumericInput
								id={`D_${selectedNode.id}`}
								value={selectedNode.data.D || 0.1}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										D: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`roughness_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Roughness (m)
							</label>
							<NumericInput
								id={`roughness_${selectedNode.id}`}
								value={selectedNode.data.roughness || 1e-5}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										roughness: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Regime
							</label>
							<select
								className="p-2 border rounded bg-white text-xs"
								value={selectedNode.data.regime || "default"}
								onChange={(e) =>
									updateNodeData(selectedNode.id, { regime: e.target.value })
								}
							>
								<option value="default">Default (Global)</option>
								<option value="incompressible">Forced Incompressible</option>
								<option value="compressible">Forced Compressible</option>
							</select>
						</div>
						<InitialGuessEditor node={selectedNode} />
					</>
				)}

				{selectedNode.type === "orifice" && (
					<>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`area_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Area (m²)
							</label>
							<NumericInput
								id={`area_${selectedNode.id}`}
								value={selectedNode.data.area || 0.01}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										area: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex items-center justify-between gap-2">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Auto Cd (correlation)
							</label>
							<input
								type="checkbox"
								id={`auto_Cd_${selectedNode.id}`}
								checked={selectedNode.data.auto_Cd !== false}
								onChange={(e) =>
									updateNodeData(selectedNode.id, { auto_Cd: e.target.checked })
								}
								className="w-4 h-4 accent-blue-500"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`Cd_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
								title={
									selectedNode.data.auto_Cd !== false
										? "Used as fallback when correlation is unavailable"
										: ""
								}
							>
								{selectedNode.data.auto_Cd !== false ? "Cd (fallback)" : "Cd"}
							</label>
							<NumericInput
								id={`Cd_${selectedNode.id}`}
								value={selectedNode.data.Cd || 0.6}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										Cd: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Regime
							</label>
							<select
								className="p-2 border rounded bg-white text-xs"
								value={selectedNode.data.regime || "default"}
								onChange={(e) =>
									updateNodeData(selectedNode.id, { regime: e.target.value })
								}
							>
								<option value="default">Default (Global)</option>
								<option value="incompressible">Forced Incompressible</option>
								<option value="compressible">Forced Compressible</option>
							</select>
						</div>
						<InitialGuessEditor node={selectedNode} />
					</>
				)}

				{selectedNode.type === "combustor" && (
					<div className="flex flex-col gap-4">
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`method_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Combustion Method
							</label>
							<select
								id={`method_${selectedNode.id}`}
								className="p-2 border rounded bg-white text-sm"
								value={selectedNode.data.method || "complete"}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										method: e.target.value,
									})
								}
							>
								<option value="complete">Complete (Fast)</option>
								<option value="equilibrium">Chemical Equilibrium</option>
							</select>
						</div>

						<div className="grid grid-cols-2 gap-3">
							<div className="flex flex-col gap-2">
								<label
									htmlFor={`area_comb_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									Area (m²)
								</label>
								<NumericInput
									id={`area_comb_${selectedNode.id}`}
									value={selectedNode.data.area || 0.1}
									onChange={(val) =>
										updateNodeData(selectedNode.id, {
											area: val,
										})
									}
									className="p-2 border rounded"
								/>
							</div>
							<div className="flex flex-col gap-2">
								<label
									htmlFor={`Dh_comb_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									Dh (m)
								</label>
								<NumericInput
									id={`Dh_comb_${selectedNode.id}`}
									value={selectedNode.data.Dh}
									placeholder="Auto"
									onChange={(val) =>
										updateNodeData(selectedNode.id, {
											Dh: val || undefined,
										})
									}
									className="p-2 border rounded"
								/>
							</div>
						</div>
						<p className="text-[9px] text-gray-400 italic -mt-2">
							Dh defaults to sqrt(4*Area/π) if omitted.
						</p>

						<InitialGuessEditor node={selectedNode} />
					</div>
				)}

				{selectedNode.type === "momentum_chamber" && (
					<div className="flex flex-col gap-4">
						<div className="grid grid-cols-2 gap-3">
							<div className="flex flex-col gap-2">
								<label
									htmlFor={`area_mom_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									Area (m²)
								</label>
								<NumericInput
									id={`area_mom_${selectedNode.id}`}
									value={selectedNode.data.area || 0.1}
									onChange={(val) =>
										updateNodeData(selectedNode.id, {
											area: val,
										})
									}
									className="p-2 border rounded"
								/>
							</div>
							<div className="flex flex-col gap-2">
								<label
									htmlFor={`Dh_mom_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									Dh (m)
								</label>
								<NumericInput
									id={`Dh_mom_${selectedNode.id}`}
									value={selectedNode.data.Dh}
									placeholder="Auto"
									onChange={(val) =>
										updateNodeData(selectedNode.id, {
											Dh: val || undefined,
										})
									}
									className="p-2 border rounded"
								/>
							</div>
						</div>
						<p className="text-[9px] text-gray-400 italic -mt-2">
							Dh defaults to sqrt(4*Area/π) if omitted.
						</p>

						<InitialGuessEditor node={selectedNode} />
					</div>
				)}

				{selectedNode.data.result && (
					<div className="mt-4 p-3 bg-stone-50 border border-stone-200 rounded">
						<div className="flex justify-between items-center mb-3">
							<h3 className="text-xs font-bold text-stone-500 uppercase flex items-center gap-2">
								<Activity size={14} />
								Live Telemetry
							</h3>
						</div>

						{selectedNode.data.result.state ? (
							<div className="flex flex-col gap-3">
								{selectedNode.data.result.state.phi !== undefined && (
									<div className="mb-1 pb-3 border-b border-stone-200">
										<p className="text-[10px] uppercase font-bold text-stone-400 mb-2 mt-[-2px]">
											Combustion
										</p>
										<div className="grid grid-cols-2 gap-y-2 text-xs">
											<div className="flex flex-col">
												<span className="text-stone-400 text-[10px] uppercase font-bold text-nowrap">
													Equivalence Ratio
												</span>
												<span
													className={`font-mono font-bold ${
														selectedNode.data.result.state.phi > 1.05
															? "text-orange-600"
															: selectedNode.data.result.state.phi < 0.95
															? "text-blue-600"
															: "text-green-600"
													}`}
												>
													Φ = {selectedNode.data.result.state.phi.toFixed(3)}
												</span>
											</div>
											{selectedNode.data.result.state.theta !== undefined && (
												<div className="flex flex-col">
													<span className="text-stone-400 text-[10px] uppercase font-bold text-nowrap">
														Temp Rise Ratio
													</span>
													<span className="font-mono font-bold">
														θ = {selectedNode.data.result.state.theta.toFixed(3)}
													</span>
												</div>
											)}
										</div>
									</div>
								)}
								<div className="grid grid-cols-2 gap-y-2 text-xs">
									<div className="flex flex-col">
										<span className="text-stone-400 text-[10px] uppercase font-bold text-nowrap">
											{selectedNode.data.result.state.P_total !== undefined &&
											Math.abs(
												selectedNode.data.result.state.P -
													selectedNode.data.result.state.P_total,
											) > 1e-1
												? "St/Tot Pressure"
												: "Pressure"}
										</span>
										<span className="font-mono font-bold flex gap-1 items-baseline">
											{unitPreferences.pressure === "Pa" &&
												`${selectedNode.data.result.state.P.toFixed(0)} Pa`}
											{unitPreferences.pressure === "kPa" &&
												`${(selectedNode.data.result.state.P / 1e3).toFixed(2)} kPa`}
											{unitPreferences.pressure === "MPa" &&
												`${(selectedNode.data.result.state.P / 1e6).toFixed(4)} MPa`}

											{selectedNode.data.result.state.P_total !== undefined &&
												Math.abs(
													selectedNode.data.result.state.P -
														selectedNode.data.result.state.P_total,
												) > 1e-1 && (
													<span className="text-[10px] text-stone-500 font-normal">
														/{" "}
														{unitPreferences.pressure === "Pa" &&
															`${selectedNode.data.result.state.P_total.toFixed(0)}`}
														{unitPreferences.pressure === "kPa" &&
															`${(selectedNode.data.result.state.P_total / 1e3).toFixed(2)}`}
														{unitPreferences.pressure === "MPa" &&
															`${(selectedNode.data.result.state.P_total / 1e6).toFixed(4)}`}
													</span>
												)}
										</span>
									</div>
									<div className="flex flex-col">
										<span className="text-stone-400 text-[10px] uppercase font-bold">
											Temperature
										</span>
										<span className="font-mono font-bold">
											{selectedNode.data.result.state.T.toFixed(2)} K
										</span>
									</div>
									<div className="flex flex-col">
										<span className="text-stone-400 text-[10px] uppercase font-bold">
											Density
										</span>
										<span className="font-mono">
											{selectedNode.data.result.state.rho?.toFixed(4)} kg/m³
										</span>
									</div>
									<div className="flex flex-col">
										{selectedNode.type !== "plenum" && (
											<div className="flex flex-col gap-1">
												<span className="text-stone-400 text-[10px] uppercase font-bold">
													Mach
												</span>
												<span className="font-mono">
													{selectedNode.data.result.state.mach?.toFixed(4) ||
														"0.0000"}
												</span>
											</div>
										)}
									</div>
								</div>

								<div className="border-t border-stone-200 pt-3">
									<div className="text-[10px] font-bold text-stone-400 uppercase mb-2">
										Thermodynamics
									</div>
									<div className="grid grid-cols-2 gap-y-1 text-xs font-mono">
										<span className="text-stone-500">h:</span>
										<span className="text-right">
											{(selectedNode.data.result.state.h / 1e3).toFixed(2)}{" "}
											kJ/kg
										</span>
										<span className="text-stone-500">s:</span>
										<span className="text-right">
											{(selectedNode.data.result.state.s / 1e3).toFixed(4)}{" "}
											kJ/kg-K
										</span>
										{selectedNode.data.result.state.cp !== undefined && (
											<>
												<span className="text-stone-500">Cp:</span>
												<span className="text-right">
													{selectedNode.data.result.state.cp.toFixed(2)} J/kg-K
												</span>
												<span className="text-stone-500">gamma:</span>
												<span className="text-right">
													{selectedNode.data.result.state.gamma.toFixed(4)}
												</span>
												<span className="text-stone-500">a:</span>
												<span className="text-right">
													{selectedNode.data.result.state.a.toFixed(2)} m/s
												</span>
											</>
										)}
									</div>
								</div>

								{selectedNode.data.result.state.mu !== undefined && (
									<div className="border-t border-stone-200 pt-3">
										<div className="text-[10px] font-bold text-stone-400 uppercase mb-2">
											Transport
										</div>
										<div className="grid grid-cols-2 gap-y-1 text-xs font-mono">
											<span className="text-stone-500">μ (viscosity):</span>
											<span className="text-right">
												{selectedNode.data.result.state.mu.toExponential(2)}{" "}
												Pa·s
											</span>
											<span className="text-stone-500">k (conductivity):</span>
											<span className="text-right">
												{selectedNode.data.result.state.k.toExponential(2)}{" "}
												W/m-K
											</span>
											<span className="text-stone-500">Pr:</span>
											<span className="text-right">
												{selectedNode.data.result.state.Pr.toFixed(3)}
											</span>
										</div>
									</div>
								)}

								{selectedNode.data.result.state.X && (
									<div className="border-t border-stone-200 pt-3">
										<div className="text-[10px] font-bold text-stone-400 uppercase mb-2">
											Composition (Mole %)
										</div>
										<div className="max-h-40 overflow-y-auto pr-1 flex flex-col gap-1">
											{selectedNode.data.result.state.X.map(
												(x: number, i: number) => ({
													x,
													name: speciesMetadata?.names[i] || `Species ${i}`,
												}),
											)
												.filter((item: any) => item.x > 1e-6)
												.map((item: any, i: number) => (
													<div
														key={i}
														className="flex justify-between font-mono text-[10px] bg-white p-1 rounded border border-stone-100"
													>
														<span className="text-stone-500">{item.name}</span>
														<span className="font-bold">
															{(item.x * 100).toFixed(4)}%
														</span>
													</div>
												))}
										</div>
									</div>
								)}
							</div>
						) : (
							<div className="flex flex-col gap-3">
								<div className="flex flex-col gap-1">
									<span className="text-stone-400 text-[10px] uppercase font-bold">
										Mass Flow
									</span>
									<span className="font-mono text-lg font-bold">
										{selectedNode.data.result.m_dot?.toFixed(4)} kg/s
									</span>
								</div>

								{/* Element Diagnostics (Mach, P-ratio) */}
								<div className="border-t border-stone-200 pt-3 flex flex-col gap-2">
									<div className="text-[10px] font-bold text-stone-400 uppercase">
										Diagnostics
									</div>
									<div className="grid grid-cols-2 gap-2">
										{selectedNode.data.result.mach !== undefined && (
											<div className="flex flex-col">
												<span className="text-stone-400 text-[8px] uppercase">
													Mach (Throat)
												</span>
												<span className="font-mono text-xs font-bold">
													{selectedNode.data.result.mach.toFixed(4)}
												</span>
											</div>
										)}
										{selectedNode.data.result.mach_in !== undefined && (
											<div className="flex flex-col">
												<span className="text-stone-400 text-[8px] uppercase">
													Mach In
												</span>
												<span className="font-mono text-xs font-bold">
													{selectedNode.data.result.mach_in.toFixed(4)}
												</span>
											</div>
										)}
										{selectedNode.data.result.mach_out !== undefined && (
											<div className="flex flex-col">
												<span className="text-stone-400 text-[8px] uppercase">
													Mach Out
												</span>
												<span className="font-mono text-xs font-bold">
													{selectedNode.data.result.mach_out.toFixed(4)}
												</span>
											</div>
										)}
										{(selectedNode.data.result.p_ratio !== undefined ||
											selectedNode.data.result.p_ratio_total !== undefined) && (
											<div className="flex flex-col">
												<span className="text-stone-400 text-[8px] uppercase">
													P Ratio
												</span>
												<span className="font-mono text-xs font-bold">
													{(
														selectedNode.data.result.p_ratio ||
														selectedNode.data.result.p_ratio_total
													).toFixed(4)}
												</span>
											</div>
										)}
										{selectedNode.data.result.Cd !== undefined && (
											<div className="flex flex-col col-span-2">
												<span className="text-stone-400 text-[8px] uppercase">
													Cd (used)
												</span>
												<span className="font-mono text-xs font-bold">
													{selectedNode.data.result.Cd.toFixed(4)}
													{selectedNode.data.result.auto_Cd === 1 && (
														<span className="ml-1 text-[8px] text-blue-500 font-normal">
															(corr.)
														</span>
													)}
												</span>
											</div>
										)}
									</div>
								</div>

								{selectedNode.data.result.mu !== undefined && (
									<div className="border-t border-stone-200 pt-3 flex flex-col gap-2">
										<div className="text-[10px] font-bold text-stone-400 uppercase">
											Transport
										</div>
										<div className="grid grid-cols-2 gap-x-2 gap-y-1 text-xs font-mono">
											{selectedNode.data.result.Re !== undefined && (
												<>
													<div className="text-stone-500">Re:</div>
													<div className="text-right font-bold">
														{selectedNode.data.result.Re.toLocaleString("en", {
															maximumFractionDigits: 0,
														})}
													</div>
												</>
											)}
											<div className="text-stone-500">μ (viscosity):</div>
											<div className="text-right font-bold">
												{selectedNode.data.result.mu.toExponential(2)} Pa·s
											</div>
											<div className="text-stone-500">k:</div>
											<div className="text-right font-bold">
												{selectedNode.data.result.k.toExponential(2)} W/m-K
											</div>
											<div className="text-stone-500">Pr (Prandtl):</div>
											<div className="text-right font-bold">
												{selectedNode.data.result.Pr.toFixed(3)}
											</div>
										</div>
									</div>
								)}
							</div>
						)}
					</div>
				)}
			</aside>
		);
	}

	// Edge Inspector
	if (selectedEdge && selectedEdge.data?.type === "thermal") {
		return (
			<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto">
				<h2 className="text-lg font-bold border-b pb-2 uppercase text-orange-600">
					Thermal Wall
				</h2>
				<div className="flex flex-col gap-4">
					<div className="flex flex-col gap-2">
						<label className="text-xs font-bold text-gray-500 uppercase">
							Thickness (m)
						</label>
						<NumericInput
							value={selectedEdge.data.thickness || 0.003}
							onChange={(val) =>
								updateEdgeData(selectedEdge.id, { thickness: val })
							}
							className="p-2 border rounded"
						/>
					</div>
					<div className="flex flex-col gap-2">
						<label className="text-xs font-bold text-gray-500 uppercase">
							Conductivity (W/mK)
						</label>
						<NumericInput
							value={selectedEdge.data.conductivity || 20.0}
							onChange={(val) =>
								updateEdgeData(selectedEdge.id, { conductivity: val })
							}
							className="p-2 border rounded"
						/>
					</div>
					<div className="flex flex-col gap-2">
						<label className="text-xs font-bold text-gray-500 uppercase">
							Area (m²)
						</label>
						<NumericInput
							value={selectedEdge.data.area || 0.05}
							onChange={(val) => updateEdgeData(selectedEdge.id, { area: val })}
							className="p-2 border rounded"
						/>
					</div>
				</div>
			</aside>
		);
	}

	return (
		<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 italic text-gray-400">
			Connection is structural. Topographic link only.
		</aside>
	);
};

const InitialGuessEditor = ({ node }: { node: any }) => {
	const { updateNodeData } = useStore();
	const initialGuess = node.data.initial_guess || {};

	const updateGuess = (key: string, val: number | undefined) => {
		const newGuess = { ...initialGuess };
		if (val === undefined || Number.isNaN(val)) {
			delete newGuess[key];
		} else {
			newGuess[key] = val;
		}
		updateNodeData(node.id, { initial_guess: newGuess });
	};

	const guessKeys = ["P", "T", "m_dot"];

	return (
		<div className="mt-4 border-t pt-4">
			<h3 className="text-xs font-bold text-stone-500 uppercase mb-2">
				Initial Guess Overrides
			</h3>
			<div className="flex flex-col gap-2">
				{guessKeys.map((key) => (
					<div key={key} className="flex items-center justify-between gap-2">
						<span className="text-[10px] font-mono text-gray-500">{key}</span>
						<NumericInput
							placeholder="Auto"
							value={initialGuess[key]}
							onChange={(val) => updateGuess(key, val || undefined)}
							className="p-1 border rounded text-xs w-24"
						/>
					</div>
				))}
			</div>
		</div>
	);
};

export default Inspector;
