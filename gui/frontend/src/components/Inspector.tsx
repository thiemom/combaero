import { Activity, Crosshair, RotateCw } from "lucide-react";
import React from "react";
import useStore from "../store/useStore";
import {
	BASIC_STATE_KEYS,
	QUANTITY_CATALOGUE,
	TRANSPORT_KEYS,
} from "../utils/quantities";
import { validateNetwork } from "../utils/validation";
import AreaInput from "./AreaInput";
import CompositionEditor from "./CompositionEditor";
import LengthInput from "./LengthInput";
import NumericInput from "./NumericInput";
import SurfaceEnhancementInspector from "./SurfaceEnhancementInspector";
import UnitInput from "./UnitInput";
import WallLayersEditor from "./WallLayersEditor";

const Inspector = () => {
	const {
		nodes,
		edges,
		updateNodeData,
		updateEdgeData,
		speciesMetadata,
		unitPreferences,
		isExporting,
		exportNetworkResults,
	} = useStore();

	const validationErrors = validateNetwork(nodes);

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
		const targetElement =
			nodes.find((n) => n.id === probeData.target_id) ||
			edges.find((e) => e.id === probeData.target_id);
		const result = targetElement?.data?.result;

		// Build available quantity list
		const availableKeys: { key: string; label: string }[] = [];
		const seenKeys = new Set<string>();

		const addKey = (k: string, isState = false) => {
			const fullKey = isState ? `state.${k}` : k;
			if (seenKeys.has(fullKey)) return;
			seenKeys.add(fullKey);

			const meta = QUANTITY_CATALOGUE[k];
			availableKeys.push({
				key: fullKey,
				label: meta ? `${meta.label} (${k})` : k,
			});
		};

		if (result) {
			// Populate from actual results
			if (result.state) {
				for (const k of Object.keys(result.state)) {
					if (typeof result.state[k] === "number") addKey(k, true);
				}
			}
			for (const k of Object.keys(result)) {
				if (k !== "state" && typeof result[k] === "number") addKey(k);
			}
		} else {
			// Fallback: show common state and transport quantities if no result yet
			for (const k of [...BASIC_STATE_KEYS, ...TRANSPORT_KEYS]) {
				addKey(k);
			}
			// Specifically add area change ones if target is area change
			if (targetElement?.type === "area_change") {
				addKey("zeta");
				addKey("mach_small");
				addKey("mach_in");
				addKey("mach_out");
				addKey("ratio");
			}
		}

		// Sort by label
		availableKeys.sort((a, b) => a.label.localeCompare(b.label));

		const slotSelect = (
			slot: "slot1_key" | "slot2_key",
			current: string | undefined,
		) => (
			<select
				className="p-1.5 border rounded text-xs bg-white w-full"
				value={current ?? ""}
				onChange={(e) =>
					selectedNode &&
					updateNodeData(selectedNode.id, {
						[slot]: e.target.value || undefined,
					})
				}
			>
				<option value="">— select quantity —</option>
				{availableKeys.map((item) => (
					<option key={item.key} value={item.key}>
						{item.label}
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
						className="p-1.5 border rounded text-xs w-full focus:ring-1 focus:ring-blue-500 outline-none"
					/>
				</div>

				{/* Target */}
				<div className="flex flex-col gap-1">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Target Element
					</label>
					<select
						className="p-1.5 border rounded text-xs bg-white w-full focus:ring-1 focus:ring-blue-500 outline-none"
						value={probeData.target_id ?? ""}
						onChange={(e) =>
							updateNodeData(selectedNode.id, {
								target_id: e.target.value || undefined,
							})
						}
					>
						<optgroup label="Network Nodes">
							<option value="">— select target —</option>
							{nodes
								.filter((n) => n.id !== selectedNode?.id && n.type !== "probe")
								.map((n) => (
									<option key={n.id} value={n.id}>
										{n.data?.label || n.type?.replace(/_/g, " ") || n.id}
									</option>
								))}
						</optgroup>
						<optgroup label="Thermal Walls">
							{edges
								.filter((e) => e.data?.type === "thermal")
								.map((e) => (
									<option key={e.id} value={e.id}>
										{e.data?.label || `Thermal Wall (${e.id})`}
									</option>
								))}
						</optgroup>
					</select>
					{probeData.target_id && !targetElement && (
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
							onClick={exportNetworkResults}
							disabled={isExporting}
							className="text-[10px] bg-stone-100 hover:bg-stone-200 px-2 py-1 rounded border border-stone-300 transition-colors uppercase font-bold disabled:opacity-50 min-w-[80px]"
						>
							{isExporting ? "Exporting..." : "Export CSV"}
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
						placeholder="e.g. Custom Label"
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
								value={selectedNode.data.m_dot ?? 1.0}
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
								value={selectedNode.data.T_total ?? 300}
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
							value={selectedNode.data.P_total ?? 101325}
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
								value={selectedNode.data.T_total ?? 300}
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

				{selectedNode.type === "channel" && (
					<>
						<LengthInput
							id={`L_${selectedNode.id}`}
							label="Length Flow Path"
							value={selectedNode.data.L ?? 1.0}
							onChange={(val) => updateNodeData(selectedNode.id, { L: val })}
						/>
						<LengthInput
							id={`D_${selectedNode.id}`}
							label="Channel Inner Diameter"
							value={selectedNode.data.D ?? 0.1}
							onChange={(val) => updateNodeData(selectedNode.id, { D: val })}
						/>
						<LengthInput
							id={`roughness_${selectedNode.id}`}
							label="Surface Roughness"
							value={selectedNode.data.roughness ?? 1e-5}
							onChange={(val) =>
								updateNodeData(selectedNode.id, { roughness: val })
							}
						/>
						<div className="grid grid-cols-2 gap-3">
							<div className="flex flex-col gap-2">
								<label className="text-xs font-bold text-gray-500 uppercase">
									Nu Multiplier
								</label>
								<NumericInput
									value={selectedNode.data.Nu_multiplier ?? 1.0}
									onChange={(val) =>
										updateNodeData(selectedNode.id, { Nu_multiplier: val })
									}
									className="p-2 border rounded"
								/>
							</div>
							<div className="flex flex-col gap-2">
								<label className="text-xs font-bold text-gray-500 uppercase">
									F Multiplier
								</label>
								<NumericInput
									value={selectedNode.data.f_multiplier ?? 1.0}
									onChange={(val) =>
										updateNodeData(selectedNode.id, { f_multiplier: val })
									}
									className="p-2 border rounded"
								/>
							</div>
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
						<SurfaceEnhancementInspector
							surface={selectedNode.data.surface || { type: "smooth" }}
							onChange={(surface) =>
								updateNodeData(selectedNode.id, { surface })
							}
						/>
						<InitialGuessEditor node={selectedNode} />
					</>
				)}

				{selectedNode.type === "orifice" && (
					<>
						<LengthInput
							id={`diameter_${selectedNode.id}`}
							label="Bore Diameter"
							value={selectedNode.data.diameter ?? 0.08}
							onChange={(val) =>
								updateNodeData(selectedNode.id, { diameter: val })
							}
						/>

						<div className="flex flex-col gap-2">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Discharge Model (Cd)
							</label>
							<select
								className="p-2 border rounded bg-white text-xs border-stone-200"
								value={selectedNode.data.correlation || "ReaderHarrisGallagher"}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										correlation: e.target.value,
									})
								}
							>
								<option value="ReaderHarrisGallagher">
									Reader-Harris/Gallagher (Sharp)
								</option>
								<option value="Stolz">Stolz (Corner Taps)</option>
								<option value="Miller">Miller (Simplified)</option>
								<option value="ThickPlate">Thick Plate (Sharp Edge)</option>
								<option value="RoundedEntry">Rounded Entry</option>
								<option value="fixed">Manual / Fixed Value</option>
							</select>
						</div>

						{/* Conditional Inputs based on correlation */}
						{selectedNode.data.correlation === "ThickPlate" && (
							<LengthInput
								id={`plate_thickness_${selectedNode.id}`}
								label="Plate Thickness (t)"
								value={selectedNode.data.plate_thickness ?? 0.0}
								onChange={(val) =>
									updateNodeData(selectedNode.id, { plate_thickness: val })
								}
							/>
						)}

						{selectedNode.data.correlation === "RoundedEntry" && (
							<LengthInput
								id={`edge_radius_${selectedNode.id}`}
								label="Inlet Edge Radius (r)"
								value={selectedNode.data.edge_radius ?? 0.0}
								onChange={(val) =>
									updateNodeData(selectedNode.id, { edge_radius: val })
								}
							/>
						)}

						{/* Calculated Result (for correlation models) */}
						{selectedNode.data.correlation !== "fixed" && (
							<div className="flex flex-col gap-1 mb-2 bg-blue-50/50 p-2 rounded border border-blue-100/50">
								<label className="text-[10px] font-bold text-blue-500 uppercase tracking-wider">
									Calculated Cd
								</label>
								<div className="text-sm font-mono font-bold text-blue-700">
									{selectedNode.data.result?.Cd?.toFixed(4) ||
										"— (Pending Solve)"}
								</div>
							</div>
						)}

						{/* Manual Entry (only for Fixed model) */}
						{selectedNode.data.correlation === "fixed" && (
							<div className="flex flex-col gap-1 mt-1">
								<label
									htmlFor={`Cd_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									Fixed Cd Value
								</label>
								<NumericInput
									id={`Cd_${selectedNode.id}`}
									value={selectedNode.data.Cd ?? 0.6}
									onChange={(val) =>
										updateNodeData(selectedNode.id, {
											Cd: val,
										})
									}
									className="p-1.5 h-8 text-sm border rounded bg-white"
									placeholder="0.6"
								/>
							</div>
						)}

						<div className="flex flex-col gap-2 mt-2">
							<label className="text-xs font-bold text-gray-500 uppercase tracking-wide">
								Flow Regime
							</label>
							<select
								className="p-2 border rounded bg-white text-xs border-stone-200"
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

				{selectedNode.type === "area_change" && (
					<>
						<div className="flex flex-col gap-2">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Model Type
							</label>
							<select
								className="p-2 border rounded bg-white text-xs border-stone-200"
								value={selectedNode.data.model_type || "sharp"}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										model_type: e.target.value,
									})
								}
							>
								<option value="sharp">Sharp-Edged (Sudden)</option>
								<option value="conical">Conical (Gradual)</option>
							</select>
						</div>

						{(selectedNode.data.F0 <= 0 ||
							selectedNode.data.F1 <= 0 ||
							(selectedNode.data.D_h !== undefined &&
								selectedNode.data.D_h !== null &&
								selectedNode.data.D_h < 0)) && (
							<div className="bg-amber-50 border border-amber-200 p-2 rounded text-[10px] text-amber-800 flex items-center gap-2 mb-2">
								<Activity size={14} className="shrink-0" />
								<div>
									<p className="font-bold uppercase leading-none mb-1">
										Invalid Geometry
									</p>
									<p>
										Area (F) must be &gt; 0 and Hydraulic Diameter (D_h) must be
										&ge; 0. Solve will fail.
									</p>
								</div>
							</div>
						)}

						<AreaInput
							id={`F0_${selectedNode.id}`}
							label="Upstream Area (F0)"
							value={selectedNode.data.F0 ?? 0.01}
							onChange={(val) => updateNodeData(selectedNode.id, { F0: val })}
						/>

						<AreaInput
							id={`F1_${selectedNode.id}`}
							label="Downstream Area (F1)"
							value={selectedNode.data.F1 ?? 0.02}
							onChange={(val) => updateNodeData(selectedNode.id, { F1: val })}
						/>

						{selectedNode.data.model_type === "conical" && (
							<LengthInput
								id={`length_${selectedNode.id}`}
								label="Axial Length"
								value={selectedNode.data.length ?? 0.1}
								onChange={(val) =>
									updateNodeData(selectedNode.id, { length: val })
								}
							/>
						)}

						{selectedNode.data.model_type !== "conical" && (
							<div className="flex flex-col gap-1 mt-1">
								<LengthInput
									id={`D_h_${selectedNode.id}`}
									label="Hydraulic Diameter (D_h)"
									value={selectedNode.data.D_h}
									placeholder="Circular (Auto)"
									onChange={(val) =>
										updateNodeData(selectedNode.id, { D_h: val })
									}
									onClear={() =>
										updateNodeData(selectedNode.id, { D_h: undefined })
									}
								/>
								<p className="text-[9px] text-gray-400 italic">
									Set explicitly for non-circular ducts to compute correct
									Reynolds number.
								</p>
							</div>
						)}

						<InitialGuessEditor node={selectedNode} />
					</>
				)}

				{selectedNode.type === "discrete_loss" &&
					(() => {
						// Compute default area from upstream node
						const upstreamEdge = edges.find(
							(e) => e.target === selectedNode.id && !e.data?.type,
						);
						const upstreamNode = upstreamEdge
							? nodes.find((n) => n.id === upstreamEdge.source)
							: undefined;
						const inheritedArea: number | undefined = upstreamNode?.data?.area;

						const userSetArea =
							selectedNode.data.area !== null &&
							selectedNode.data.area !== undefined;
						const displayArea: number = userSetArea
							? selectedNode.data.area
							: (inheritedArea ?? 0.1);

						// Compute available has_theta nodes (combustors)
						const thetaNodes = nodes.filter(
							(n) => n.type === "combustor" && n.id !== selectedNode.id,
						);

						// Determine auto-default: first downstream combustor, then upstream
						const downstreamEdge = edges.find(
							(e) => e.source === selectedNode.id && !e.data?.type,
						);
						const downstreamNode = downstreamEdge
							? nodes.find((n) => n.id === downstreamEdge.target)
							: undefined;
						const autoDefault =
							downstreamNode?.type === "combustor"
								? downstreamNode.id
								: upstreamNode?.type === "combustor"
									? upstreamNode.id
									: null;

						const thetaSource: string | null =
							selectedNode.data.theta_source ?? null;
						const hasThetaSource =
							thetaSource !== null && thetaSource !== "none";
						const corrType: string =
							selectedNode.data.correlation_type ?? "constant_fraction";
						const isThetaCorr =
							corrType === "linear_theta_fraction" ||
							corrType === "linear_theta_head";

						const areaWarning =
							(userSetArea ? selectedNode.data.area : inheritedArea) != null &&
							(userSetArea ? selectedNode.data.area : inheritedArea) <= 0;

						return (
							<div className="flex flex-col gap-4">
								{/* Correlation type */}
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Correlation Type
									</label>
									<select
										className="p-2 border rounded bg-white text-sm"
										value={corrType}
										onChange={(e) => {
											const t = e.target.value;
											const needsTheta =
												t === "linear_theta_fraction" ||
												t === "linear_theta_head";
											if (
												needsTheta &&
												!hasThetaSource &&
												thetaNodes.length > 0
											) {
												updateNodeData(selectedNode.id, {
													correlation_type: t,
													theta_source: autoDefault ?? thetaNodes[0].id,
												});
											} else if (needsTheta && thetaNodes.length === 0) {
												return;
											} else {
												updateNodeData(selectedNode.id, {
													correlation_type: t,
												});
											}
										}}
									>
										<option value="constant_fraction">
											Constant fraction — dP/P = ξ
										</option>
										<option value="constant_head">
											Constant head — dP = ζ · q
										</option>
										{thetaNodes.length > 0 || hasThetaSource ? (
											<option value="linear_theta_fraction">
												Linear Θ fraction — dP/P = k·Θ + ξ₀
											</option>
										) : null}
										{thetaNodes.length > 0 || hasThetaSource ? (
											<option value="linear_theta_head">
												Linear Θ head — dP = (k·Θ + ζ₀) · q
											</option>
										) : null}
									</select>
									{isThetaCorr && !hasThetaSource && (
										<p className="text-[9px] text-amber-500 italic">
											Select a theta source below to enable this correlation.
										</p>
									)}
								</div>

								{/* Flow Area — always shown, inherited from upstream by default */}
								<div className="flex flex-col gap-2">
									<div className="flex items-center justify-between">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Flow Area (m²)
										</label>
										{userSetArea && (
											<button
												type="button"
												className="text-[9px] text-blue-500 hover:underline"
												onClick={() =>
													updateNodeData(selectedNode.id, { area: null })
												}
											>
												Reset to inherited
											</button>
										)}
									</div>
									<NumericInput
										value={displayArea}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { area: val })
										}
										onClear={() =>
											// Clearing the field is semantically the same as clicking
											// "Reset to inherited": we drop the user override so the
											// backend auto-infers from the upstream node again.
											updateNodeData(selectedNode.id, { area: null })
										}
										min={0}
										className={`p-2 border rounded ${
											userSetArea
												? "text-gray-900 font-semibold"
												: "text-gray-400"
										}`}
										placeholder={
											inheritedArea != null
												? `${inheritedArea.toFixed(4)} (from upstream)`
												: "0.1"
										}
									/>
									{!userSetArea && inheritedArea != null && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Inherited from upstream ({upstreamNode?.id}). Edit to
											override.
										</p>
									)}
									{!userSetArea && inheritedArea == null && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Connect an upstream node to auto-discover area.
										</p>
									)}
									{areaWarning && (
										<div className="text-[10px] text-amber-600 bg-amber-50 border border-amber-200 rounded px-2 py-1">
											⚠ Area ≤ 0 — head loss will be undefined.
										</div>
									)}
								</div>

								{/* Parameters */}
								{(corrType === "constant_fraction" ||
									corrType === "linear_theta_fraction") && (
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											{corrType === "constant_fraction"
												? "ξ — loss fraction [−]"
												: "ξ₀ — cold-flow base fraction [−]"}
										</label>
										<NumericInput
											value={
												corrType === "constant_fraction"
													? (selectedNode.data.xi ?? 0.03)
													: (selectedNode.data.xi0 ?? 0.02)
											}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													[corrType === "constant_fraction" ? "xi" : "xi0"]:
														val,
												})
											}
											min={0}
											className="p-2 border rounded"
										/>
										{corrType === "linear_theta_fraction" && (
											<>
												<label className="text-xs font-bold text-gray-500 uppercase">
													k — Theta sensitivity [−/−]
												</label>
												<NumericInput
													value={selectedNode.data.k ?? 0.001}
													onChange={(val) =>
														updateNodeData(selectedNode.id, { k: val })
													}
													min={0}
													className="p-2 border rounded"
												/>
												<p className="text-[9px] text-gray-400 italic">
													P_out = P_in · (1 − k·Θ − ξ₀). Θ = T_burned/T_in − 1
												</p>
											</>
										)}
										{corrType === "constant_fraction" && (
											<p className="text-[9px] text-gray-400 italic">
												P_out = P_in · (1 − ξ). e.g. 0.03 = 3%
											</p>
										)}
									</div>
								)}

								{(corrType === "constant_head" ||
									corrType === "linear_theta_head") && (
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											{corrType === "constant_head"
												? "ζ — Euler loss coeff. [−]"
												: "ζ₀ — cold-flow Euler coeff. [−]"}
										</label>
										<NumericInput
											value={
												corrType === "constant_head"
													? (selectedNode.data.zeta ?? 1.0)
													: (selectedNode.data.zeta0 ?? 1.0)
											}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													[corrType === "constant_head" ? "zeta" : "zeta0"]:
														val,
												})
											}
											min={0}
											className="p-2 border rounded"
										/>
										{corrType === "linear_theta_head" && (
											<>
												<label className="text-xs font-bold text-gray-500 uppercase">
													k — Theta sensitivity [−/−]
												</label>
												<NumericInput
													value={selectedNode.data.k ?? 0.001}
													onChange={(val) =>
														updateNodeData(selectedNode.id, { k: val })
													}
													min={0}
													className="p-2 border rounded"
												/>
												<p className="text-[9px] text-gray-400 italic">
													dP = (k·Θ + ζ₀) · ½ρv². Θ = T_burned/T_in − 1
												</p>
											</>
										)}
									</div>
								)}

								{/* Theta Source */}
								<div className="flex flex-col gap-2 border-t pt-3">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Theta Source
									</label>
									<select
										className="p-2 border rounded bg-white text-sm"
										value={thetaSource ?? ""}
										onChange={(e) => {
											const val = e.target.value || null;
											updateNodeData(selectedNode.id, { theta_source: val });
										}}
									>
										<option value="">
											Auto{autoDefault ? ` (→ ${autoDefault})` : " (cold-flow)"}
										</option>
										<option value="none">None (cold-flow)</option>
										{thetaNodes.map((n) => (
											<option key={n.id} value={n.id}>
												{n.data?.label
													? `${n.data.label} (${n.id})`
													: `Combustor: ${n.id}`}
											</option>
										))}
									</select>
									{thetaNodes.length === 0 && (
										<p className="text-[9px] text-gray-400 italic">
											No combustor nodes in network. Only ξ/ζ correlations
											available.
										</p>
									)}
								</div>

								<div className="border-t pt-3">
									<SurfaceEnhancementInspector
										surface={selectedNode.data.surface || { type: "smooth" }}
										onChange={(surface) =>
											updateNodeData(selectedNode.id, { surface })
										}
									/>
								</div>

								<InitialGuessEditor node={selectedNode} />
							</div>
						);
					})()}

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

						<div className="grid grid-cols-2 gap-2 pb-2">
							<AreaInput
								id={`area_comb_${selectedNode.id}`}
								label="Area"
								value={selectedNode.data.area ?? 0.1}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										area: val,
									})
								}
							/>
							<div className="flex flex-col gap-1">
								<LengthInput
									id={`Dh_comb_${selectedNode.id}`}
									label="Dh"
									value={selectedNode.data.Dh}
									placeholder={Math.sqrt(
										(4 * (selectedNode.data.area ?? 0.1)) / Math.PI,
									)}
									onChange={(val) =>
										updateNodeData(selectedNode.id, {
											Dh: val || undefined,
										})
									}
								/>
							</div>
						</div>
						<p className="text-[9px] text-gray-400 italic -mt-2 mb-2">
							Dh defaults to sqrt(4*Area/π) if omitted.
						</p>
						<div className="flex flex-col gap-2">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Nu Multiplier
							</label>
							<NumericInput
								id={`Nu_multiplier_comb_${selectedNode.id}`}
								value={selectedNode.data.Nu_multiplier ?? 1.0}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										Nu_multiplier: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<SurfaceEnhancementInspector
							surface={selectedNode.data.surface || { type: "smooth" }}
							onChange={(surface) =>
								updateNodeData(selectedNode.id, { surface })
							}
						/>

						<InitialGuessEditor node={selectedNode} />
					</div>
				)}

				{selectedNode.type === "momentum_chamber" && (
					<div className="flex flex-col gap-4">
						<div className="grid grid-cols-2 gap-2 pb-2">
							<AreaInput
								id={`area_mom_${selectedNode.id}`}
								label="Area"
								value={selectedNode.data.area ?? 0.1}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										area: val,
									})
								}
							/>
							<div className="flex flex-col gap-1">
								<LengthInput
									id={`Dh_mom_${selectedNode.id}`}
									label="Dh"
									value={selectedNode.data.Dh}
									placeholder={Math.sqrt(
										(4 * (selectedNode.data.area ?? 0.1)) / Math.PI,
									)}
									onChange={(val) =>
										updateNodeData(selectedNode.id, {
											Dh: val || undefined,
										})
									}
								/>
							</div>
						</div>
						<p className="text-[9px] text-gray-400 italic -mt-2 mb-2">
							Dh defaults to sqrt(4*Area/π) if omitted.
						</p>
						<div className="flex flex-col gap-2">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Nu Multiplier
							</label>
							<NumericInput
								id={`Nu_multiplier_mom_${selectedNode.id}`}
								value={selectedNode.data.Nu_multiplier ?? 1.0}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										Nu_multiplier: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<SurfaceEnhancementInspector
							surface={selectedNode.data.surface || { type: "smooth" }}
							onChange={(surface) =>
								updateNodeData(selectedNode.id, { surface })
							}
						/>

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
								{/* 1. Basic State Section (Always visible) */}
								<div className="grid grid-cols-2 gap-x-2 gap-y-3">
									{BASIC_STATE_KEYS.map((key) => {
										const val = selectedNode.data.result.state[key];
										if (val === undefined || val === null) return null;
										const meta = QUANTITY_CATALOGUE[key];
										return (
											<div key={key} className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold">
													{meta?.label || key}
												</span>
												<span className="font-mono text-xs font-bold whitespace-nowrap">
													{meta?.format(val) || val.toFixed(4)}{" "}
													<span className="text-[9px] font-normal text-stone-400 ml-0.5">
														{meta?.unit}
													</span>
												</span>
											</div>
										);
									})}
								</div>

								{/* 2. Advanced/Transport Section (Conditional) */}
								{TRANSPORT_KEYS.some(
									(k) => selectedNode.data.result.state[k] !== undefined,
								) && (
									<div className="border-t border-stone-100 pt-3">
										<div className="text-[10px] font-bold text-stone-400 uppercase mb-2">
											Advanced Properties
										</div>
										<div className="grid grid-cols-2 gap-x-2 gap-y-2">
											{TRANSPORT_KEYS.map((key) => {
												const val = selectedNode.data.result.state[key];
												if (val === undefined || val === null) return null;
												const meta = QUANTITY_CATALOGUE[key];
												return (
													<div
														key={key}
														className="flex justify-between items-baseline border-b border-dotted border-stone-100 pb-1"
													>
														<span className="text-stone-500 text-[10px]">
															{meta?.label || key}
														</span>
														<span className="font-mono text-[10px] font-bold">
															{meta?.format(val) || val.toFixed(4)}
														</span>
													</div>
												);
											})}
										</div>
									</div>
								)}

								{/* 3. Composition Section */}
								{selectedNode.data.result.state.X && (
									<div className="border-t border-stone-100 pt-3">
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
														className="flex justify-between font-mono text-[9px] bg-white p-1 rounded border border-stone-50"
													>
														<span className="text-stone-500 text-ellipsis overflow-hidden shrink pr-2">
															{item.name}
														</span>
														<span className="font-bold flex-shrink-0">
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
								<div className="text-[10px] font-bold text-stone-400 uppercase mb-1">
									Element Diagnostics
								</div>
								<div className="grid grid-cols-2 gap-x-2 gap-y-3">
									{Object.keys(QUANTITY_CATALOGUE)
										.filter((k) => !k.startsWith("state."))
										.map((key) => {
											const val = selectedNode.data.result[key];
											if (val === undefined || val === null) return null;
											const meta = QUANTITY_CATALOGUE[key];

											// Special styling for primary metrics like m_dot
											if (key === "m_dot") {
												return (
													<div
														key={key}
														className="col-span-2 flex flex-col mb-1 border-b border-stone-100 pb-2"
													>
														<span className="text-stone-400 text-[10px] uppercase font-bold">
															{meta.label}
														</span>
														<span className="font-mono text-lg font-bold text-orange-600">
															{meta.format(val)}{" "}
															<span className="text-xs font-normal text-stone-400">
																{meta.unit}
															</span>
														</span>
													</div>
												);
											}

											return (
												<div key={key} className="flex flex-col">
													<span className="text-stone-400 text-[9px] font-bold">
														{meta?.label || key}
													</span>
													<span className="font-mono text-xs font-bold whitespace-nowrap">
														{meta?.format(val) || val.toFixed(4)}{" "}
														<span className="text-[9px] font-normal text-stone-400 ml-0.5">
															{meta?.unit}
														</span>
													</span>
												</div>
											);
										})}
								</div>
							</div>
						)}
					</div>
				)}
			</aside>
		);
	}

	// Edge Inspector
	if (selectedEdge && selectedEdge.data?.type === "thermal") {
		let probeTemp: number | null = null;
		let probeLabel = "Probe Temp";
		if (selectedEdge.data.result?.T_interface) {
			const targetX_manual = selectedEdge.data.probe_depth ?? 0;
			const tInt = selectedEdge.data.result.T_interface as number[];

			// Fallback layers for legacy/initial elements
			const layers = selectedEdge.data.layers || [
				{
					thickness: selectedEdge.data.thickness || 0.003,
					conductivity: selectedEdge.data.conductivity || 20.0,
				},
			];

			let targetX = targetX_manual;
			if (
				selectedEdge.data.probe_mode === "preset" &&
				selectedEdge.data.probe_preset
			) {
				const { type, index } = selectedEdge.data.probe_preset;
				// Ensure index is valid for current layers
				const safeIndex = Math.min(index, layers.length - 1);
				const layer = layers[safeIndex];

				let runningX = 0;
				for (let i = 0; i < safeIndex; i++) {
					runningX += layers[i].thickness;
				}

				const L = layer?.thickness || 0;
				if (type === "hot") {
					targetX = runningX;
					probeLabel = `L${safeIndex + 1} Hot Side`;
				} else if (type === "avg") {
					targetX = runningX + L / 2;
					probeLabel = `L${safeIndex + 1} Average`;
				} else {
					targetX = runningX + L;
					probeLabel = `L${safeIndex + 1} Cold Side`;
				}
			} else {
				const scale = unitPreferences.length === "mm" ? 1000 : 1;
				probeLabel = `Custom Depth d=${(targetX * scale).toFixed(unitPreferences.length === "mm" ? 1 : 3)}${unitPreferences.length}`;
			}

			if (tInt.length === layers.length + 1) {
				if (targetX <= 0) {
					probeTemp = tInt[0];
				} else {
					let found = false;
					let iterX = 0;
					for (let i = 0; i < layers.length; i++) {
						const t = layers[i].thickness;
						const nextX = iterX + t;
						if (targetX <= nextX + 1e-9) {
							const frac = t > 0 ? (targetX - iterX) / t : 0;
							probeTemp = tInt[i] + frac * (tInt[i + 1] - tInt[i]);
							found = true;
							break;
						}
						iterX = nextX;
					}
					if (!found) probeTemp = tInt[tInt.length - 1]; // Beyond cold side
				}
			} else {
				probeLabel = "(solve needed)";
				probeTemp = null;
			}
		}

		// Calculate derived HTCs for diagnostics
		const layers = selectedEdge.data.layers || [];
		const R_wall =
			layers.reduce(
				(acc: number, l: any) =>
					acc + (l.thickness || 0) / (l.conductivity || 1),
				0,
			) + (selectedEdge.data.R_fouling || 0);

		const h_wall = R_wall > 0 ? 1 / R_wall : Infinity;
		const h_a = selectedEdge.data.result?.h_a ?? Infinity;
		const h_b = selectedEdge.data.result?.h_b ?? Infinity;

		const h_total_inv =
			(h_a > 0 ? 1 / h_a : 0) + R_wall + (h_b > 0 ? 1 / h_b : 0);
		const h_total = h_total_inv > 0 ? 1 / h_total_inv : 0;

		return (
			<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto">
				<h2 className="text-lg font-bold border-b pb-2 uppercase text-orange-600">
					Thermal Wall
				</h2>

				{/* Wall Name Label */}
				<div className="flex flex-col gap-1">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Wall Name
					</label>
					<input
						type="text"
						value={selectedEdge.data.label ?? ""}
						onChange={(e) =>
							updateEdgeData(selectedEdge.id, { label: e.target.value })
						}
						placeholder="e.g. Custom Wall"
						className="p-1.5 border rounded text-xs w-full focus:ring-1 focus:ring-orange-500 outline-none"
					/>
				</div>

				{/* Show label on graph */}
				<label className="flex items-center gap-2 text-xs text-gray-600 cursor-pointer select-none">
					<input
						type="checkbox"
						checked={!!selectedEdge.data?.show_label}
						onChange={(e) =>
							updateEdgeData(selectedEdge.id, { show_label: e.target.checked })
						}
						className="rounded"
					/>
					Show label in graph
				</label>
				<AreaInput
					id={`area_edge_${selectedEdge.id}`}
					label="Contact Area"
					value={selectedEdge.data.area ?? 0.05}
					onChange={(val) => updateEdgeData(selectedEdge.id, { area: val })}
				/>

				<div className="flex flex-col gap-2">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Fouling Resistance (m²K/W)
					</label>
					<NumericInput
						value={selectedEdge.data.R_fouling ?? 0}
						onChange={(val) =>
							updateEdgeData(selectedEdge.id, { R_fouling: val })
						}
						className="p-1.5 border rounded text-sm bg-stone-50 focus:bg-white"
					/>
				</div>

				<div className="border-t pt-4">
					<WallLayersEditor
						layers={
							selectedEdge.data.layers || [
								{
									thickness: selectedEdge.data.thickness ?? 0.003,
									conductivity: selectedEdge.data.conductivity ?? 20.0,
									material: "generic",
								},
							]
						}
						onChange={(layers) => updateEdgeData(selectedEdge.id, { layers })}
					/>
				</div>

				<div className="flex flex-col gap-3 pt-4 border-t border-stone-100">
					<div className="flex flex-col gap-1">
						<label className="text-xs font-bold text-gray-500 uppercase">
							Probe Location
						</label>
						<select
							className="p-1.5 border rounded text-xs bg-white w-full border-stone-200 outline-none focus:ring-1 focus:ring-orange-500"
							value={
								selectedEdge.data.probe_mode === "preset" &&
								selectedEdge.data.probe_preset
									? `preset:${selectedEdge.data.probe_preset.type}:${selectedEdge.data.probe_preset.index}`
									: "custom"
							}
							onChange={(e) => {
								const val = e.target.value;
								if (val === "custom") {
									updateEdgeData(selectedEdge.id, { probe_mode: "custom" });
								} else {
									const [_, type, idx] = val.split(":");
									updateEdgeData(selectedEdge.id, {
										probe_mode: "preset",
										probe_preset: { type, index: Number.parseInt(idx, 10) },
									});
								}
							}}
						>
							<option value="custom">Custom Position (d=x)</option>
							{/* Predefined locations for Layer 1 */}
							<option value="preset:hot:0">L1 Hot Side</option>
							<option value="preset:avg:0">L1 Average</option>
							<option value="preset:cold:0">L1 Cold Side</option>
							{/* Dynamic locations for additional layers */}
							{(selectedEdge.data.layers || [])
								.slice(1)
								.map((_: any, i: number) => (
									<React.Fragment key={i + 1}>
										<option value={`preset:avg:${i + 1}`}>
											L{i + 2} Average
										</option>
										<option value={`preset:cold:${i + 1}`}>
											L{i + 2} Cold Side
										</option>
									</React.Fragment>
								))}
						</select>
					</div>

					{selectedEdge.data.probe_mode === "custom" && (
						<LengthInput
							id={`probe_depth_${selectedEdge.id}`}
							label="Custom Depth d"
							value={selectedEdge.data.probe_depth ?? 0}
							onChange={(val) => {
								const totalThickness =
									(selectedEdge.data.layers || []).reduce(
										(acc: number, l: any) => acc + (l.thickness || 0),
										0,
									) || 0.003;
								const clamped = Math.min(Math.max(0, val), totalThickness);
								updateEdgeData(selectedEdge.id, { probe_depth: clamped });
							}}
							placeholder="0.00"
						/>
					)}
					<span className="text-[9px] text-gray-400 font-normal italic text-right -mt-1">
						0 = Hot, L = Cold
					</span>
				</div>

				{selectedEdge.data.result && (
					<div className="mt-4 border-t pt-4 flex flex-col gap-3">
						<h3 className="text-xs font-bold text-orange-600 uppercase">
							Solve Results
						</h3>
						<div className="bg-stone-50 rounded p-3 flex flex-col gap-2 border border-stone-100 shadow-sm">
							<div className="flex justify-between items-center">
								<span className="text-[10px] text-stone-500 font-bold uppercase">
									Heat Flow (Q)
								</span>
								<span className="font-mono text-xs font-black text-orange-700">
									{selectedEdge.data.result.Q?.toLocaleString()} W
								</span>
							</div>
							<div className="flex justify-between items-center">
								<span className="text-[10px] text-stone-500 font-bold uppercase">
									Hot Side Wall Temp
								</span>
								<span className="font-mono text-xs font-black">
									{selectedEdge.data.result.T_hot?.toFixed(1)} K
								</span>
							</div>
							<div className="grid grid-cols-4 gap-1 pt-2 border-t border-stone-100">
								<div className="flex flex-col">
									<span className="text-[7px] text-stone-400 uppercase">
										h (Side A)
									</span>
									<span className="text-[10px] font-mono font-bold">
										{h_a < 1e6 ? h_a.toFixed(1) : "∞"}
									</span>
								</div>
								<div className="flex flex-col">
									<span className="text-[7px] text-stone-400 uppercase">
										h (Wall)
									</span>
									<span className="text-[10px] font-mono font-bold text-orange-600">
										{h_wall < 1e6 ? h_wall.toFixed(1) : "∞"}
									</span>
								</div>
								<div className="flex flex-col">
									<span className="text-[7px] text-stone-400 uppercase">
										h (Side B)
									</span>
									<span className="text-[10px] font-mono font-bold">
										{h_b < 1e6 ? h_b.toFixed(1) : "∞"}
									</span>
								</div>
								<div className="flex flex-col">
									<span className="text-[7px] text-stone-400 uppercase">
										h (Total)
									</span>
									<span className="text-[10px] font-mono font-bold text-blue-600">
										{h_total.toFixed(1)}
									</span>
								</div>
							</div>
							{probeTemp !== null && (
								<div className="flex justify-between items-center pt-2 mt-2 border-t border-stone-200 bg-orange-50/50 p-2 rounded border">
									<span className="text-[10px] text-orange-600 font-bold uppercase tracking-wider">
										{probeLabel}
									</span>
									<span className="font-mono text-sm font-black text-orange-600">
										{probeTemp.toFixed(1)} K
									</span>
								</div>
							)}
						</div>
					</div>
				)}
			</aside>
		);
	}

	if (selectedEdge) {
		return (
			<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto">
				<h2 className="text-lg font-bold border-b pb-2 uppercase text-slate-600">
					Connection
				</h2>

				{/* Custom label */}
				<div className="flex flex-col gap-1">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Custom Label (optional)
					</label>
					<input
						type="text"
						value={selectedEdge.data?.label ?? ""}
						onChange={(e) =>
							updateEdgeData(selectedEdge.id, {
								label: e.target.value || undefined,
							})
						}
						placeholder="e.g. Main Feed"
						className="p-1.5 border rounded text-xs w-full focus:ring-1 focus:ring-slate-400 outline-none"
					/>
				</div>

				{/* Show label on graph */}
				<label className="flex items-center gap-2 text-xs text-gray-600 cursor-pointer select-none">
					<input
						type="checkbox"
						checked={!!selectedEdge.data?.show_label}
						onChange={(e) =>
							updateEdgeData(selectedEdge.id, { show_label: e.target.checked })
						}
						className="rounded"
					/>
					Show label in graph
				</label>

				{/* Auto ID for reference */}
				<div className="flex flex-col gap-1 mt-2 border-t pt-3">
					<span className="text-[10px] font-bold text-stone-400 uppercase">
						Internal ID
					</span>
					<span className="font-mono text-[10px] text-stone-400 break-all">
						{selectedEdge.id}
					</span>
				</div>
			</aside>
		);
	}

	return null;
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

	// Type-aware guess field selection. Elements carry only a solvable m_dot
	// unknown; nodes carry P / T. Pressure-boundary nodes are fully specified
	// (no unknowns) and are handled by returning an empty field list below.
	const FIELDS_BY_TYPE: Record<string, { key: string; unit: string }[]> = {
		// Volume nodes (P and T are solver unknowns)
		plenum: [
			{ key: "P", unit: "Pa" },
			{ key: "T", unit: "K" },
		],
		combustor: [
			{ key: "P", unit: "Pa" },
			{ key: "T", unit: "K" },
		],
		momentum_chamber: [
			{ key: "P", unit: "Pa" },
			{ key: "T", unit: "K" },
		],
		// Boundaries: P_total is solved (mass-flow inlet) or fixed (pressure
		// inlet). Short key "P" is broadcast to both by the backend; we show
		// it only where at least one of the two is actually an unknown.
		mass_boundary: [{ key: "P", unit: "Pa" }],
		// Elements (edges in React Flow but still nodes here): m_dot only.
		channel: [{ key: "m_dot", unit: "kg/s" }],
		orifice: [{ key: "m_dot", unit: "kg/s" }],
		discrete_loss: [{ key: "m_dot", unit: "kg/s" }],
	};
	const guessFields = FIELDS_BY_TYPE[node.type as string] ?? [];
	if (guessFields.length === 0) {
		return null;
	}

	return (
		<div className="mt-4 border-t pt-4">
			<h3 className="text-xs font-bold text-stone-500 uppercase mb-2">
				Initial Guess Overrides
			</h3>
			<div className="flex flex-col gap-2">
				{guessFields.map(({ key, unit }) => (
					<div key={key} className="flex items-center justify-between gap-2">
						<span className="text-[10px] font-mono text-gray-500">
							{key}{" "}
							<span className="text-[9px] text-gray-400 font-normal">
								[{unit}]
							</span>
						</span>
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
