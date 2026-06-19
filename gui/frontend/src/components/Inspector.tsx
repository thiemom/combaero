import { Activity, Crosshair, RotateCw } from "lucide-react";
import React from "react";
import useStore from "../store/useStore";
import { resolveField } from "../utils/diagnostics";
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
		setRotSpeedUnit,
		isExporting,
		exportNetworkResults,
		displaySettings,
		solverSettings,
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
			const unit = QUANTITY_CATALOGUE[k]?.unit;
			availableKeys.push({
				key: fullKey,
				label: unit ? `${k} [${unit}]` : k,
			});
		};

		if (result) {
			// Only offer keys whose value is actually a finite number in this
			// target's result — mirrors discoverFields to avoid NaN / Infinity slots.
			const isFiniteNum = (v: unknown): v is number =>
				typeof v === "number" && Number.isFinite(v);
			if (result.state) {
				for (const k of Object.keys(result.state)) {
					if (isFiniteNum(result.state[k])) addKey(k, true);
				}
			}
			for (const k of Object.keys(result)) {
				if (k !== "state" && isFiniteNum(result[k])) addKey(k);
			}
		} else {
			// Fallback: show common state and transport quantities if no result yet
			for (const k of [...BASIC_STATE_KEYS, ...TRANSPORT_KEYS]) {
				addKey(k);
			}
			// Specifically add area change ones if target is area change
			if (targetElement?.type === "area_change") {
				addKey("zeta");
				addKey("mach_loss_ref");
				addKey("mach_in");
				addKey("mach_out");
				addKey("area_ratio");
			}
		}

		// Sort by label
		availableKeys.sort((a, b) => a.label.localeCompare(b.label));

		// availableKeys is already filtered to finite values for this target
		// (via isFiniteNum when building the list). No further displaySettings
		// filtering here — that would hide fields like phi/theta that exist on
		// the target but aren't in the current displaySettings selection.

		const slotSelect = (
			slot: "slot1_key" | "slot2_key",
			current: string | undefined,
		) => (
			<select
				className="p-1.5 border rounded text-xs bg-white w-full"
				value={current ?? ""}
				onChange={(e) =>
					selectedNode &&
					updateNodeData(
						selectedNode.id,
						{ [slot]: e.target.value || undefined },
						true,
					)
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
							updateNodeData(selectedNode.id, { label: e.target.value }, true)
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
								updateNodeData(
									selectedNode.id,
									{ rotation: (currentRotation + 90) % 360 },
									true,
								);
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
							updateNodeData(selectedNode.id, { label: e.target.value }, true)
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
								htmlFor={`Tt_mb_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Tt (K)
							</label>
							<NumericInput
								id={`Tt_mb_${selectedNode.id}`}
								value={selectedNode.data.Tt ?? 300}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										Tt: val,
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
							label="Pt"
							value={selectedNode.data.Pt ?? 101325}
							onChange={(val) => updateNodeData(selectedNode.id, { Pt: val })}
						/>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`Tt_pb_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Tt (K)
							</label>
							<NumericInput
								id={`Tt_pb_${selectedNode.id}`}
								value={selectedNode.data.Tt ?? 300}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										Tt: val,
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

				{selectedNode.type === "wall" && (
					<div className="text-xs text-stone-400 italic px-1">
						Closed end — zero mass flow. No configurable parameters.
					</div>
				)}

				{selectedNode.type === "channel" &&
					(() => {
						const upstreamEdge = edges.find(
							(e) => e.target === selectedNode.id && !e.data?.type,
						);
						const downstreamEdge = edges.find(
							(e) => e.source === selectedNode.id && !e.data?.type,
						);
						const upstreamNode = upstreamEdge
							? nodes.find((n) => n.id === upstreamEdge.source)
							: undefined;
						const downstreamNode = downstreamEdge
							? nodes.find((n) => n.id === downstreamEdge.target)
							: undefined;
						// Returns {d, dh} to propagate both circular-equivalent diameter
						// and hydraulic diameter through chains. channel: dh = Dh ?? D.
						// area_change: dh = D_h (if set) ?? sqrt(4F/pi). Tees: dh = d.
						const walkGeometry = (
							n: typeof upstreamNode,
							edge: typeof upstreamEdge,
							side: "upstream" | "downstream",
							depth = 0,
						): { d: number; dh: number } | undefined => {
							if (!n || depth > 20) return undefined;
							if (n.type === "channel") {
								const d = n.data.D as number | null | undefined;
								const dh = n.data.Dh as number | null | undefined;
								if (d != null) return { d, dh: dh ?? d };
								const hopEdge =
									side === "upstream"
										? edges.find((e) => e.target === n.id && !e.data?.type)
										: edges.find((e) => e.source === n.id && !e.data?.type);
								const hopNode = hopEdge
									? nodes.find(
											(nn) =>
												nn.id ===
												(side === "upstream" ? hopEdge.source : hopEdge.target),
										)
									: undefined;
								return walkGeometry(hopNode, hopEdge, side, depth + 1);
							}
							if (n.type === "area_change") {
								const dhAc = n.data.D_h as number | null | undefined;
								const f =
									side === "upstream"
										? (n.data.F1 as number | null | undefined)
										: (n.data.F0 as number | null | undefined);
								if (f != null) {
									const d = Math.sqrt((4 * f) / Math.PI);
									return { d, dh: dhAc != null && dhAc > 0 ? dhAc : d };
								}
								return undefined;
							}
							if (n.type === "tee_junction") {
								const handle =
									side === "upstream" ? edge?.sourceHandle : edge?.targetHandle;
								const isBranch =
									side === "upstream"
										? handle === "port-branch-source"
										: handle === "port-branch-target";
								let d: number | undefined;
								if (isBranch) {
									const fb = n.data.F_branch as number | null | undefined;
									const fc = n.data.F_C as number | null | undefined;
									const psi = (n.data.psi as number | undefined) ?? 1.0;
									const fBranch =
										fb != null ? fb : fc != null ? fc / psi : undefined;
									d =
										fBranch != null
											? Math.sqrt((4 * fBranch) / Math.PI)
											: undefined;
								} else {
									const fc = n.data.F_C as number | null | undefined;
									d = fc != null ? Math.sqrt((4 * fc) / Math.PI) : undefined;
								}
								return d != null ? { d, dh: d } : undefined;
							}
							return undefined;
						};
						const inheritedGeomUp = walkGeometry(
							upstreamNode,
							upstreamEdge,
							"upstream",
						);
						const inheritedGeomDown = walkGeometry(
							downstreamNode,
							downstreamEdge,
							"downstream",
						);
						const inheritedD = inheritedGeomUp?.d ?? inheritedGeomDown?.d;
						const inheritedDh = inheritedGeomUp?.dh ?? inheritedGeomDown?.dh;
						const inheritSide =
							inheritedGeomUp != null
								? "upstream"
								: inheritedGeomDown != null
									? "downstream"
									: undefined;
						const inheritSourceId =
							inheritSide === "upstream"
								? upstreamNode?.id
								: downstreamNode?.id;
						const userSetD =
							selectedNode.data.D !== null && selectedNode.data.D !== undefined;
						const userSetDh =
							selectedNode.data.Dh !== null &&
							selectedNode.data.Dh !== undefined;
						const effectiveD = userSetD
							? selectedNode.data.D
							: (inheritedD ?? 0.1);
						// True when the upstream source provides a non-circular Dh
						const dhIsNonCircular =
							inheritedDh != null &&
							(inheritedD == null ||
								Math.abs(inheritedDh - inheritedD) > 1e-10);
						return (
							<>
								<LengthInput
									id={`L_${selectedNode.id}`}
									label="Length Flow Path"
									value={selectedNode.data.L ?? 1.0}
									onChange={(val) =>
										updateNodeData(selectedNode.id, { L: val })
									}
								/>
								<div className="flex flex-col gap-2">
									<div className="flex items-center justify-between">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Channel Inner Diameter
										</label>
										{userSetD && (
											<button
												type="button"
												className="text-[9px] text-blue-500 hover:underline"
												onClick={() =>
													updateNodeData(selectedNode.id, { D: null })
												}
											>
												Reset to inherited
											</button>
										)}
									</div>
									<LengthInput
										id={`D_${selectedNode.id}`}
										label=""
										value={userSetD ? selectedNode.data.D : (inheritedD ?? 0.1)}
										placeholder={inheritedD ?? 0.1}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { D: val })
										}
										onClear={() => updateNodeData(selectedNode.id, { D: null })}
									/>
									{!userSetD && inheritedD != null && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Inherited from {inheritSide} ({inheritSourceId}). Edit to
											override.
										</p>
									)}
								</div>
								<div className="flex flex-col gap-2">
									<div className="flex items-center justify-between">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Hydraulic Diameter Dh
										</label>
										{userSetDh && (
											<button
												type="button"
												className="text-[9px] text-blue-500 hover:underline"
												onClick={() =>
													updateNodeData(selectedNode.id, { Dh: null })
												}
											>
												{dhIsNonCircular
													? "Reset to inherited"
													: "Reset to circular"}
											</button>
										)}
									</div>
									<LengthInput
										id={`Dh_${selectedNode.id}`}
										label=""
										value={
											userSetDh
												? selectedNode.data.Dh
												: (inheritedDh ?? effectiveD)
										}
										placeholder={inheritedDh ?? effectiveD}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { Dh: val })
										}
										onClear={() =>
											updateNodeData(selectedNode.id, { Dh: null })
										}
										inputClassName={
											userSetDh ? "font-semibold" : "text-gray-400"
										}
									/>
									{!userSetDh && dhIsNonCircular && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Inherited non-circular Dh from {inheritSide}. Edit to
											override.
										</p>
									)}
									{!userSetDh && !dhIsNonCircular && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Dh = D (circular). Override for non-circular ducts.
										</p>
									)}
								</div>
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Friction Model
									</label>
									<select
										className="p-2 border rounded bg-white text-xs"
										value={selectedNode.data.friction_model || "haaland"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												friction_model: e.target.value,
											})
										}
									>
										<option value="haaland">
											Haaland (default, ~2–3% vs Colebrook)
										</option>
										<option value="serghides">
											Serghides (&lt;0.3% vs Colebrook)
										</option>
										<option value="colebrook">
											Colebrook-White (reference)
										</option>
										<option value="petukhov">
											Petukhov (smooth pipe, pairs with Petukhov Nu)
										</option>
									</select>
								</div>
								<LengthInput
									id={`roughness_${selectedNode.id}`}
									label="Surface Roughness"
									value={selectedNode.data.roughness ?? 1e-5}
									onChange={(val) =>
										updateNodeData(selectedNode.id, { roughness: val })
									}
									disabled={selectedNode.data.friction_model === "petukhov"}
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
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.Nu_multiplier ?? 1}`}
										</span>
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
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.f_multiplier ?? 1}`}
										</span>
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
											updateNodeData(selectedNode.id, {
												regime: e.target.value,
											})
										}
									>
										<option value="default">Default (Global)</option>
										<option value="incompressible">
											Forced Incompressible
										</option>
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
						);
					})()}

				{selectedNode.type === "orifice" &&
					(() => {
						const upstreamEdge = edges.find(
							(e) => e.target === selectedNode.id && !e.data?.type,
						);
						const upstreamNode = upstreamEdge
							? nodes.find((n) => n.id === upstreamEdge.source)
							: undefined;
						const upstreamChannelD: number | null | undefined =
							upstreamNode?.type === "channel"
								? (upstreamNode.data.D as number | null | undefined)
								: undefined;
						const boreDiameter: number = selectedNode.data.diameter ?? 0.08;
						const boreArea = (Math.PI / 4) * boreDiameter * boreDiameter;
						return (
							<>
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Bore Diameter
									</label>
									<LengthInput
										id={`diameter_${selectedNode.id}`}
										label=""
										value={boreDiameter}
										placeholder={0.08}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { diameter: val })
										}
										onClear={() =>
											updateNodeData(selectedNode.id, { diameter: null })
										}
									/>
									{(() => {
										const pipeD =
											upstreamChannelD != null && upstreamChannelD > 0
												? upstreamChannelD
												: null;
										const beta =
											pipeD != null && boreDiameter < pipeD
												? boreDiameter / pipeD
												: null;
										const beta2 = beta != null ? beta ** 2 : null;
										const E =
											beta != null ? 1.0 / Math.sqrt(1 - beta ** 4) : null;
										return (
											<div className="text-[9px] text-gray-400 -mt-1 flex flex-col gap-0.5">
												<span>
													Bore area:{" "}
													<span className="font-mono">
														{boreArea.toExponential(3)} m²
													</span>
													{beta2 != null && (
														<>
															{" · "}Area ratio β²:{" "}
															<span className="font-mono">
																{beta2.toFixed(4)}
															</span>
														</>
													)}
												</span>
												{pipeD != null && (
													<span>
														Pipe D:{" "}
														<span className="font-mono">
															{pipeD.toFixed(4)} m
														</span>
														{beta != null && E != null ? (
															<>
																{" · "}β = {beta.toFixed(4)} → E ={" "}
																<span className="font-mono">
																	{E.toFixed(4)}
																</span>{" "}
																(velocity-of-approach)
															</>
														) : (
															<span className="text-amber-500">
																{" · "}β ≥ 1 — correction skipped
															</span>
														)}
													</span>
												)}
											</div>
										);
									})()}
								</div>

								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Discharge Model (Cd)
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={
											selectedNode.data.correlation || "ReaderHarrisGallagher"
										}
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
											updateNodeData(selectedNode.id, {
												regime: e.target.value,
											})
										}
									>
										<option value="default">Default (Global)</option>
										<option value="incompressible">
											Forced Incompressible
										</option>
										<option value="compressible">Forced Compressible</option>
									</select>
								</div>
								<InitialGuessEditor node={selectedNode} />
							</>
						);
					})()}

				{selectedNode.type === "area_change" &&
					(() => {
						const upstreamEdge = edges.find(
							(e) => e.target === selectedNode.id && !e.data?.type,
						);
						const downstreamEdge = edges.find(
							(e) => e.source === selectedNode.id && !e.data?.type,
						);
						const upstreamNode = upstreamEdge
							? nodes.find((n) => n.id === upstreamEdge.source)
							: undefined;
						const downstreamNode = downstreamEdge
							? nodes.find((n) => n.id === downstreamEdge.target)
							: undefined;
						const areaFromD = (d: number | null | undefined) =>
							d != null ? (Math.PI / 4) * d * d : undefined;
						const walkArea = (
							nodeId: string,
							dir: "upstream" | "downstream",
							depth = 0,
						): number | undefined => {
							if (depth > 20) return undefined;
							const e =
								dir === "upstream"
									? edges.find((ee) => ee.target === nodeId && !ee.data?.type)
									: edges.find((ee) => ee.source === nodeId && !ee.data?.type);
							const n = e
								? nodes.find(
										(nn) =>
											nn.id === (dir === "upstream" ? e.source : e.target),
									)
								: undefined;
							if (!n) return undefined;
							if (n.type === "channel") {
								const d = n.data.D as number | null | undefined;
								if (d != null) return areaFromD(d);
								return walkArea(n.id, dir, depth + 1);
							}
							// For area_change: walking upstream reads its exit face (F1);
							// walking downstream reads its entry face (F0).
							if (n.type === "area_change") {
								const f = (dir === "upstream" ? n.data.F1 : n.data.F0) as
									| number
									| null
									| undefined;
								if (f != null && f > 0) return f;
								return walkArea(n.id, dir, depth + 1);
							}
							if (n.type === "tee_junction") {
								const fc = n.data.F_C as number | null | undefined;
								if (fc != null && fc > 0) return fc;
							}
							const a = n.data?.area as number | null | undefined;
							return a != null && a > 0 ? a : undefined;
						};
						const inheritedF0: number | undefined = walkArea(
							selectedNode.id,
							"upstream",
						);
						const inheritedF1: number | undefined = walkArea(
							selectedNode.id,
							"downstream",
						);
						const userSetF0 =
							selectedNode.data.F0 !== null &&
							selectedNode.data.F0 !== undefined;
						const userSetF1 =
							selectedNode.data.F1 !== null &&
							selectedNode.data.F1 !== undefined;
						const effectiveF0 = userSetF0
							? selectedNode.data.F0
							: (inheritedF0 ?? 0.01);
						const effectiveF1 = userSetF1
							? selectedNode.data.F1
							: (inheritedF1 ?? 0.02);
						return (
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

								<div className="flex flex-col gap-2">
									<div className="flex items-center justify-between">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Upstream Area (F0)
										</label>
										{userSetF0 && (
											<button
												type="button"
												className="text-[9px] text-blue-500 hover:underline"
												onClick={() =>
													updateNodeData(selectedNode.id, { F0: null })
												}
											>
												Reset to inherited
											</button>
										)}
									</div>
									<AreaInput
										id={`F0_${selectedNode.id}`}
										label=""
										value={effectiveF0}
										placeholder={
											inheritedF0 != null ? inheritedF0.toFixed(5) : "0.01"
										}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { F0: val })
										}
										onClear={() =>
											updateNodeData(selectedNode.id, { F0: null })
										}
										inputClassName={
											userSetF0 ? "font-semibold" : "text-gray-400"
										}
									/>
									{!userSetF0 && inheritedF0 != null && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Inherited from upstream ({upstreamNode?.id}). Edit to
											override.
										</p>
									)}
								</div>

								<div className="flex flex-col gap-2">
									<div className="flex items-center justify-between">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Downstream Area (F1)
										</label>
										{userSetF1 && (
											<button
												type="button"
												className="text-[9px] text-blue-500 hover:underline"
												onClick={() =>
													updateNodeData(selectedNode.id, { F1: null })
												}
											>
												Reset to inherited
											</button>
										)}
									</div>
									<AreaInput
										id={`F1_${selectedNode.id}`}
										label=""
										value={effectiveF1}
										placeholder={
											inheritedF1 != null ? inheritedF1.toFixed(5) : "0.02"
										}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { F1: val })
										}
										onClear={() =>
											updateNodeData(selectedNode.id, { F1: null })
										}
										inputClassName={
											userSetF1 ? "font-semibold" : "text-gray-400"
										}
									/>
									{!userSetF1 && inheritedF1 != null && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Inherited from downstream ({downstreamNode?.id}). Edit to
											override.
										</p>
									)}
								</div>

								{(() => {
									const ratio =
										effectiveF0 > 0 ? effectiveF1 / effectiveF0 : null;
									const d0 = Math.sqrt((4 * effectiveF0) / Math.PI);
									const d1 = Math.sqrt((4 * effectiveF1) / Math.PI);
									const label =
										ratio == null
											? null
											: ratio > 1.005
												? "expansion"
												: ratio < 0.995
													? "contraction"
													: "straight";
									return ratio != null ? (
										<div className="rounded bg-stone-50 border border-stone-200 px-3 py-2 text-[10px] text-stone-600 flex flex-col gap-0.5">
											<div className="flex justify-between">
												<span className="font-medium text-stone-500 uppercase tracking-wide text-[8px]">
													Area ratio
												</span>
												<span
													className={
														label === "expansion"
															? "text-blue-600 font-semibold"
															: label === "contraction"
																? "text-amber-600 font-semibold"
																: "text-stone-500"
													}
												>
													{label}
												</span>
											</div>
											<div className="flex justify-between mt-0.5">
												<span className="text-stone-400">F₁/F₀</span>
												<span className="font-mono font-semibold">
													{ratio.toFixed(4)}
												</span>
											</div>
											<div className="flex justify-between">
												<span className="text-stone-400">D₀ → D₁</span>
												<span className="font-mono">
													{(
														d0 * (unitPreferences.length === "mm" ? 1000 : 1)
													).toFixed(1)}{" "}
													→{" "}
													{(
														d1 * (unitPreferences.length === "mm" ? 1000 : 1)
													).toFixed(1)}{" "}
													{unitPreferences.length}
												</span>
											</div>
										</div>
									) : null;
								})()}

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
											Defaults to circular diameter:{" "}
											{(
												Math.sqrt(
													(4 * Math.min(effectiveF0, effectiveF1)) / Math.PI,
												) * (unitPreferences.length === "mm" ? 1000 : 1)
											).toFixed(2)}{" "}
											{unitPreferences.length}
										</p>
									</div>
								)}

								<InitialGuessEditor node={selectedNode} />
							</>
						);
					})()}

				{selectedNode.type === "tee_junction" &&
					(() => {
						const teeType = selectedNode.data.tee_type ?? "merging";

						// Helper: find the channel node connected to a named tee port.
						// For merging: common port is outgoing, straight/branch are incoming.
						// For branching: common port is incoming, straight/branch are outgoing.
						const findPortChannel = (portName: string) => {
							const isOutgoing =
								(teeType === "merging" && portName === "common") ||
								(teeType === "branching" && portName !== "common");
							let neighbour: (typeof nodes)[0] | undefined;
							if (isOutgoing) {
								const edge = edges.find(
									(e) =>
										e.source === selectedNode.id &&
										e.sourceHandle === `port-${portName}-source` &&
										!e.data?.type,
								);
								neighbour = edge
									? nodes.find((n) => n.id === edge.target)
									: undefined;
							} else {
								const edge = edges.find(
									(e) =>
										e.target === selectedNode.id &&
										e.targetHandle === `port-${portName}-target` &&
										!e.data?.type,
								);
								neighbour = edge
									? nodes.find((n) => n.id === edge.source)
									: undefined;
							}
							if (neighbour?.type === "channel") return neighbour;
							// Neighbour is a plenum/boundary; look one hop further
							if (isOutgoing) {
								return nodes.find(
									(n) =>
										n.type === "channel" &&
										edges.some(
											(e) =>
												e.source === neighbour?.id &&
												e.target === n.id &&
												!e.data?.type,
										),
								);
							}
							return nodes.find(
								(n) =>
									n.type === "channel" &&
									edges.some(
										(e) =>
											e.target === neighbour?.id &&
											e.source === n.id &&
											!e.data?.type,
									),
							);
						};

						// Walk through null-D channels (and through tee/area-change nodes)
						// to find the first explicit D in the chain.
						const walkChD = (
							ch: (typeof nodes)[0] | undefined,
							isOutgoing: boolean,
							depth = 0,
						): number | undefined => {
							if (!ch || depth > 20) return undefined;
							const d = ch.data.D as number | null | undefined;
							if (d != null) return d;
							const e = isOutgoing
								? edges.find((ee) => ee.source === ch.id && !ee.data?.type)
								: edges.find((ee) => ee.target === ch.id && !ee.data?.type);
							if (!e) return undefined;
							const next = nodes.find(
								(n) => n.id === (isOutgoing ? e.target : e.source),
							);
							if (!next) return undefined;
							if (next.type === "channel")
								return walkChD(next, isOutgoing, depth + 1);
							// Non-channel: extract area-equivalent D from tee/area-change
							if (next.type === "tee_junction") {
								const fc = next.data.F_C as number | null | undefined;
								if (fc != null && fc > 0) return Math.sqrt((4 * fc) / Math.PI);
							}
							if (next.type === "area_change") {
								const f = isOutgoing
									? (next.data.F0 as number | null | undefined)
									: (next.data.F1 as number | null | undefined);
								if (f != null && f > 0) return Math.sqrt((4 * f) / Math.PI);
							}
							return undefined;
						};

						const portD = (portName: string) => {
							const isOutgoing =
								(teeType === "merging" && portName === "common") ||
								(teeType === "branching" && portName !== "common");
							return walkChD(findPortChannel(portName), isOutgoing);
						};

						const commonCh = findPortChannel("common");
						const straightCh = findPortChannel("straight");
						const branchCh = findPortChannel("branch");

						const areaFromD = (d: number | null | undefined) =>
							d != null ? (Math.PI / 4) * d * d : undefined;

						// F_C: common arm (= straight arm per Bassett); prefer common over straight channel
						const inheritedFC: number | undefined =
							areaFromD(portD("common")) ?? areaFromD(portD("straight"));
						const inheritedFBranch: number | undefined = areaFromD(
							portD("branch"),
						);

						const userSetFC = selectedNode.data.F_C != null;
						const userSetFBranch = selectedNode.data.F_branch != null;

						const effectiveFC: number = userSetFC
							? (selectedNode.data.F_C as number)
							: (inheritedFC ?? 0.01);
						const effectiveFBranch: number = userSetFBranch
							? (selectedNode.data.F_branch as number)
							: (inheritedFBranch ?? effectiveFC);
						const psi =
							effectiveFBranch > 0 ? effectiveFC / effectiveFBranch : null;
						const dFromF = (f: number) => Math.sqrt((4 * f) / Math.PI);
						const lu = unitPreferences.length === "mm" ? 1000 : 1;
						const lu_label = unitPreferences.length;
						return (
							<div className="flex flex-col gap-4">
								{/* Port legend */}
								{(() => {
									const tm = teeType === "merging";
									return (
										<div className="rounded border border-stone-200 bg-stone-50 px-3 py-2 text-[9px] leading-snug text-stone-600">
											<div className="font-bold uppercase text-stone-400 mb-1">
												Port Guide
											</div>
											<div className="grid grid-cols-3 gap-1 text-center font-mono mb-1.5">
												<div>
													<span
														className="font-extrabold"
														style={{ color: tm ? "#3b82f6" : "#f59e0b" }}
													>
														S
													</span>{" "}
													straight
												</div>
												<div>
													<span
														className="font-extrabold"
														style={{ color: tm ? "#f59e0b" : "#3b82f6" }}
													>
														C
													</span>{" "}
													common
												</div>
												<div>
													<span
														className="font-extrabold"
														style={{ color: tm ? "#3b82f6" : "#f59e0b" }}
													>
														B
													</span>{" "}
													branch
												</div>
											</div>
											<div className="text-stone-500">
												{tm ? (
													<>
														<span style={{ color: "#3b82f6" }}>S, B</span> →{" "}
														<span style={{ color: "#f59e0b" }}>C</span>
														{"  "}(two inlets, one outlet at C)
													</>
												) : (
													<>
														<span style={{ color: "#3b82f6" }}>C</span> →{" "}
														<span style={{ color: "#f59e0b" }}>S, B</span>
														{"  "}(one inlet at C, two outlets)
													</>
												)}
											</div>
										</div>
									);
								})()}

								{/* Tee type */}
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Tee Type
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={selectedNode.data.tee_type ?? "merging"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												tee_type: e.target.value,
											})
										}
									>
										<option value="merging">
											Merging (2 inlets, 1 outlet)
										</option>
										<option value="branching">
											Branching (1 inlet, 2 outlets)
										</option>
									</select>
									<p className="text-[9px] text-gray-400 italic">
										Declares intended flow direction; the Bassett model blends
										smoothly if flow reverses.
									</p>
								</div>

								{/* Branch angle */}
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Branch Angle (deg)
									</label>
									<NumericInput
										value={selectedNode.data.theta_deg ?? 90}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { theta_deg: val })
										}
										className="p-2 border rounded"
										placeholder="90"
									/>
									<p className="text-[9px] text-gray-400 italic">
										Angle between branch and common arm. Sign is ignored
										(symmetric). Bassett (2001) validated range: 30–90 deg.
									</p>
								</div>

								{/* ── Arm areas (3-arm, mirrors AreaChange) ── */}

								{/* Common arm (F_C) */}
								<div className="flex flex-col gap-2">
									<div className="flex items-center justify-between">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Common Arm Area (F_C)
										</label>
										{userSetFC && (
											<button
												type="button"
												className="text-[9px] text-blue-500 hover:underline"
												onClick={() =>
													updateNodeData(selectedNode.id, { F_C: null })
												}
											>
												Reset to inherited
											</button>
										)}
									</div>
									<AreaInput
										id={`F_C_${selectedNode.id}`}
										label=""
										value={effectiveFC}
										placeholder={
											inheritedFC != null ? inheritedFC.toFixed(5) : "0.01"
										}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { F_C: val })
										}
										onClear={() =>
											updateNodeData(selectedNode.id, { F_C: null })
										}
										inputClassName={
											userSetFC ? "font-semibold" : "text-gray-400"
										}
									/>
									{!userSetFC && inheritedFC != null && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Inherited from{" "}
											{commonCh?.id ?? straightCh?.id ?? "connected channel"}.
											Edit to override.
										</p>
									)}
								</div>

								{/* Straight arm (read-only = F_C per Bassett) */}
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Straight Arm Area (F_S)
									</label>
									<div className="rounded bg-stone-50 border border-stone-200 px-3 py-2 text-[10px] text-stone-500 flex justify-between items-center">
										<span className="font-mono">
											{effectiveFC.toExponential(3)} m²
											{"  "}(D = {(dFromF(effectiveFC) * lu).toFixed(2)}{" "}
											{lu_label})
										</span>
										<span className="italic text-[9px] text-stone-400">
											= F_C (Bassett)
										</span>
									</div>
									<p className="text-[9px] text-stone-400 italic">
										The Bassett (2001) model assumes F_S = F_C. If the physical
										straight arm differs, add an AreaChange element.
									</p>
								</div>

								{/* Branch arm (F_branch) */}
								<div className="flex flex-col gap-2">
									<div className="flex items-center justify-between">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Branch Arm Area (F_B)
										</label>
										{userSetFBranch && (
											<button
												type="button"
												className="text-[9px] text-blue-500 hover:underline"
												onClick={() =>
													updateNodeData(selectedNode.id, {
														F_branch: null,
													})
												}
											>
												Reset to inherited
											</button>
										)}
									</div>
									<AreaInput
										id={`F_branch_${selectedNode.id}`}
										label=""
										value={effectiveFBranch}
										placeholder={
											inheritedFBranch != null
												? inheritedFBranch.toFixed(5)
												: inheritedFC != null
													? inheritedFC.toFixed(5)
													: "0.01"
										}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { F_branch: val })
										}
										onClear={() =>
											updateNodeData(selectedNode.id, { F_branch: null })
										}
										inputClassName={
											userSetFBranch ? "font-semibold" : "text-gray-400"
										}
									/>
									{!userSetFBranch && inheritedFBranch != null && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Inherited from {branchCh?.id}. Edit to override.
										</p>
									)}
								</div>

								{/* Derived ratio summary card */}
								{psi != null && (
									<div className="rounded bg-stone-50 border border-stone-200 px-3 py-2 text-[10px] text-stone-600 flex flex-col gap-0.5">
										<div className="flex justify-between">
											<span className="font-medium text-stone-500 uppercase tracking-wide text-[8px]">
												Area ratio
											</span>
											<span
												className={
													psi > 1.005
														? "text-amber-600 font-semibold"
														: psi < 0.995
															? "text-blue-600 font-semibold"
															: "text-stone-500"
												}
											>
												{psi > 1.005
													? "C larger than B"
													: psi < 0.995
														? "B larger than C"
														: "equal areas"}
											</span>
										</div>
										<div className="flex justify-between mt-0.5">
											<span className="text-stone-400">psi = F_C / F_B</span>
											<span className="font-mono font-semibold">
												{psi.toFixed(4)}
											</span>
										</div>
										<div className="flex justify-between">
											<span className="text-stone-400">D_C → D_B</span>
											<span className="font-mono">
												{(dFromF(effectiveFC) * lu).toFixed(1)} →{" "}
												{(dFromF(effectiveFBranch) * lu).toFixed(1)} {lu_label}
											</span>
										</div>
									</div>
								)}

								<InitialGuessEditor node={selectedNode} />

								{/* Post-solve diagnostics */}
								{selectedNode.data.result && (
									<div className="mt-4 p-3 bg-stone-50 border border-stone-200 rounded">
										<h3 className="text-xs font-bold text-stone-500 uppercase mb-3">
											Flow Distribution
										</h3>
										{selectedNode.data.result.correlation_extrapolated > 0 && (
											<div className="text-[10px] text-amber-700 bg-amber-50 border border-amber-200 rounded px-2 py-1 mb-3">
												Inputs outside Bassett (2001) validated range -- results
												are extrapolated.
											</div>
										)}
										<div className="grid grid-cols-2 gap-x-2 gap-y-3">
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													ṁ common
												</span>
												<span className="font-mono text-xs font-bold">
													{(selectedNode.data.result.m_dot_com ?? 0).toFixed(4)}{" "}
													<span className="text-[9px] font-normal text-stone-400">
														kg/s
													</span>
												</span>
											</div>
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													ṁ straight
												</span>
												<span className="font-mono text-xs font-bold">
													{(
														selectedNode.data.result.m_dot_straight ?? 0
													).toFixed(4)}{" "}
													<span className="text-[9px] font-normal text-stone-400">
														kg/s
													</span>
												</span>
											</div>
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													ṁ branch
												</span>
												<span className="font-mono text-xs font-bold">
													{(selectedNode.data.result.m_dot_branch ?? 0).toFixed(
														4,
													)}{" "}
													<span className="text-[9px] font-normal text-stone-400">
														kg/s
													</span>
												</span>
											</div>
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													ṁ straight / ṁ common
												</span>
												<span className="font-mono text-xs font-bold">
													{(
														selectedNode.data.result.mass_flow_ratio ?? 0
													).toFixed(3)}{" "}
													<span className="text-[9px] font-normal text-stone-400">
														kg/kg
													</span>
												</span>
											</div>
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													K straight
												</span>
												<span className="font-mono text-xs font-bold">
													{(selectedNode.data.result.K_straight ?? 0).toFixed(
														3,
													)}
												</span>
											</div>
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													K branch
												</span>
												<span className="font-mono text-xs font-bold">
													{(selectedNode.data.result.K_branch ?? 0).toFixed(3)}
												</span>
											</div>
										</div>
									</div>
								)}
							</div>
						);
					})()}

				{selectedNode.type === "mpce_tee" && (
					<div className="flex flex-col gap-4">
						{/* Flow direction */}
						<div className="flex flex-col gap-2">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Flow Direction
							</label>
							<select
								className="p-2 border rounded bg-white text-xs border-stone-200"
								value={selectedNode.data.flow_direction ?? "branch"}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										flow_direction: e.target.value,
									})
								}
							>
								<option value="branch">Branch (1 inlet, 2 outlets)</option>
								<option value="merge">Merge (2 inlets, 1 outlet)</option>
							</select>
							<p className="text-[9px] text-gray-400 italic">
								Constrained: the Mynard residual asserts at solve time and
								errors if the observed flow direction disagrees with this
								declaration.
							</p>
						</div>

						{/* Branch angle */}
						<div className="flex flex-col gap-2">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Branch Angle (deg)
							</label>
							<NumericInput
								value={selectedNode.data.theta_deg ?? 90}
								onChange={(val) =>
									updateNodeData(selectedNode.id, { theta_deg: val })
								}
								className="p-2 border rounded"
								placeholder="90"
							/>
							<p className="text-[9px] text-gray-400 italic">
								Branch arm angle measured from the main axis. Mynard's smooth K
								relation is valid across the full 0-180 deg range.
							</p>
						</div>
					</div>
				)}

				{selectedNode.type === "vortex" && (
					<div className="flex flex-col gap-4">
						<div className="flex flex-col gap-1">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Rotational Speed
							</label>
							<div className="flex items-center gap-1">
								<NumericInput
									id={`omega_rpm_${selectedNode.id}`}
									className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200 flex-1"
									value={
										selectedNode.data.omega_rpm != null
											? selectedNode.data.omega_rpm *
												(unitPreferences.rotSpeed === "Hz" ? 1 / 60 : 1)
											: null
									}
									onChange={(val) =>
										updateNodeData(selectedNode.id, {
											omega_rpm:
												val * (unitPreferences.rotSpeed === "Hz" ? 60 : 1),
										})
									}
									onClear={() =>
										updateNodeData(selectedNode.id, { omega_rpm: null })
									}
									min={0}
									placeholder={`global: ${(
										(solverSettings.omega_rpm ?? 0) *
											(unitPreferences.rotSpeed === "Hz" ? 1 / 60 : 1)
									).toFixed(unitPreferences.rotSpeed === "Hz" ? 3 : 0)}`}
								/>
								<select
									value={unitPreferences.rotSpeed}
									onChange={(e) =>
										setRotSpeedUnit(e.target.value as "rpm" | "Hz")
									}
									className="w-14 h-[30px] border rounded text-[9px] bg-white outline-none shrink-0"
								>
									<option value="rpm">rpm</option>
									<option value="Hz">Hz</option>
								</select>
							</div>
						</div>
						<LengthInput
							id={`r_c_${selectedNode.id}`}
							label="Core Radius r_c"
							value={selectedNode.data.r_c ?? 0.02}
							onChange={(val) => updateNodeData(selectedNode.id, { r_c: val })}
						/>
						<LengthInput
							id={`r_out_${selectedNode.id}`}
							label="Outer Radius r_out"
							value={selectedNode.data.r_out ?? 0.1}
							onChange={(val) =>
								updateNodeData(selectedNode.id, { r_out: val })
							}
						/>
						<LengthInput
							id={`r_in_${selectedNode.id}`}
							label="Inner Radius r_in"
							value={selectedNode.data.r_in ?? 0.0}
							onChange={(val) => updateNodeData(selectedNode.id, { r_in: val })}
						/>
						<div className="flex flex-col gap-1">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Shape Parameter n
							</label>
							<NumericInput
								id={`n_${selectedNode.id}`}
								className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200"
								value={selectedNode.data.n ?? 2.0}
								onChange={(val) => updateNodeData(selectedNode.id, { n: val })}
								min={1}
								placeholder="2.0"
							/>
							<span className="text-[9px] text-stone-400">
								Vatistas n ≥ 1, default n=2
							</span>
						</div>
						{selectedNode.data.result && (
							<div className="rounded border border-stone-200 bg-stone-50 px-3 py-2 text-xs">
								<div className="font-bold uppercase text-stone-400 mb-1 text-[9px]">
									Result
								</div>
								<div className="flex justify-between">
									<span className="text-stone-500">ΔP vortex</span>
									<span className="font-mono font-bold">
										{((selectedNode.data.result.dP_vortex ?? 0) / 1000).toFixed(
											2,
										)}{" "}
										kPa
									</span>
								</div>
								<div className="flex justify-between">
									<span className="text-stone-500">ṁ</span>
									<span className="font-mono font-bold">
										{(selectedNode.data.result.m_dot ?? 0).toFixed(4)} kg/s
									</span>
								</div>
							</div>
						)}
					</div>
				)}

				{selectedNode.type === "discrete_loss" &&
					(() => {
						// Compute default area from upstream node or channel element
						const upstreamEdge = edges.find(
							(e) => e.target === selectedNode.id && !e.data?.type,
						);
						const upstreamNode = upstreamEdge
							? nodes.find((n) => n.id === upstreamEdge.source)
							: undefined;
						// Walk one more hop upstream to reach a channel when there is an
						// auto-inserted junction node between the channel and this element.
						const upstreamChannelNode =
							upstreamNode?.type === "channel"
								? upstreamNode
								: (() => {
										if (!upstreamNode) return undefined;
										const hopEdge = edges.find(
											(e) => e.target === upstreamNode.id && !e.data?.type,
										);
										const hopNode = hopEdge
											? nodes.find((n) => n.id === hopEdge.source)
											: undefined;
										return hopNode?.type === "channel" ? hopNode : undefined;
									})();
						const inheritedArea: number | undefined =
							upstreamNode?.data?.area && upstreamNode.data.area > 0
								? upstreamNode.data.area
								: upstreamChannelNode?.data?.D != null
									? (Math.PI / 4) * upstreamChannelNode.data.D ** 2
									: undefined;

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
											Inherited from upstream (
											{upstreamChannelNode?.id ?? upstreamNode?.id}). Edit to
											override.
										</p>
									)}
									{!userSetArea && inheritedArea == null && (
										<p className="text-[9px] text-gray-400 italic -mt-1">
											Connect an upstream channel or node to auto-discover area.
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

								<div className="grid grid-cols-2 gap-3 border-t pt-3">
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Nu Multiplier
										</label>
										<NumericInput
											value={selectedNode.data.Nu_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													Nu_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.Nu_multiplier ?? 1}`}
										</span>
									</div>
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											F Multiplier
										</label>
										<NumericInput
											value={selectedNode.data.f_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													f_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.f_multiplier ?? 1}`}
										</span>
									</div>
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

				{selectedNode.type === "combustor" &&
					(() => {
						const upstreamEdge = edges.find(
							(e) => e.target === selectedNode.id && !e.data?.type,
						);
						const upstreamNode = upstreamEdge
							? nodes.find((n) => n.id === upstreamEdge.source)
							: undefined;
						const walkChannelD = (
							nodeId: string,
							depth = 0,
						): number | undefined => {
							if (depth > 20) return undefined;
							const e = edges.find(
								(ee) => ee.target === nodeId && !ee.data?.type,
							);
							const n = e ? nodes.find((nn) => nn.id === e.source) : undefined;
							if (!n) return undefined;
							if (n.type === "channel") {
								const d = n.data.D as number | null | undefined;
								const dh = n.data.Dh as number | null | undefined;
								if (d != null) return dh ?? d;
								return walkChannelD(n.id, depth + 1);
							}
							if (n.type === "area_change") {
								const dhAc = n.data.D_h as number | null | undefined;
								if (dhAc != null && dhAc > 0) return dhAc;
								const f = n.data.F1 as number | null | undefined;
								if (f != null && f > 0) return Math.sqrt((4 * f) / Math.PI);
								return walkChannelD(n.id, depth + 1);
							}
							if (n.type === "tee_junction") {
								const fc = n.data.F_C as number | null | undefined;
								if (fc != null && fc > 0) return Math.sqrt((4 * fc) / Math.PI);
							}
							return (
								(n.data?.Dh as number | undefined) ??
								(n.data?.area != null
									? Math.sqrt((4 * n.data.area) / Math.PI)
									: undefined)
							);
						};
						const inheritedDhComb: number | undefined = walkChannelD(
							selectedNode.id,
						);
						const userSetDhComb =
							selectedNode.data.Dh !== null &&
							selectedNode.data.Dh !== undefined;
						const circularDhComb = Math.sqrt(
							(4 * (selectedNode.data.area ?? 0.1)) / Math.PI,
						);
						return (
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
										<div className="flex items-center justify-between">
											<label className="text-xs font-bold text-gray-500 uppercase">
												Dh
											</label>
											{userSetDhComb && (
												<button
													type="button"
													className="text-[9px] text-blue-500 hover:underline"
													onClick={() =>
														updateNodeData(selectedNode.id, {
															Dh: null,
														})
													}
												>
													Reset to inherited
												</button>
											)}
										</div>
										<LengthInput
											id={`Dh_comb_${selectedNode.id}`}
											label=""
											value={
												userSetDhComb
													? selectedNode.data.Dh
													: (inheritedDhComb ?? circularDhComb)
											}
											placeholder={inheritedDhComb ?? circularDhComb}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													Dh: val,
												})
											}
											onClear={() =>
												updateNodeData(selectedNode.id, {
													Dh: null,
												})
											}
											inputClassName={
												userSetDhComb ? "font-semibold" : "text-gray-400"
											}
										/>
										{!userSetDhComb && inheritedDhComb != null && (
											<p className="text-[9px] text-gray-400 italic -mt-1">
												Inherited from upstream ({upstreamNode?.id}).
											</p>
										)}
										{!userSetDhComb && inheritedDhComb == null && (
											<p className="text-[9px] text-gray-400 italic -mt-1">
												Circular from area (D ={" "}
												{(
													circularDhComb *
													(unitPreferences.length === "mm" ? 1000 : 1)
												).toFixed(2)}{" "}
												{unitPreferences.length}).
											</p>
										)}
									</div>
								</div>
								<div className="grid grid-cols-2 gap-3">
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
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.Nu_multiplier ?? 1}`}
										</span>
									</div>
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											F Multiplier
										</label>
										<NumericInput
											id={`f_multiplier_comb_${selectedNode.id}`}
											value={selectedNode.data.f_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													f_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.f_multiplier ?? 1}`}
										</span>
									</div>
								</div>
								<SurfaceEnhancementInspector
									surface={selectedNode.data.surface || { type: "smooth" }}
									onChange={(surface) =>
										updateNodeData(selectedNode.id, { surface })
									}
								/>

								<InitialGuessEditor node={selectedNode} />
							</div>
						);
					})()}

				{selectedNode.type === "momentum_chamber" &&
					(() => {
						const upstreamEdge = edges.find(
							(e) => e.target === selectedNode.id && !e.data?.type,
						);
						const upstreamNode = upstreamEdge
							? nodes.find((n) => n.id === upstreamEdge.source)
							: undefined;
						const walkChannelD = (
							nodeId: string,
							depth = 0,
						): number | undefined => {
							if (depth > 20) return undefined;
							const e = edges.find(
								(ee) => ee.target === nodeId && !ee.data?.type,
							);
							const n = e ? nodes.find((nn) => nn.id === e.source) : undefined;
							if (!n) return undefined;
							if (n.type === "channel") {
								const d = n.data.D as number | null | undefined;
								const dh = n.data.Dh as number | null | undefined;
								if (d != null) return dh ?? d;
								return walkChannelD(n.id, depth + 1);
							}
							if (n.type === "area_change") {
								const dhAc = n.data.D_h as number | null | undefined;
								if (dhAc != null && dhAc > 0) return dhAc;
								const f = n.data.F1 as number | null | undefined;
								if (f != null && f > 0) return Math.sqrt((4 * f) / Math.PI);
								return walkChannelD(n.id, depth + 1);
							}
							if (n.type === "tee_junction") {
								const fc = n.data.F_C as number | null | undefined;
								if (fc != null && fc > 0) return Math.sqrt((4 * fc) / Math.PI);
							}
							return (
								(n.data?.Dh as number | undefined) ??
								(n.data?.area != null
									? Math.sqrt((4 * n.data.area) / Math.PI)
									: undefined)
							);
						};
						const inheritedDhMom: number | undefined = walkChannelD(
							selectedNode.id,
						);
						const userSetDhMom =
							selectedNode.data.Dh !== null &&
							selectedNode.data.Dh !== undefined;
						const circularDh = Math.sqrt(
							(4 * (selectedNode.data.area ?? 0.1)) / Math.PI,
						);
						return (
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
										<div className="flex items-center justify-between">
											<label className="text-xs font-bold text-gray-500 uppercase">
												Dh
											</label>
											{userSetDhMom && (
												<button
													type="button"
													className="text-[9px] text-blue-500 hover:underline"
													onClick={() =>
														updateNodeData(selectedNode.id, {
															Dh: null,
														})
													}
												>
													Reset to inherited
												</button>
											)}
										</div>
										<LengthInput
											id={`Dh_mom_${selectedNode.id}`}
											label=""
											value={
												userSetDhMom
													? selectedNode.data.Dh
													: (inheritedDhMom ?? circularDh)
											}
											placeholder={inheritedDhMom ?? circularDh}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													Dh: val,
												})
											}
											onClear={() =>
												updateNodeData(selectedNode.id, {
													Dh: null,
												})
											}
											inputClassName={
												userSetDhMom ? "font-semibold" : "text-gray-400"
											}
										/>
										{!userSetDhMom && inheritedDhMom != null && (
											<p className="text-[9px] text-gray-400 italic -mt-1">
												Inherited from upstream ({upstreamNode?.id}).
											</p>
										)}
										{!userSetDhMom && inheritedDhMom == null && (
											<p className="text-[9px] text-gray-400 italic -mt-1">
												Circular from area (D ={" "}
												{(
													circularDh *
													(unitPreferences.length === "mm" ? 1000 : 1)
												).toFixed(2)}{" "}
												{unitPreferences.length}).
											</p>
										)}
									</div>
								</div>
								<div className="grid grid-cols-2 gap-3">
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
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.Nu_multiplier ?? 1}`}
										</span>
									</div>
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											F Multiplier
										</label>
										<NumericInput
											id={`f_multiplier_mom_${selectedNode.id}`}
											value={selectedNode.data.f_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													f_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.f_multiplier ?? 1}`}
										</span>
									</div>
								</div>
								<SurfaceEnhancementInspector
									surface={selectedNode.data.surface || { type: "smooth" }}
									onChange={(surface) =>
										updateNodeData(selectedNode.id, { surface })
									}
								/>

								<InitialGuessEditor node={selectedNode} />
							</div>
						);
					})()}

				{selectedNode.data.result && (
					<div className="mt-4 p-3 bg-stone-50 border border-stone-200 rounded">
						<div className="flex justify-between items-center mb-3">
							<h3 className="text-xs font-bold text-stone-500 uppercase flex items-center gap-2">
								<Activity size={14} />
								Live Telemetry
							</h3>
						</div>

						{displaySettings.length === 0 ? (
							<p className="text-[10px] text-stone-400 italic text-center py-2">
								Select fields in the Display Fields panel to show diagnostics
								here.
							</p>
						) : (
							<div className="flex flex-col gap-3">
								<div className="grid grid-cols-2 gap-x-2 gap-y-3">
									{displaySettings.map((key) => {
										const val = resolveField(selectedNode.data.result, key);
										if (val === undefined) return null;
										const meta = QUANTITY_CATALOGUE[key];
										return (
											<div key={key} className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold">
													{meta?.label ?? key}
												</span>
												<span className="font-mono text-xs font-bold whitespace-nowrap">
													{meta ? meta.format(val) : val.toFixed(4)}{" "}
													<span className="text-[9px] font-normal text-stone-400 ml-0.5">
														{meta?.unit ?? ""}
													</span>
												</span>
											</div>
										);
									})}
								</div>

								{/* Composition — not a scalar field, always shown when present */}
								{selectedNode.data.result.state?.X && (
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
						)}
					</div>
				)}
			</aside>
		);
	}

	// Edge Inspector
	if (selectedEdge && selectedEdge.data?.type === "thermal") {
		const edgeQ: number | undefined = selectedEdge.data.result?.Q;
		const flowsForward = edgeQ === undefined || edgeQ >= 0;

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
			let hotSideShortcut = false;
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
					targetX = flowsForward ? runningX : runningX + L;
					probeLabel = `L${safeIndex + 1} ${flowsForward ? "Hot" : "Cold"} Side`;
				} else if (type === "avg") {
					targetX = runningX + L / 2;
					probeLabel = `L${safeIndex + 1} Average`;
				} else {
					targetX = flowsForward ? runningX + L : runningX;
					probeLabel = `L${safeIndex + 1} ${flowsForward ? "Cold" : "Hot"} Side`;
				}
			} else if (!selectedEdge.data.probe_mode) {
				hotSideShortcut = true;
				probeLabel = "Hot Side";
			} else {
				const scale = unitPreferences.length === "mm" ? 1000 : 1;
				probeLabel = `Custom Depth d=${(targetX * scale).toFixed(unitPreferences.length === "mm" ? 1 : 3)}${unitPreferences.length}`;
			}

			if (tInt.length === layers.length + 1) {
				if (hotSideShortcut) {
					probeTemp = flowsForward ? tInt[0] : tInt[tInt.length - 1];
				} else if (targetX <= 0) {
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
					if (!found) probeTemp = tInt[tInt.length - 1];
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
							updateEdgeData(selectedEdge.id, { label: e.target.value }, true)
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
							updateEdgeData(
								selectedEdge.id,
								{ show_label: e.target.checked },
								true,
							)
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
									updateEdgeData(
										selectedEdge.id,
										{ probe_mode: "custom" },
										true,
									);
								} else {
									const [_, type, idx] = val.split(":");
									updateEdgeData(
										selectedEdge.id,
										{
											probe_mode: "preset",
											probe_preset: { type, index: Number.parseInt(idx, 10) },
										},
										true,
									);
								}
							}}
						>
							<option value="custom">Custom Position (d=x)</option>
							{/* Predefined locations for Layer 1 */}
							<option value="preset:hot:0">
								L1 {flowsForward ? "Hot" : "Cold"} Side
							</option>
							<option value="preset:avg:0">L1 Average</option>
							<option value="preset:cold:0">
								L1 {flowsForward ? "Cold" : "Hot"} Side
							</option>
							{/* Dynamic locations for additional layers */}
							{(selectedEdge.data.layers || [])
								.slice(1)
								.map((_: any, i: number) => (
									<React.Fragment key={i + 1}>
										<option value={`preset:avg:${i + 1}`}>
											L{i + 2} Average
										</option>
										<option value={`preset:cold:${i + 1}`}>
											L{i + 2} {flowsForward ? "Cold" : "Hot"} Side
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
								updateEdgeData(selectedEdge.id, { probe_depth: clamped }, true);
							}}
							placeholder="0.00"
						/>
					)}
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
									{(() => {
										const tInt = selectedEdge.data.result.T_interface as
											| number[]
											| undefined;
										if (tInt && tInt.length >= 2) {
											return flowsForward
												? tInt[0].toFixed(1)
												: tInt[tInt.length - 1].toFixed(1);
										}
										return selectedEdge.data.result.T_hot?.toFixed(1);
									})()} K
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
							updateEdgeData(
								selectedEdge.id,
								{ label: e.target.value || undefined },
								true,
							)
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
							updateEdgeData(
								selectedEdge.id,
								{ show_label: e.target.checked },
								true,
							)
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
		// Boundaries: Pt is solved (mass-flow inlet) or fixed (pressure
		// inlet). Short key "P" is broadcast to both by the backend; we show
		// it only where at least one of the two is actually an unknown.
		mass_boundary: [{ key: "P", unit: "Pa" }],
		// Elements (edges in React Flow but still nodes here): m_dot only.
		channel: [{ key: "m_dot", unit: "kg/s" }],
		orifice: [{ key: "m_dot", unit: "kg/s" }],
		discrete_loss: [{ key: "m_dot", unit: "kg/s" }],
		area_change: [{ key: "m_dot", unit: "kg/s" }],
		// Tee junction: two independent mass-flow unknowns.
		tee_junction: [
			{ key: "m_dot_com", unit: "kg/s" },
			{ key: "m_dot_branch", unit: "kg/s" },
		],
		vortex: [{ key: "m_dot", unit: "kg/s" }],
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
