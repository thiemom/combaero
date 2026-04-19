import { useEffect, useState } from "react";
import useStore from "../store/useStore";
import NumericInput from "./NumericInput";
import UnitInput from "./UnitInput";

interface CompositionEditorProps {
	nodeId: string;
	data: any;
}

const DEFAULT_COMPOSITION = {
	source: "dry_air",
	mode: "mole",
	custom_fractions: {},
	relative_humidity: 0.6,
	ambient_T: 288.15,
	ambient_P: 101325.0,
};

const CompositionEditor: React.FC<CompositionEditorProps> = ({
	nodeId,
	data,
}) => {
	const { speciesMetadata, fetchSpeciesMetadata, updateNodeData } = useStore();
	const [localComp, setLocalComp] = useState(() => ({
		...DEFAULT_COMPOSITION,
		...(data.composition || {}),
	}));

	useEffect(() => {
		if (!speciesMetadata) {
			fetchSpeciesMetadata();
		}
	}, [speciesMetadata, fetchSpeciesMetadata]);

	const updateComp = (updates: any) => {
		const newComp = { ...localComp, ...updates };
		setLocalComp(newComp);
		updateNodeData(nodeId, { composition: newComp });
	};

	const handleFractionChange = (species: string, val: number) => {
		const custom_fractions = { ...localComp.custom_fractions, [species]: val };
		updateComp({ custom_fractions });
	};

	return (
		<div className="flex flex-col gap-3 p-3 bg-stone-50 rounded-lg border border-stone-200">
			<div className="text-xs font-bold text-stone-500 uppercase">
				Composition
			</div>

			<div className="grid grid-cols-2 gap-2">
				<div className="flex flex-col gap-1">
					<label className="text-[10px] text-stone-400 uppercase font-bold">
						Source
					</label>
					<select
						className="p-1 text-sm border rounded bg-white"
						value={localComp.source}
						onChange={(e) => {
							const newSource = e.target.value;
							const updates: any = { source: newSource };

							if (newSource === "custom" && speciesMetadata?.names?.length) {
								// Smart initialization: copy from current resolved state if available
								let initFractions: Record<string, number> = {};

								if (data.result?.state) {
									const currentFractions =
										localComp.mode === "mass"
											? data.result.state.Y
											: data.result.state.X;
									if (currentFractions && currentFractions.length > 0) {
										currentFractions.forEach((val: number, i: number) => {
											if (val > 1e-6)
												initFractions[speciesMetadata.names[i]] = val;
										});
									}
								}

								// Fallback to first species equal to 1.0 if no resolved state or it failed
								if (Object.keys(initFractions).length === 0) {
									initFractions = { [speciesMetadata.names[0]]: 1.0 };
								}

								// Actually set the initial fractions if they are missing or empty
								// OR if we are switching from a non-custom source, overwrite them
								if (
									!localComp.custom_fractions ||
									Object.keys(localComp.custom_fractions).length === 0 ||
									localComp.source !== "custom"
								) {
									updates.custom_fractions = initFractions;
								}
							}
							updateComp(updates);
						}}
					>
						<option value="dry_air">Dry Air</option>
						<option value="humid_air">Humid Air</option>
						<option value="fuel">Fuel (CH4)</option>
						<option value="custom">Custom</option>
					</select>
				</div>

				<div className="flex flex-col gap-1">
					<label className="text-[10px] text-stone-400 uppercase font-bold">
						Mode
					</label>
					<select
						className="p-1 text-sm border rounded bg-white"
						value={localComp.mode}
						onChange={(e) => updateComp({ mode: e.target.value })}
					>
						<option value="mole">Mole Fractions</option>
						<option value="mass">Mass Fractions</option>
					</select>
				</div>
			</div>

			{localComp.source === "humid_air" && (
				<div className="flex flex-col gap-2 border-t pt-2 border-stone-200">
					<div className="grid grid-cols-2 gap-2">
						<div className="flex flex-col gap-1">
							<label className="text-[10px] text-stone-400 uppercase font-bold">
								Ambient T (K)
							</label>
							<NumericInput
								value={localComp.ambient_T}
								onChange={(val) => updateComp({ ambient_T: val })}
								className="p-1 text-sm border rounded bg-white"
							/>
						</div>
						<div className="flex flex-col gap-1">
							<label className="text-[10px] text-stone-400 uppercase font-bold">
								Rel. Humidity (0-1)
							</label>
							<NumericInput
								value={localComp.relative_humidity}
								onChange={(val) => updateComp({ relative_humidity: val })}
								className="p-1 text-sm border rounded bg-white"
							/>
						</div>
					</div>
					<UnitInput
						label="Ambient Pressure"
						value={localComp.ambient_P}
						onChange={(val) => updateComp({ ambient_P: val })}
					/>
				</div>
			)}

			{localComp.source === "custom" && speciesMetadata && (
				<div className="flex flex-col gap-2 mt-2 border-t pt-2 border-stone-200">
					<div className="flex items-center justify-between">
						<div className="text-[10px] text-stone-400 uppercase font-bold">
							Species Fractions
						</div>
						<div
							className={`text-[10px] font-mono px-1 rounded ${
								Math.abs(
									(
										Object.values(localComp.custom_fractions || {}) as number[]
									).reduce((a: number, b: number) => a + (b || 0), 0) - 1.0,
								) > 1e-4
									? "bg-amber-100 text-amber-700"
									: "bg-emerald-100 text-emerald-700"
							}`}
						>
							Sum:{" "}
							{(Object.values(localComp.custom_fractions || {}) as number[])
								.reduce((a: number, b: number) => a + (b || 0), 0)
								.toFixed(4)}
						</div>
					</div>

					<div className="max-h-48 overflow-y-auto pr-1 flex flex-col gap-1">
						{speciesMetadata.names.map((name) => (
							<div
								key={name}
								className="flex items-center justify-between gap-2"
							>
								<span className="text-xs font-mono">{name}</span>
								<NumericInput
									placeholder="0.00"
									className="w-20 p-1 text-xs border rounded bg-white text-right"
									value={localComp.custom_fractions?.[name] || 0}
									onChange={(val) => handleFractionChange(name, val)}
								/>
							</div>
						))}
					</div>

					<button
						type="button"
						onClick={() => {
							const values = Object.values(
								localComp.custom_fractions || {},
							) as number[];
							const sum = values.reduce(
								(a: number, b: number) => a + (b || 0),
								0,
							);
							if (sum > 0) {
								const normalized: any = {};
								Object.entries(localComp.custom_fractions || {}).forEach(
									([k, v]) => {
										normalized[k] = (v as number) / sum;
									},
								);
								updateComp({ custom_fractions: normalized });
							}
						}}
						className="mt-1 p-1 text-[10px] bg-stone-200 hover:bg-stone-300 text-stone-700 rounded uppercase font-bold"
					>
						Normalize Fractions
					</button>
				</div>
			)}
		</div>
	);
};

export default CompositionEditor;
