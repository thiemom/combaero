import { ArrowDown, ArrowUp, Plus, Trash2 } from "lucide-react";
import type React from "react";
import { useEffect, useState } from "react";
import { fetchMaterials } from "../api";
import LengthInput from "./LengthInput";
import NumericInput from "./NumericInput";

interface WallLayer {
	thickness: number;
	conductivity: number;
	material: string;
}

interface WallLayersEditorProps {
	layers: WallLayer[];
	onChange: (layers: WallLayer[]) => void;
}

const WallLayersEditor: React.FC<WallLayersEditorProps> = ({
	layers,
	onChange,
}) => {
	const [materialList, setMaterialList] = useState<string[]>([]);

	useEffect(() => {
		fetchMaterials().then(setMaterialList).catch(console.error);
	}, []);

	const addLayer = () => {
		const newLayer: WallLayer = {
			thickness: 0.001,
			conductivity: 20.0,
			material: "custom",
		};
		onChange([...layers, newLayer]);
	};

	const removeLayer = (index: number) => {
		const newLayers = layers.filter((_, i) => i !== index);
		onChange(newLayers);
	};

	const updateLayer = (index: number, updates: Partial<WallLayer>) => {
		const newLayers = layers.map((layer, i) =>
			i === index ? { ...layer, ...updates } : layer,
		);
		onChange(newLayers);
	};

	const moveLayer = (index: number, direction: "up" | "down") => {
		const newLayers = [...layers];
		const targetIndex = direction === "up" ? index - 1 : index + 1;
		if (targetIndex >= 0 && targetIndex < newLayers.length) {
			[newLayers[index], newLayers[targetIndex]] = [
				newLayers[targetIndex],
				newLayers[index],
			];
			onChange(newLayers);
		}
	};

	return (
		<div className="flex flex-col gap-3">
			<div className="flex justify-between items-center">
				<label className="text-xs font-bold text-stone-500 uppercase">
					Wall Layers
				</label>
				<button
					type="button"
					onClick={addLayer}
					className="p-1 rounded hover:bg-stone-100 text-stone-600 flex items-center gap-1 text-[10px] font-bold uppercase transition-colors"
					title="Add layer"
				>
					<Plus size={14} /> Add Layer
				</button>
			</div>

			<div className="flex flex-col gap-4">
				{layers.map((layer, index) => {
					const isCustom =
						layer.material.toLowerCase() === "generic" ||
						layer.material.toLowerCase() === "custom";

					return (
						<div
							key={index}
							className="p-3 border rounded-lg bg-stone-50/50 flex flex-col gap-3 relative group shadow-sm"
						>
							<div className="flex justify-between items-center pb-2 border-b border-stone-100">
								<div className="flex items-center gap-2">
									<span className="text-[9px] font-black text-stone-400 bg-white border px-1.5 py-0.5 rounded shadow-sm">
										L{index + 1}
									</span>
									<select
										value={
											layer.material.toLowerCase() === "generic"
												? "custom"
												: layer.material
										}
										onChange={(e) =>
											updateLayer(index, { material: e.target.value })
										}
										className="text-[10px] bg-white border rounded px-1.5 py-0.5 outline-none font-bold text-stone-600 focus:ring-1 focus:ring-orange-500 min-w-[100px]"
									>
										<option value="custom">Custom Material</option>
										{materialList.map((m) => (
											<option key={m} value={m}>
												{m.toUpperCase()}
											</option>
										))}
									</select>
								</div>
								<div className="flex gap-1 opacity-0 group-hover:opacity-100 transition-opacity">
									<button
										type="button"
										onClick={() => moveLayer(index, "up")}
										disabled={index === 0}
										className="p-1 h-6 w-6 flex items-center justify-center rounded hover:bg-white disabled:opacity-30"
									>
										<ArrowUp size={12} />
									</button>
									<button
										type="button"
										onClick={() => moveLayer(index, "down")}
										disabled={index === layers.length - 1}
										className="p-1 h-6 w-6 flex items-center justify-center rounded hover:bg-white disabled:opacity-30"
									>
										<ArrowDown size={12} />
									</button>
									<button
										type="button"
										onClick={() => removeLayer(index)}
										className="p-1 h-6 w-6 flex items-center justify-center rounded hover:bg-red-50 text-red-500"
									>
										<Trash2 size={12} />
									</button>
								</div>
							</div>

							<div className="grid grid-cols-2 gap-3 items-end">
								<LengthInput
									label="Thickness"
									value={layer.thickness}
									onChange={(val) => updateLayer(index, { thickness: val })}
								/>
								<div className="flex flex-col gap-1">
									<label className="text-[10px] font-medium text-stone-500 flex justify-between">
										<span>Conductivity</span>
										{!isCustom && (
											<span className="text-[9px] text-orange-600 font-bold uppercase">
												DB [k(T)]
											</span>
										)}
									</label>
									<div className="relative">
										<NumericInput
											value={layer.conductivity}
											onChange={(val) =>
												isCustom && updateLayer(index, { conductivity: val })
											}
											className={`p-1.5 border rounded text-sm w-full outline-none transition-all ${
												isCustom
													? "bg-white focus:ring-1 focus:ring-orange-500"
													: "bg-stone-100 text-stone-400 cursor-not-allowed border-transparent"
											}`}
											disabled={!isCustom}
										/>
										{!isCustom && (
											<div
												className="absolute inset-0 cursor-not-allowed"
												title="Conductivity is determined by the material database k(T)"
											/>
										)}
									</div>
								</div>
							</div>
						</div>
					);
				})}

				{layers.length === 0 && (
					<div className="flex flex-col items-center justify-center py-8 border-2 border-dashed border-stone-200 rounded-lg text-stone-400 gap-2">
						<span className="text-xs italic">No layers defined</span>
						<button
							type="button"
							onClick={addLayer}
							className="text-[10px] font-bold text-orange-600 hover:text-orange-700 underline"
						>
							Initialize with a default layer
						</button>
					</div>
				)}
			</div>
		</div>
	);
};

export default WallLayersEditor;
