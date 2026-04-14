import type React from "react";
import AreaInput from "./AreaInput";
import LengthInput from "./LengthInput";
import NumericInput from "./NumericInput";

interface SurfaceModelData {
	type: "smooth" | "ribbed" | "dimpled" | "pin_fin" | "impingement";
	e_D?: number;
	pitch_to_height?: number;
	alpha_deg?: number;
	d_Dh?: number;
	h_d?: number;
	S_d?: number;
	pin_diameter?: number;
	S_D?: number;
	X_D?: number;
	N_rows?: number;
	is_staggered?: boolean;
	channel_height?: number;
	d_jet?: number;
	z_D?: number;
	x_D?: number;
	y_D?: number;
	A_target?: number;
	Cd_jet?: number;
}

interface Props {
	surface: SurfaceModelData;
	onChange: (data: SurfaceModelData) => void;
}

const SurfaceEnhancementInspector: React.FC<Props> = ({
	surface,
	onChange,
}) => {
	const currentType = surface.type || "smooth";

	const handleTypeChange = (type: SurfaceModelData["type"]) => {
		const defaults: Record<string, any> = {
			smooth: { type: "smooth" },
			ribbed: {
				type: "ribbed",
				e_D: 0.05,
				pitch_to_height: 10.0,
				alpha_deg: 90.0,
			},
			dimpled: { type: "dimpled", d_Dh: 0.2, h_d: 0.15, S_d: 2.0 },
			pin_fin: {
				type: "pin_fin",
				pin_diameter: 0.005,
				channel_height: 0.01,
				S_D: 2.5,
				X_D: 2.5,
				N_rows: 10,
				is_staggered: true,
			},
			impingement: {
				type: "impingement",
				d_jet: 0.002,
				z_D: 4.0,
				x_D: 6.0,
				y_D: 6.0,
				A_target: 0.01,
				Cd_jet: 0.8,
			},
		};
		onChange(defaults[type]);
	};

	const updateFields = (fields: Partial<SurfaceModelData>) => {
		onChange({ ...surface, ...fields });
	};

	return (
		<div className="flex flex-col gap-3 border-t pt-3 mt-2">
			<div className="text-[10px] font-bold text-stone-400 uppercase">
				Surface Geometry
			</div>

			<div className="flex flex-col gap-1">
				<label className="text-[10px] text-stone-500 font-bold uppercase">
					Enhancement Type
				</label>
				<select
					className="p-1.5 border rounded text-xs bg-white w-full"
					value={currentType}
					onChange={(e) =>
						handleTypeChange(e.target.value as SurfaceModelData["type"])
					}
				>
					<option value="smooth">Smooth (Default)</option>
					<option value="ribbed">Ribbed</option>
					<option value="dimpled">Dimpled</option>
					<option value="pin_fin">Pin Fin Array</option>
					<option value="impingement">Jet Impingement</option>
				</select>
			</div>

			{currentType === "ribbed" && (
				<div className="grid grid-cols-1 gap-3 bg-stone-50 p-2 rounded border border-stone-100">
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							e/Dh (Rib Height Ratio)
						</label>
						<NumericInput
							value={surface.e_D || 0.05}
							onChange={(val) => updateFields({ e_D: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							P/e (Rib Pitch Ratio)
						</label>
						<NumericInput
							value={surface.pitch_to_height || 10.0}
							onChange={(val) => updateFields({ pitch_to_height: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">Angle (deg)</label>
						<NumericInput
							value={surface.alpha_deg || 90.0}
							onChange={(val) => updateFields({ alpha_deg: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
				</div>
			)}

			{currentType === "dimpled" && (
				<div className="grid grid-cols-1 gap-3 bg-stone-50 p-2 rounded border border-stone-100">
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							d/Dh (Dimple Dia Ratio)
						</label>
						<NumericInput
							value={surface.d_Dh || 0.2}
							onChange={(val) => updateFields({ d_Dh: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							h/d (Depth Ratio)
						</label>
						<NumericInput
							value={surface.h_d || 0.15}
							onChange={(val) => updateFields({ h_d: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							S/d (Spacing Ratio)
						</label>
						<NumericInput
							value={surface.S_d || 2.0}
							onChange={(val) => updateFields({ S_d: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
				</div>
			)}

			{currentType === "pin_fin" && (
				<div className="grid grid-cols-1 gap-3 bg-stone-50 p-2 rounded border border-stone-100">
					<LengthInput
						label="Pin Diameter"
						value={surface.pin_diameter || 0.005}
						onChange={(val) => updateFields({ pin_diameter: val })}
					/>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							S/D (Transverse Spacing)
						</label>
						<NumericInput
							value={surface.S_D || 2.5}
							onChange={(val) => updateFields({ S_D: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							X/D (Streamwise Spacing)
						</label>
						<NumericInput
							value={surface.X_D || 2.5}
							onChange={(val) => updateFields({ X_D: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">Number of Rows</label>
						<NumericInput
							value={surface.N_rows || 10}
							onChange={(val) =>
								updateFields({ N_rows: Math.max(1, Math.round(val)) })
							}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<LengthInput
						label="Channel Height"
						value={surface.channel_height || 0.01}
						onChange={(val) => updateFields({ channel_height: val })}
					/>
					<div className="flex items-center justify-between gap-2 mt-1 px-1">
						<label className="text-[10px] font-bold text-gray-500 uppercase">
							Staggered
						</label>
						<input
							type="checkbox"
							checked={surface.is_staggered !== false}
							onChange={(e) => updateFields({ is_staggered: e.target.checked })}
							className="w-3 h-3 accent-blue-500"
						/>
					</div>
				</div>
			)}

			{currentType === "impingement" && (
				<div className="grid grid-cols-1 gap-3 bg-stone-50 p-2 rounded border border-stone-100">
					<LengthInput
						label="Jet Diameter"
						value={surface.d_jet || 0.002}
						onChange={(val) => updateFields({ d_jet: val })}
					/>
					<AreaInput
						label="Target Area"
						value={surface.A_target || 0.01}
						onChange={(val) => updateFields({ A_target: val })}
					/>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							z/D (Jet Dist Ratio)
						</label>
						<NumericInput
							value={surface.z_D || 4.0}
							onChange={(val) => updateFields({ z_D: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							x/D (Streamwise Spacing)
						</label>
						<NumericInput
							value={surface.x_D || 6.0}
							onChange={(val) => updateFields({ x_D: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							y/D (Spanwise Spacing)
						</label>
						<NumericInput
							value={surface.y_D || 6.0}
							onChange={(val) => updateFields({ y_D: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">Cd (Discharge)</label>
						<NumericInput
							value={surface.Cd_jet || 0.8}
							onChange={(val) => updateFields({ Cd_jet: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
				</div>
			)}
		</div>
	);
};

export default SurfaceEnhancementInspector;
