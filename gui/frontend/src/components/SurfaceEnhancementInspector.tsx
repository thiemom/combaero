import type React from "react";
import LengthInput from "./LengthInput";
import NumericInput from "./NumericInput";

interface SurfaceModelData {
	type:
		| "smooth"
		| "ribbed"
		| "dimpled"
		| "pin_fin"
		| "impingement"
		| "single_jet_impingement";
	// Ribbed, rebuilt in 0.8.0 on a provenanced correlation set.
	e_D?: number;
	p_e?: number;
	alpha_deg?: number;
	W_H?: number;
	n_ribbed_walls?: number;
	smooth_wall_Nu_multiplier?: number;
	// Kept for pre-0.7.0 networks whose rib nodes used the old field name.
	pitch_to_height?: number;
	d_Dh?: number;
	h_d?: number;
	S_d?: number;
	pin_diameter?: number;
	S_D?: number;
	X_D?: number;
	N_rows?: number;
	is_staggered?: boolean;
	channel_height?: number;
	// Impingement (array), rebuilt in 0.9.0 on Florschuetz, Truman and Metzger
	// (1981) with a real Gc/Gj crossflow term -- see issue #337. d_jet is
	// shared with the single-jet block below.
	d_jet?: number;
	xn_d?: number;
	yn_d?: number;
	z_d?: number;
	row?: number;
	C_D?: number;
	// Impingement, single jet (Goldstein, Behbahani and Heppelmann 1986).
	bc?: "constant_heat_flux" | "constant_wall_temperature";
	L_D?: number;
	R_D?: number;
	// Kept for pre-0.7.0 networks whose impingement nodes used the old,
	// crossflow-less field names -- see #332 for why that correlation was
	// removed. Unreachable from this dropdown's new "impingement" defaults.
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
				e_D: 0.06,
				p_e: 10.0,
				alpha_deg: 90.0,
				W_H: 1.0,
				n_ribbed_walls: 2,
				smooth_wall_Nu_multiplier: 1.0,
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
				xn_d: 8.0,
				yn_d: 6.0,
				z_d: 2.0,
				row: 1,
				C_D: 0.79,
			},
			single_jet_impingement: {
				type: "single_jet_impingement",
				bc: "constant_heat_flux",
				d_jet: 0.003,
				L_D: 7.75,
				R_D: 5.0,
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
					<option value="impingement">Impingement (Jet Array)</option>
					<option value="single_jet_impingement">
						Impingement (Single Jet)
					</option>
					{/* Dimpled and pin-fin were removed in 0.7.0: their correlations
					    could not be traced to their cited sources, and remain
					    deferred (issue #339). Ribbed returned in 0.8.0 on a
					    provenanced correlation set; impingement returned in 0.9.0
					    on Florschuetz (1981) with a real crossflow term (#337).

					    Dimpled/pin-fin's parameter blocks below are deliberately
					    KEPT. They are unreachable from this dropdown, so no new one
					    can be created -- but a network saved before 0.7.0 still
					    carries one, and removing the blocks would hide its
					    parameters and risk losing them on the next save. The user
					    can see what is there; the backend explains on solve why it
					    will not run. */}
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
							value={surface.p_e || 10.0}
							onChange={(val) => updateFields({ p_e: val })}
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
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							W/H (Channel Aspect Ratio)
						</label>
						<NumericInput
							value={surface.W_H || 1.0}
							onChange={(val) => updateFields({ W_H: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							Ribbed walls (1, 2 or 4)
						</label>
						<select
							value={surface.n_ribbed_walls ?? 2}
							onChange={(e) =>
								updateFields({ n_ribbed_walls: Number(e.target.value) })
							}
							className="p-1 border rounded text-xs"
						>
							<option value={1}>1</option>
							<option value={2}>2 (opposite walls)</option>
							<option value={4}>4 (all walls)</option>
						</select>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							Smooth-wall Nu multiplier
						</label>
						<NumericInput
							value={surface.smooth_wall_Nu_multiplier || 1.0}
							onChange={(val) =>
								updateFields({ smooth_wall_Nu_multiplier: val })
							}
							className="p-1 border rounded text-xs"
						/>
						<span className="text-[9px] text-stone-400">
							Plain smooth walls under-predict the channel average by about 20%
							against Han&apos;s data, because ribs enhance the adjacent smooth
							wall too. ~1.67 reproduces Han.
						</span>
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
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							xn/d (Streamwise Spacing Ratio)
						</label>
						<NumericInput
							value={surface.xn_d || 8.0}
							onChange={(val) => updateFields({ xn_d: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							yn/d (Spanwise Spacing Ratio)
						</label>
						<NumericInput
							value={surface.yn_d || 6.0}
							onChange={(val) => updateFields({ yn_d: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							z/d (Channel Height Ratio)
						</label>
						<NumericInput
							value={surface.z_d || 2.0}
							onChange={(val) => updateFields({ z_d: val })}
							className="p-1 border rounded text-xs"
						/>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							Row (1 = first, no crossflow)
						</label>
						<NumericInput
							value={surface.row ?? 1}
							onChange={(val) =>
								updateFields({ row: Math.max(1, Math.round(val)) })
							}
							className="p-1 border rounded text-xs"
						/>
						<span className="text-[9px] text-stone-400">
							One row per element. A full array is a chain of these, each with
							its own row number -- downstream rows see progressively more
							crossflow from the rows upstream of them.
						</span>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							C_D (Jet Plate Discharge Coefficient)
						</label>
						<NumericInput
							value={surface.C_D ?? 0.79}
							onChange={(val) => updateFields({ C_D: val })}
							className="p-1 border rounded text-xs"
						/>
						<span className="text-[9px] text-stone-400">
							0.79 is Florschuetz&apos;s own recommended default absent a
							measured value; combaero has no discharge-coefficient correlation
							of its own for a jet-plate array yet (issue #375).
						</span>
					</div>
				</div>
			)}

			{currentType === "single_jet_impingement" && (
				<div className="grid grid-cols-1 gap-3 bg-stone-50 p-2 rounded border border-stone-100">
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							Boundary Condition
						</label>
						<select
							value={surface.bc || "constant_heat_flux"}
							onChange={(e) =>
								updateFields({
									bc: e.target.value as SurfaceModelData["bc"],
								})
							}
							className="p-1 border rounded text-xs"
						>
							<option value="constant_heat_flux">Constant Heat Flux</option>
							<option value="constant_wall_temperature">
								Constant Wall Temperature
							</option>
						</select>
					</div>
					<LengthInput
						label="Jet Diameter"
						value={surface.d_jet || 0.003}
						onChange={(val) => updateFields({ d_jet: val })}
					/>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							L/D (Jet-to-Plate Spacing Ratio)
						</label>
						<NumericInput
							value={surface.L_D ?? 7.75}
							onChange={(val) => updateFields({ L_D: val })}
							className="p-1 border rounded text-xs"
						/>
						<span className="text-[9px] text-stone-400">
							7.75 is the correlation&apos;s own optimum spacing.
						</span>
					</div>
					<div className="flex flex-col gap-1">
						<label className="text-[10px] text-stone-500">
							R/D (Radial Position Ratio)
						</label>
						<NumericInput
							value={surface.R_D ?? 5.0}
							onChange={(val) => updateFields({ R_D: val })}
							className="p-1 border rounded text-xs"
						/>
						<span className="text-[9px] text-stone-400">
							Nu is local to this radial position, not area-averaged -- pick a
							representative value for the target patch this surface covers.
						</span>
					</div>
				</div>
			)}
		</div>
	);
};

export default SurfaceEnhancementInspector;
