import resultKeysJson from "./resultKeys.generated.json";

/**
 * Quantity catalogue: mapping result keys to labels, units, and formatting logic.
 */
export const QUANTITY_CATALOGUE: Record<
	string,
	{
		label: string;
		unit: string;
		format: (v: number) => string;
		/** Shorter label for the narrow node cards; Telemetry keeps `label`. */
		cardLabel?: string;
	}
> = {
	// Scalar results (accessed via result.*)
	m_dot: { label: "ṁ", unit: "kg/s", format: (v) => v.toFixed(4) },
	mach: { label: "Ma", unit: "—", format: (v) => v.toFixed(4) },
	mach_in: { label: "Ma in", unit: "—", format: (v) => v.toFixed(4) },
	mach_out: { label: "Ma out", unit: "—", format: (v) => v.toFixed(4) },
	mach_loss_ref: {
		label: "Ma (loss ref.)",
		unit: "—",
		format: (v) => v.toFixed(4),
	},
	mach_throat: { label: "Ma throat", unit: "—", format: (v) => v.toFixed(4) },
	zeta: { label: "ζ", unit: "—", format: (v) => v.toFixed(4) },
	sigma_norm: {
		label: "σ/(R ṁ) (2nd law)",
		unit: "—",
		format: (v) => v.toExponential(2),
	},
	phi: { label: "Phi", unit: "—", format: (v) => v.toFixed(3) },
	theta: { label: "Theta", unit: "—", format: (v) => v.toFixed(3) },
	pr_total: { label: "P ratio", unit: "—", format: (v) => v.toFixed(4) },
	pr_static: {
		label: "P ratio (static)",
		unit: "—",
		format: (v) => v.toFixed(4),
	},
	// Area change only: F1/F0 (#773). The effusion plate's and the junction
	// ports' ratios have their own keys below.
	area_ratio: { label: "Area ratio", unit: "—", format: (v) => v.toFixed(3) },
	approach_to_hole_area_ratio: {
		label: "Approach / hole area",
		unit: "—",
		format: (v) => v.toFixed(3),
	},
	dPt: {
		label: "dP (total)",
		cardLabel: "ΔPt",
		unit: "bar",
		format: (v) => (v / 1e5).toFixed(4),
	},
	dP: {
		label: "dP (static)",
		cardLabel: "ΔP",
		unit: "bar",
		format: (v) => (v / 1e5).toFixed(4),
	},
	Re: {
		label: "Re",
		unit: "—",
		format: (v) => v.toLocaleString("en", { maximumFractionDigits: 0 }),
	},
	Cd: { label: "Cd", unit: "—", format: (v) => v.toFixed(4) },
	Nu: { label: "Nu", unit: "—", format: (v) => v.toFixed(2) },
	htc: { label: "htc", unit: "W/m²/K", format: (v) => v.toFixed(1) },
	// Impingement array (#465): the array's own summary of its rows.
	Re_j_min: { label: "Re_j min", unit: "—", format: (v) => v.toFixed(0) },
	Re_j_max: { label: "Re_j max", unit: "—", format: (v) => v.toFixed(0) },
	Gc_Gj_max: { label: "Gc/Gj max", unit: "—", format: (v) => v.toFixed(3) },
	Nu_mean: { label: "Nu mean", unit: "—", format: (v) => v.toFixed(1) },
	htc_mean: { label: "h mean", unit: "W/m²/K", format: (v) => v.toFixed(0) },
	jet_flow_nonuniformity: {
		label: "Jet flow max/min",
		unit: "—",
		format: (v) => v.toFixed(3),
	},
	// Impingement array with an initial (bypass) crossflow (#467)
	m_dot_bypass: {
		label: "Bypass ṁ",
		unit: "kg/s",
		format: (v) => v.toFixed(4),
	},
	bypass_to_jet_flow: {
		label: "Bypass / jets",
		unit: "—",
		format: (v) => v.toFixed(3),
	},
	bypass_Tt_minus_supply: {
		label: "Bypass Tt - supply",
		unit: "K",
		format: (v) => v.toFixed(1),
	},
	heat_transfer_outside_source: {
		label: "eta_r outside data",
		unit: "—",
		format: (v) => (v > 0 ? "yes" : "no"),
	},
	bypass_reversed: {
		label: "Bypass reversed",
		unit: "—",
		format: (v) => (v > 0 ? "yes" : "no"),
	},
	n_holes: { label: "Holes", unit: "—", format: (v) => v.toFixed(0) },
	Cd_in_range: {
		label: "Cd in range",
		unit: "—",
		format: (v) => (v > 0 ? "yes" : "no"),
	},
	// Rayleigh omission of a wall-heated compressible channel (#551)
	rayleigh_omitted_frac: {
		label: "Heating dPt omitted",
		unit: "% Pt",
		format: (v) => (v * 100).toFixed(2),
	},
	rayleigh_omitted: {
		label: "Heating not in momentum",
		unit: "—",
		format: (v) => (v > 0 ? "yes" : "no"),
	},
	heat_area_ratio: {
		label: "Heat area / wetted area",
		unit: "—",
		format: (v) => v.toFixed(2),
	},
	heat_area_exceeds_wetted: {
		label: "Heat area beyond wetted",
		unit: "—",
		format: (v) => (v > 0 ? "yes" : "no"),
	},
	rayleigh_in_momentum: {
		label: "Heating in momentum",
		unit: "—",
		format: (v) => (v > 0 ? "yes (Fanno-Rayleigh)" : "no"),
	},
	// Orifice velocity of approach and ISO 5167 taps (#649/#645/#651)
	beta: { label: "β (declared)", unit: "—", format: (v) => v.toFixed(4) },
	beta_flow: { label: "β (flow)", unit: "—", format: (v) => v.toFixed(4) },
	E: { label: "E (approach)", unit: "—", format: (v) => v.toFixed(4) },
	Cd_pipe_ratio: {
		label: "Cd factor in pipe",
		unit: "—",
		format: (v) => v.toFixed(4),
	},
	P_approach: {
		label: "p1 (approach pipe)",
		unit: "bar",
		format: (v) => (v / 1e5).toFixed(4),
	},
	dP_tap: {
		label: "dP at taps",
		unit: "bar",
		format: (v) => (v / 1e5).toFixed(4),
	},
	dP_permanent: {
		label: "dP permanent",
		unit: "bar",
		format: (v) => (v / 1e5).toFixed(4),
	},
	loss_ratio: {
		label: "Permanent / tap dP",
		unit: "—",
		format: (v) => v.toFixed(3),
	},
	pressure_recovery: {
		label: "Pressure recovery",
		unit: "—",
		format: (v) => (v > 0 ? "yes" : "no"),
	},
	recovery_in_range: {
		label: "Recovery in ISO range",
		unit: "—",
		format: (v) => (v > 0 ? "yes" : "no"),
	},
	// Effusion plate wall (#471)
	T_wall_hot: { label: "T wall (hot)", unit: "K", format: (v) => v.toFixed(0) },
	T_wall_cold: {
		label: "T wall (cold)",
		unit: "K",
		format: (v) => v.toFixed(0),
	},
	q_wall: {
		label: "q wall",
		unit: "kW/m²",
		format: (v) => (v / 1000).toFixed(1),
	},
	eta_overall: { label: "η overall", unit: "—", format: (v) => v.toFixed(3) },
	eta_film: { label: "η film", unit: "—", format: (v) => v.toFixed(3) },
	blowing_ratio: { label: "M", unit: "—", format: (v) => v.toFixed(2) },
	// Effusion liner (#471)
	bleed_fraction: {
		label: "Bleed",
		unit: "%",
		format: (v) => (100 * v).toFixed(1),
	},
	U1_over_Vi_max: {
		label: "U1/Vi max",
		unit: "—",
		format: (v) => v.toFixed(3),
	},
	bleed_nonuniformity: {
		label: "Bleed max/min",
		unit: "—",
		format: (v) => v.toFixed(3),
	},
	T_wall_hot_max: {
		label: "T wall max",
		unit: "K",
		format: (v) => v.toFixed(0),
	},
	Re_surface: {
		label: "Re (surface)",
		unit: "—",
		format: (v) => v.toLocaleString("en", { maximumFractionDigits: 0 }),
	},
	surface_extrapolated: {
		label: "Surface out of range",
		unit: "—",
		format: (v) => (v ? "yes" : "no"),
	},
	Q: { label: "Q", unit: "W", format: (v) => v.toFixed(0) },
	// State results (accessed via result.state.*). P and T are static, Pt and
	// Tt total (#773, #782): until 0.7 T held the total temperature. T is Tt
	// at a plenum and an incompressible node; a moving compressible node's T
	// is on the isentrope through (Tt, Pt). T_static is the old name of a
	// momentum chamber's T, emitted until 0.9.0.
	P: { label: "P", unit: "bar", format: (v) => (v / 1e5).toFixed(3) },
	Pt: { label: "Pt", unit: "bar", format: (v) => (v / 1e5).toFixed(3) },
	T: { label: "T (static)", unit: "K", format: (v) => v.toFixed(2) },
	T_static: {
		label: "T_static (deprecated: T)",
		unit: "K",
		format: (v) => v.toFixed(2),
	},
	Tt: { label: "Tt", unit: "K", format: (v) => v.toFixed(2) },
	// A heat source's set outlet temperature (#792): a TOTAL temperature.
	Tt_target: {
		label: "Tt target (total)",
		unit: "K",
		format: (v) => v.toFixed(2),
	},
	rho: { label: "rho", unit: "kg/m³", format: (v) => v.toFixed(4) },
	a: { label: "a", unit: "m/s", format: (v) => v.toFixed(2) },
	h: { label: "h", unit: "J/kg", format: (v) => v.toFixed(0) },
	s: { label: "s", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	cp: { label: "cp", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	cv: { label: "cv", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	gamma: { label: "gamma", unit: "—", format: (v) => v.toFixed(4) },
	mu: { label: "mu", unit: "Pa·s", format: (v) => v.toExponential(2) },
	k: { label: "k", unit: "W/m-K", format: (v) => v.toExponential(2) },
	Pr: { label: "Pr", unit: "—", format: (v) => v.toFixed(3) },
	velocity: { label: "v", unit: "m/s", format: (v) => v.toFixed(3) },
	mw: { label: "mw", unit: "g/mol", format: (v) => v.toFixed(2) },
	nu: { label: "nu", unit: "m²/s", format: (v) => v.toExponential(2) },
	u: { label: "u", unit: "J/kg", format: (v) => v.toFixed(0) },
	// Tee junction diagnostics
	m_dot_com: { label: "ṁ com", unit: "kg/s", format: (v) => v.toFixed(4) },
	m_dot_straight: { label: "ṁ str", unit: "kg/s", format: (v) => v.toFixed(4) },
	m_dot_branch: { label: "ṁ br", unit: "kg/s", format: (v) => v.toFixed(4) },
	mass_flow_ratio: {
		// Every tee closure reports the BRANCH share of the common flow.
		label: "ṁ br/ṁ com",
		cardLabel: "ṁbr/ṁcom",
		unit: "kg/kg",
		format: (v) => v.toFixed(3),
	},
	K_straight: { label: "K str", unit: "—", format: (v) => v.toFixed(3) },
	K_branch: { label: "K br", unit: "—", format: (v) => v.toFixed(3) },
	// MPCE-v2 joining-flow K aliases (Bassett K11/K12 convention)
	K11: { label: "K11", unit: "—", format: (v) => v.toFixed(3) },
	K12: { label: "K12", unit: "—", format: (v) => v.toFixed(3) },
	// Junction pressures, every junction element (#773): P static, Pt total.
	// The ejector reports P_jct only (its mixing plane has no single Pt).
	Pt_jct: { label: "Pt jct", unit: "bar", format: (v) => (v / 1e5).toFixed(3) },
	P_jct: { label: "P jct", unit: "bar", format: (v) => (v / 1e5).toFixed(3) },
	// Ejector (#773): Huang's omega is reported as the entrainment ratio.
	entrainment_ratio: {
		label: "ṁs/ṁp (entrainment)",
		cardLabel: "ṁs/ṁp",
		unit: "kg/kg",
		format: (v) => v.toFixed(4),
	},
	m_dot_primary: {
		label: "ṁ primary",
		cardLabel: "ṁp",
		unit: "kg/s",
		format: (v) => v.toFixed(4),
	},
	m_dot_secondary: {
		label: "ṁ secondary",
		cardLabel: "ṁs",
		unit: "kg/s",
		format: (v) => v.toFixed(4),
	},
	n_reversed_ports: {
		label: "reversed ports",
		unit: "—",
		format: (v) => v.toFixed(0),
	},
	// Element end states: P_in/P_out and T_in/T_out the end nodes' static
	// state (T_in/T_out were the total temperature until 0.7, #782).
	P_in: { label: "Pi", unit: "bar", format: (v) => (v / 1e5).toFixed(3) },
	P_out: { label: "Po", unit: "bar", format: (v) => (v / 1e5).toFixed(3) },
	T_in: { label: "Ti (static)", unit: "K", format: (v) => v.toFixed(2) },
	T_out: { label: "To (static)", unit: "K", format: (v) => v.toFixed(2) },
	// Element-specific diagnostics (Total)
	Pt_in: { label: "Pti", unit: "bar", format: (v) => (v / 1e5).toFixed(3) },
	Pt_out: { label: "Pto", unit: "bar", format: (v) => (v / 1e5).toFixed(3) },
	Tt_in: { label: "Tti", unit: "K", format: (v) => v.toFixed(2) },
	Tt_out: { label: "Tto", unit: "K", format: (v) => v.toFixed(2) },
	// Thermal diagnostics
	T_hot: { label: "Th", unit: "K", format: (v) => v.toFixed(1) },
	T_aw: { label: "Taw", unit: "K", format: (v) => v.toFixed(1) },
};

// ---------------------------------------------------------------------------
// What the backend emits (derived, not hand-kept)
// ---------------------------------------------------------------------------
//
// resultKeys.generated.json is written by scripts/gui_result_keys.py: it
// solves a network per palette kind (and the variants that switch on optional
// diagnostics) and lists every finite numeric key the node's result carries,
// with its SI unit from gui/backend/units.py. The catalogue above only adds a
// display label, a display unit and a format; a key without an entry is still
// shown, under its own name, in its SI unit. quantities.test.ts fails when a
// catalogue entry is never emitted or its unit disagrees with the backend's.

type GeneratedResultKeys = {
	kinds: Record<string, string[]>;
	units: Record<string, string>;
};
const RESULT_KEYS = resultKeysJson as GeneratedResultKeys;

/** The pseudo-kind under which thermal-wall (edge) results are listed. */
export const THERMAL_WALL_KIND = "thermal_wall";

/** Every result key a solve can give an element or node of this kind, as the
 *  backend emits it (`state.<key>` for a node's state). Empty for a kind the
 *  backend does not solve (the probe itself). */
export function expectedResultKeys(kind: string | undefined): string[] {
	return (kind && RESULT_KEYS.kinds[kind]) || [];
}

/** Every kind the generated list covers (palette kinds and the wall). */
export function resultKinds(): string[] {
	return Object.keys(RESULT_KEYS.kinds);
}

/** The backend's SI unit of a result key ("-" dimensionless), if known. */
export function siUnit(key: string): string | undefined {
	return RESULT_KEYS.units[key];
}

/** A readable number for a key with no catalogue format. */
export function formatGeneric(v: number): string {
	if (v === 0) return "0";
	const a = Math.abs(v);
	if (a >= 1e6 || a < 1e-3) return v.toExponential(3);
	return String(Number.parseFloat(v.toPrecision(5)));
}

export type QuantityMeta = {
	label: string;
	/** Shorter label for node cards, where it differs from `label`. */
	cardLabel?: string;
	unit: string;
	format: (v: number) => string;
	/** False for a key with no catalogue entry (label = the key itself). */
	catalogued: boolean;
};

/** Label, display unit and format of a result key: the catalogue entry, or the
 *  key under its own name in the backend's SI unit. */
export function quantityMeta(key: string): QuantityMeta {
	const raw = key.startsWith("state.") ? key.slice(6) : key;
	const cat = QUANTITY_CATALOGUE[raw];
	if (cat) return { ...cat, catalogued: true };
	const si = siUnit(raw);
	return {
		label: raw,
		unit: si === undefined || si === "-" ? "—" : si,
		format: formatGeneric,
		catalogued: false,
	};
}
