/**
 * Quantity catalogue: mapping result keys to labels, units, and formatting logic.
 */
export const QUANTITY_CATALOGUE: Record<
	string,
	{ label: string; unit: string; format: (v: number) => string }
> = {
	// Scalar results (accessed via result.*)
	m_dot: { label: "ṁ", unit: "kg/s", format: (v) => v.toFixed(4) },
	mach: { label: "Ma", unit: "—", format: (v) => v.toFixed(4) },
	mach_in: { label: "Ma in", unit: "—", format: (v) => v.toFixed(4) },
	mach_out: { label: "Ma out", unit: "—", format: (v) => v.toFixed(4) },
	phi: { label: "Φ", unit: "—", format: (v) => v.toFixed(3) },
	theta: { label: "θ", unit: "—", format: (v) => v.toFixed(3) },
	p_ratio_total: { label: "P ratio", unit: "—", format: (v) => v.toFixed(4) },
	p_ratio: { label: "P ratio", unit: "—", format: (v) => v.toFixed(4) },
	Re: {
		label: "Re",
		unit: "—",
		format: (v) => v.toLocaleString("en", { maximumFractionDigits: 0 }),
	},
	Cd: { label: "Cd", unit: "—", format: (v) => v.toFixed(4) },
	// State results (accessed via result.state.*)
	P: { label: "P", unit: "Pa", format: (v) => v.toFixed(0) },
	P_total: { label: "Pt", unit: "Pa", format: (v) => v.toFixed(0) },
	T: { label: "T", unit: "K", format: (v) => v.toFixed(2) },
	T_total: { label: "Tt", unit: "K", format: (v) => v.toFixed(2) },
	rho: { label: "ρ", unit: "kg/m³", format: (v) => v.toFixed(4) },
	a: { label: "a", unit: "m/s", format: (v) => v.toFixed(2) },
	h: { label: "h", unit: "J/kg", format: (v) => v.toFixed(0) },
	s: { label: "s", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	cp: { label: "cp", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	cv: { label: "cv", unit: "J/kg·K", format: (v) => v.toFixed(2) },
	gamma: { label: "γ", unit: "—", format: (v) => v.toFixed(4) },
	mu: { label: "μ", unit: "Pa·s", format: (v) => v.toExponential(2) },
	k: { label: "k", unit: "W/m-K", format: (v) => v.toExponential(2) },
	Pr: { label: "Pr", unit: "—", format: (v) => v.toFixed(3) },
	velocity: { label: "v", unit: "m/s", format: (v) => v.toFixed(3) },
	mw: { label: "MW", unit: "g/mol", format: (v) => v.toFixed(2) },
	nu: { label: "ν", unit: "m²/s", format: (v) => v.toExponential(2) },
	u: { label: "u", unit: "J/kg", format: (v) => v.toFixed(0) },
};

export const BASIC_STATE_KEYS = [
	"P",
	"P_total",
	"T",
	"T_total",
	"m_dot",
	"rho",
	"mach",
	"velocity",
];
export const TRANSPORT_KEYS = [
	"mu",
	"k",
	"Pr",
	"Re",
	"cp",
	"cv",
	"gamma",
	"a",
	"mw",
	"nu",
	"h",
	"s",
	"u",
];
