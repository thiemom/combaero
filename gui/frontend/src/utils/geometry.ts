import type { ResolvedField } from "./inheritance";

// One flow face's cross-section, shown the same way on every element (#559):
// diameter D, hydraulic diameter Dh, and area A, each either an input or
// computed. D and A are two views of one value (A = pi D^2 / 4; for a
// non-circular duct D is the area-equivalent diameter, which is what the
// solver's elements use), and Dh = D unless overridden.

export const areaFromDiameter = (d: number): number => (Math.PI * d * d) / 4;
export const diameterFromArea = (a: number): number =>
	Math.sqrt((4 * a) / Math.PI);

export interface GeometryInputs {
	// Which quantity the element stores ("D": a diameter field, "A": an area).
	stores: "D" | "A";
	// The stored (user) value of that quantity; null/undefined = inherit.
	userValue: number | null | undefined;
	// The solver's resolved value + provenance of that quantity (preview).
	inherited: ResolvedField | undefined;
	// The element's own Dh override, when it has one (undefined = no Dh field).
	userDh?: number | null;
	// The solver's resolved Dh, when the preview reports one.
	inheritedDh?: ResolvedField;
	// False for geometry that is the physics and never inherits (orifice bore).
	inheritable?: boolean;
}

export interface GeometryView {
	// Locked: the values are inherited and shown read-only until the user
	// opts out. Only possible when there is something to inherit.
	locked: boolean;
	// There is an inherited/defaulted value to fall back on (and so an
	// opt-out checkbox to offer).
	canInherit: boolean;
	// Effective values (user value, else inherited); null = unknown.
	D: number | null;
	A: number | null;
	Dh: number | null;
	// The user set D/A (bold), and/or Dh.
	userSet: boolean;
	userSetDh: boolean;
	// Provenance of the D/A value and of Dh, for the notes.
	source: string;
	dhSource: string;
}

const isNum = (v: number | null | undefined): v is number =>
	v !== null && v !== undefined && Number.isFinite(v);

export function geometryView(g: GeometryInputs): GeometryView {
	const userSet = isNum(g.userValue);
	const userSetDh = isNum(g.userDh);
	const inheritable = g.inheritable ?? true;
	const inheritedValue = isNum(g.inherited?.value) ? g.inherited.value : null;
	const canInherit = inheritable && inheritedValue !== null;
	// An override of either D/A or Dh opts the whole face out of inheritance;
	// an element with nothing to inherit is necessarily editable.
	const locked = canInherit && !userSet && !userSetDh;

	const stored = userSet ? (g.userValue as number) : inheritedValue;
	const D =
		stored === null
			? null
			: g.stores === "D"
				? stored
				: diameterFromArea(stored);
	const A =
		stored === null
			? null
			: g.stores === "A"
				? stored
				: areaFromDiameter(stored);

	let Dh: number | null = D;
	let dhSource = "circular";
	if (userSetDh) {
		Dh = g.userDh as number;
		dhSource = "user";
	} else if (isNum(g.inheritedDh?.value)) {
		// An element that inherits Dh on its own (a chamber takes it from the
		// upstream channel even when its area is set) reports the solver's Dh;
		// one that passes no inherited Dh is circular on D.
		Dh = g.inheritedDh.value;
		dhSource = g.inheritedDh.source || "circular";
	}

	return {
		locked,
		canInherit,
		D,
		A,
		Dh,
		userSet,
		userSetDh,
		source: userSet ? "user" : (g.inherited?.source ?? ""),
		dhSource,
	};
}

/**
 * The stored value for a D or A entered by the user, converted to the
 * quantity the element stores.
 */
export function storedFrom(
	stores: "D" | "A",
	entered: "D" | "A",
	value: number,
): number {
	if (stores === entered) return value;
	return entered === "D" ? areaFromDiameter(value) : diameterFromArea(value);
}
