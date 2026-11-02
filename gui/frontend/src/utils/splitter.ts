// The flow splitter's data (#903): one inlet to N outlets, outlet i at
// fractions[i] of the throughput. At most one fraction is null, the
// remainder (1 minus the others). `reference` is the index of the outlet
// whose valve is open (Pt_out = Pt_in); null means the remainder. The
// backend's FlowSplitter refuses a set that does not sum to 1, never
// renormalising; validateSplitter says so before a solve.

/** |sum - 1| allowed, the round-off of typed decimals:
 *  combaero.network.components.SPLITTER_FRACTION_SUM_TOL. */
export const SPLITTER_FRACTION_SUM_TOL = 1e-9;

/** Outlets a splitter may have: gui/backend/schemas.py SPLITTER_MAX_OUTLETS. */
export const SPLITTER_MAX_OUTLETS = 32;

/** Card width and the height per outlet row [px]. */
export const SPLITTER_WIDTH = 140;
export const SPLITTER_ROW_PX = 22;

export interface SplitterData {
	fractions?: (number | null)[];
	reference?: number | null;
	areas?: (number | null)[] | null;
}

export const DEFAULT_FRACTIONS: (number | null)[] = [0.5, null];

export function splitterFractions(data: SplitterData): (number | null)[] {
	return data.fractions && data.fractions.length >= 2
		? data.fractions
		: DEFAULT_FRACTIONS;
}

/** The card's height for n outlets: one row each, never below the 56 px of
 *  a two-port card. */
export function splitterHeight(n: number): number {
	return Math.max(56, n * SPLITTER_ROW_PX + 12);
}

/** Each outlet's fraction, the remainder resolved. */
export function resolvedFractions(fractions: (number | null)[]): number[] {
	const given = fractions.reduce<number>((a, f) => a + (f ?? 0), 0);
	return fractions.map((f) => (f == null ? 1 - given : f));
}

/** The reference outlet's index: the one chosen, else the remainder, else
 *  null (none: the solve is refused). */
export function splitterReference(data: SplitterData): number | null {
	const fractions = splitterFractions(data);
	if (data.reference != null) return data.reference;
	const rest = fractions.indexOf(null);
	return rest >= 0 ? rest : null;
}

/** Why the backend would refuse this splitter, one message per cause. */
export function validateSplitter(data: SplitterData): string[] {
	const errors: string[] = [];
	const fractions = splitterFractions(data);
	const n = fractions.length;
	if (n > SPLITTER_MAX_OUTLETS)
		errors.push(`at most ${SPLITTER_MAX_OUTLETS} outlets`);
	const rest = fractions.filter((f) => f == null).length;
	if (rest > 1) errors.push("at most one outlet can be the remainder");
	for (const [i, f] of fractions.entries()) {
		if (f != null && !(f > 0 && f < 1))
			errors.push(`outlet ${i + 1}: fraction must be in (0, 1)`);
	}
	const given = fractions.reduce<number>((a, f) => a + (f ?? 0), 0);
	if (rest === 1 && !(1 - given > SPLITTER_FRACTION_SUM_TOL))
		errors.push(
			`fractions sum to ${given.toPrecision(6)}: nothing left for the remainder`,
		);
	if (rest === 0 && !(Math.abs(given - 1) <= SPLITTER_FRACTION_SUM_TOL))
		errors.push(
			`fractions sum to ${given.toPrecision(6)}, not 1 (make one outlet the remainder)`,
		);
	const ref = splitterReference(data);
	if (ref == null)
		errors.push("choose the reference outlet (its valve open), or a remainder");
	else if (!(ref >= 0 && ref < n))
		errors.push(`reference outlet ${ref + 1} does not exist`);
	return errors;
}
